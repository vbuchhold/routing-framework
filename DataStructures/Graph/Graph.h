#pragma once

#include <algorithm>
#include <cassert>
#include <ostream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "DataStructures/Graph/Export/DefaultExporter.h"
#include "DataStructures/Graph/Import/XatfImporter.h"
#include "DataStructures/Utilities/Permutation.h"
#include "Tools/TemplateProgramming.h"
#include "Tools/Workarounds.h"

namespace impl {

// This struct delimits the range of edges out of a particular vertex v. It always stores the
// index of the first edge out of v. If the graph is static, the index of the last edge out of v
// is given implicitly by the index of the first edge out of (v + 1). If the graph is dynamic, we
// store the index of the last edge out of v explicitly.

struct StaticOutEdgeRange {
  int first() const { return firstEdge; }
  int last() const { return firstEdge; }

  int& first() { return firstEdge; }
  int& last() { return firstEdge; }

  int firstEdge;
};

struct DynamicOutEdgeRange {
  int first() const { return firstEdge; }
  int last() const { return lastEdge; }

  int& first() { return firstEdge; }
  int& last() { return lastEdge; }

  int firstEdge;
  int lastEdge;
};

}

// Auxiliary classes for grouping vertex and edge attributes, respectively.
template <typename ...Attributes>
class VertexAttrs {};
template <typename ...Attributes>
class EdgeAttrs {};

template <typename VertexAttributes, typename EdgeAttributes, bool dynamic>
class Graph;

// Type trait checking if the specified type is a graph.
template <typename T>
struct IsGraph : std::false_type {};
template <typename T>
struct IsGraph<const T> : IsGraph<T> {};
template <typename VertexAttributes, typename EdgeAttributes, bool dynamic>
struct IsGraph<Graph<VertexAttributes, EdgeAttributes, dynamic>> : std::true_type {};

// Graph implementation using adjacency arrays. Can be instantiated as a static or dynamic data
// structure, where the dynamic variant allows holes in the edge array in order to efficiently
// support inserting and removing edges. Vertex and edge attributes can be added as needed.
// Attributes are stored in a separate array each.
template <typename ...VertexAttributes, typename ...EdgeAttributes, bool dynamic>
class Graph<VertexAttrs<VertexAttributes...>, EdgeAttrs<EdgeAttributes...>, dynamic>
    : public VertexAttributes..., public EdgeAttributes... {
  // All specializations of this class template should be friends of each other.
  template <typename, typename, bool>
  friend class Graph;

 private:
  // The type of this template specialization.
  using GraphT = Graph<VertexAttrs<VertexAttributes...>, EdgeAttrs<EdgeAttributes...>, dynamic>;

 public:
  // Constructs an empty graph.
  Graph() {
    init();
  }

  // Constructs an empty graph. The constructed graph can hold at least the specified number of
  // vertices and edges without requiring reallocation.
  Graph(const int numVertices, const int numEdges) {
    reserve(numVertices, numEdges);
    init();
  }

  // Constructs a graph from a file. Different importers support different file formats.
  template <typename ImporterT = XatfImporter>
  explicit Graph(const std::string& filename, ImporterT im = ImporterT()) : edgeCount(0) {
    importFrom(filename, im);
  }

  // Converts an arbitrary source graph into in arbitrary destination graph. Attributes associated
  // with both the source and the destination graph are copied or, if possible, moved. Attributes
  // associated only with the destination graph are defaulted.
  // CAUTION: IF THE DESTINATION GRAPH IS STATIC, THE SOURCE GRAPH HAS TO BE DEFRAGMENTED.
  template <
      typename SourceT,
      typename = std::enable_if_t<IsGraph<std::remove_reference_t<SourceT>>::value>>
  explicit Graph(SourceT&& src) {
    // Assert that if the destination graph is static, the source graph has to be defragmented.
    assert(dynamic || src.isDefrag());
    setOutEdges(std::forward<SourceT>(src).outEdges);
    edgeHeads = std::forward<SourceT>(src).edgeHeads;
    edgeCount = src.edgeCount;
    RUN_FORALL(setAttribute<VertexAttributes>(std::forward<SourceT>(src), numVertices()));
    RUN_FORALL(setAttribute<EdgeAttributes>(std::forward<SourceT>(src), edgeHeads.size()));
  }

  // Returns a vector with the names of all attributes of the graph.
  static std::vector<std::string> getAttributeNames() {
    std::vector<std::string> names;
    names.reserve(sizeof...(VertexAttributes) + sizeof...(EdgeAttributes));
    RUN_FORALL(names.push_back(VertexAttributes::NAME));
    RUN_FORALL(names.push_back(EdgeAttributes::NAME));
    return names;
  }

  // Returns the number of vertices in the graph.
  int numVertices() const {
    return outEdges.size() - !dynamic;
  }

  // Returns the number of edges in the graph.
  int numEdges() const {
    return edgeCount;
  }

  // Returns the index of the edge with the largest index.
  int maxEdgeIndex() const {
    return edgeHeads.size() - 1;
  }

  // Returns the index of the first edge out of vertex v.
  int firstEdge(const int v) const {
    assert(v >= 0); assert(v < numVertices());
    return outEdges[v].first();
  }

  // Returns the index of the past-the-last edge out of vertex v.
  int lastEdge(const int v) const {
    assert(v >= 0); assert(v < numVertices());
    return outEdges[v + !dynamic].last();
  }

  // Returns true if edge e is an actual edge. If the graph is dynamic, edges may also be holes.
  bool isValidEdge(const int e) const {
    assert(e >= 0); assert(e < edgeHeads.size());
    return edgeHeads[e] != INVALID_EDGE;
  }

  // Returns the head vertex of edge e.
  int edgeHead(const int e) const {
    assert(e >= 0); assert(e < edgeHeads.size()); assert(isValidEdge(e));
    return edgeHeads[e];
  }

  // Returns true if the graph contains an edge from tail to head.
  bool containsEdge(const int tail, const int head) const {
    for (int e = firstEdge(tail); e != lastEdge(tail); ++e)
      if (edgeHead(e) == head)
        return true;
    return false;
  }

  // Ensures that the graph can hold at least the specified number of vertices and edges without
  // requiring reallocation.
  void reserve(const int numVertices, const int numEdges) {
    outEdges.reserve(numVertices + !dynamic);
    RUN_FORALL(VertexAttributes::reserve(numVertices));

    edgeHeads.reserve(numEdges);
    RUN_FORALL(EdgeAttributes::reserve(numEdges));
  }

  // Appends a vertex to the graph. Returns the ID of the new vertex.
  int appendVertex() {
    OutEdgeRange range;
    range.first() = edgeHeads.size();
    range.last() = range.first();
    outEdges.push_back(range);
    RUN_FORALL(VertexAttributes::values.push_back(use(VertexAttributes::DEFAULT_VALUE)));
    return numVertices() - 1;
  }

  // Appends a vertex with the specified attributes to the graph. Returns the ID of the new vertex.
  template <typename ...Attrs>
  int appendVertex(Attrs&& ...attrs) {
    appendVertex();
    RUN_FORALL(VertexAttributes::values.back() = std::forward<Attrs>(attrs));
    return numVertices() - 1;
  }

  // Appends the specified number of vertices to the graph.
  void appendVertices(const int num) {
    assert(num >= 0);
    OutEdgeRange range;
    range.first() = edgeHeads.size();
    range.last() = range.first();
    outEdges.insert(outEdges.end(), num, range);
    const int size = numVertices();
    RUN_FORALL(VertexAttributes::values.resize(size, use(VertexAttributes::DEFAULT_VALUE)));
  }

  // Inserts an edge from the last inserted vertex to v. Returns the ID of the inserted edge. Note
  // that v does not need to be already present in static graphs.
  int appendEdge(const std::enable_if_t<!dynamic, int> v) {
    assert(v >= 0);
    assert(numVertices() > 0);
    ++outEdges.back().last();
    edgeHeads.push_back(v);
    RUN_FORALL(EdgeAttributes::values.push_back(use(EdgeAttributes::DEFAULT_VALUE)));
    ++edgeCount;
    return numEdges() - 1;
  }

  // Inserts an edge from the last inserted vertex to v. Returns the ID of the inserted edge. Note
  // that v does not need to be already present in static graphs.
  int appendEdge(const std::enable_if_t<dynamic, int> v) {
    assert(v >= 0);
    assert(numVertices() > 0);
    ++outEdges.back().last();
    edgeHeads.push_back(v);
    RUN_FORALL(EdgeAttributes::values.push_back(use(EdgeAttributes::DEFAULT_VALUE)));
    ++edgeCount;
    return numEdges() - 1;
    return insertEdge(numVertices() - 1, v);
  }

  // Inserts an edge with the specified attributes from the last inserted vertex to v. Returns the
  // ID of the inserted edge. Note that v does not need to be already present in static graphs.
  template <typename ...Attrs>
  int appendEdge(const int v, Attrs&& ...attrs) {
    const int idx = appendEdge(v);
    RUN_FORALL(EdgeAttributes::values[idx] = std::forward<Attrs>(attrs));
    return idx;
  }

  // Inserts an edge from u to v. Returns the ID of the newly inserted edge. Note that this
  // operation is not supported by static graphs.
  std::enable_if_t<dynamic, int> insertEdge(const int /*u*/, const int /*v*/) { return 0; }

  // Inserts an edge with the specified attributes from u to v. Returns the ID of the newly
  // inserted edge. Note that this operation is not supported by static graphs.
  template <typename ...Attrs>
  std::enable_if_t<dynamic, int> insertEdge(const int u, const int v, Attrs&& ...attrs) {
    const int idx = insertEdge(u, v);
    RUN_FORALL(EdgeAttributes::values[idx] = std::forward<Attrs>(attrs));
    return idx;
  }

  // Sets the head of edge e to vertex v.
  void setEdgeHead(const int e, const int v) {
    assert(e >= 0); assert(e <= maxEdgeIndex());
    assert(v >= 0); assert(v < numVertices());
    edgeHeads[e] = v;
  }

  // Removes all vertices and edges from the graph.
  void clear() {
    outEdges.clear();
    RUN_FORALL(VertexAttributes::values.clear());

    edgeHeads.clear();
    RUN_FORALL(EdgeAttributes::values.clear());
    edgeCount = 0;

    // If the graph is static, add a sentinel to the vertex array.
    if (!dynamic) {
      OutEdgeRange range;
      range.first() = 0;
      outEdges.push_back(range);
    }
  }

  // Returns true if the edge arrays are sorted by tail ID and contain no holes, false otherwise.
  bool isDefrag() const {
    if (!dynamic || numVertices() == 0)
      return true;
    if (outEdges[0].first() != 0)
      return false;
    for (int v = 0; v < numVertices() - 1; ++v)
      if (outEdges[v].last() != outEdges[v + 1].first())
        return false;
    if (outEdges.back().last() != edgeHeads.size())
      return false;
    return true;
  }

  // Sorts the edge arrays by tail ID, leaving them without holes.
  void defrag() {
    if (!dynamic)
      return;
    Permutation perm(edgeHeads.size());
    int newEdgeIdx = 0;

    // Push all valid edges to the front, breaking ties by tail ID.
    for (int u = 0; u < numVertices(); ++u) {
      const int first = newEdgeIdx; // The index of the first edge out of u.
      for (int e = firstEdge(u); e < lastEdge(u); ++e)
        perm[e] = newEdgeIdx++;
      outEdges[u].first() = first;
      outEdges[u].last() = newEdgeIdx;
    }

    // Push all invalid edges to the back.
    for (int e = 0; e < edgeHeads.size(); ++e)
      if (edgeHeads[e] == INVALID_EDGE)
        perm[e] = newEdgeIdx++;
    assert(newEdgeIdx == edgeHeads.size());

    // Reorder the edge arrays according to the permutation.
    perm.applyTo(edgeHeads);
    edgeHeads.resize(numEdges());
    RUN_FORALL(perm.applyTo(EdgeAttributes::values));
    RUN_FORALL(EdgeAttributes::values.resize(numEdges()));
  }

  // Reorders the vertices according to the specified permutation.
  void permuteVertices(const Permutation& perm) {
    assert(perm.size() == numVertices()); assert(perm.validate());
    if (dynamic) {
      perm.applyTo(outEdges);
    } else {
      std::vector<OutEdgeRange> temp(outEdges.size());
      temp.back().first() = numEdges();
      Permutation inversePerm = perm.getInversePermutation();
      Permutation edgePerm(numEdges());
      int newEdgeIdx = 0;

      // Sort the edge arrays by new tail ID.
      for (int newTailIdx = 0; newTailIdx < numVertices(); ++newTailIdx) {
        const int oldTailIdx = inversePerm[newTailIdx];
        temp[newTailIdx].first() = newEdgeIdx;
        for (int e = firstEdge(oldTailIdx); e < lastEdge(oldTailIdx); ++e)
          edgePerm[e] = newEdgeIdx++;
      }

      edgePerm.applyTo(edgeHeads);
      RUN_FORALL(edgePerm.applyTo(EdgeAttributes::values));
      outEdges.swap(temp);
    }

    RUN_FORALL(perm.applyTo(VertexAttributes::values));

    // Update edge heads.
    for (auto& head : edgeHeads)
      head = perm[head];
  }

  // Removes all vertices and edges that do not lie in the vertex-induced subgraph specified by the
  // given bitmask. The bitmask must contain one bit per vertex, which should be set iff the vertex
  // belongs to the subgraph.
  void extractVertexInducedSubgraph(const boost::dynamic_bitset<> bitmask) {
    GraphT temp = getVertexInducedSubgraph(bitmask);
    *this = std::move(temp);
  }

  // Returns the vertex-induced subgraph specified by the given bitmask. The bitmask must contain
  // one bit per vertex, which should be set iff the corresponding vertex belongs to the subgraph.
  GraphT getVertexInducedSubgraph(const boost::dynamic_bitset<>& bitmask) const {
    // Assign new sequential IDs to the vertices in the subgraph.
    assert(bitmask.size() == numVertices());
    int nextId = 0;
    std::vector<int> origToNewIds(numVertices(), -1);
    for (int v = bitmask.find_first(); v != boost::dynamic_bitset<>::npos; v = bitmask.find_next(v))
      origToNewIds[v] = nextId++;

    GraphT subgraph;
    subgraph.outEdges.resize(nextId + !dynamic);
    RUN_FORALL(subgraph.VertexAttributes::values.resize(nextId));

    for (int i = 0, u = bitmask.find_first(); i != nextId; ++i, u = bitmask.find_next(u)) {
      // Copy the current vertex belonging to the subgraph.
      subgraph.outEdges[i].first() = subgraph.edgeHeads.size();
      RUN_FORALL(subgraph.VertexAttributes::values[i] = VertexAttributes::values[u]);

      // Copy the edges out of u going to vertices belonging to the subgraph.
      for (int e = firstEdge(u); e != lastEdge(u); ++e) {
        const int v = origToNewIds[edgeHeads[e]];
        if (v != -1) {
          subgraph.edgeHeads.push_back(v);
          RUN_FORALL(subgraph.EdgeAttributes::values.push_back(EdgeAttributes::values[e]));
        }
      }

      if (dynamic)
        subgraph.outEdges[i].last() = subgraph.edgeHeads.size();
    }

    subgraph.outEdges.back().last() = subgraph.edgeHeads.size();
    subgraph.edgeCount = subgraph.edgeHeads.size();
    return subgraph;
  }

  // Reads a graph from disk. Different importers support different file formats.
  template <typename ImporterT = XatfImporter>
  void importFrom(const std::string& filename, ImporterT im = ImporterT()) {
    clear();
    outEdges.clear();

    // Open the input file(s), read the header line(s), and allocate the vertex and edge arrays.
    im.init(filename);
    assert(im.numVertices() >= 0);
    assert(im.numEdges() >= 0);
    reserve(im.numVertices(), im.numEdges());

    // Read the vertices, one after another. The vertices do not have to be in any particular order.
    int vertexCount = 0;
    while (im.nextVertex()) {
      assert(im.vertexId() >= 0);
      assert(im.numVertices() == 0 || im.vertexId() < im.numVertices());

      // Ensure that the vertex arrays can accommodate the vertex we just read.
      if (im.vertexId() >= numVertices()) {
        OutEdgeRange range;
        range.first() = -1;
        outEdges.resize(im.vertexId() + 1 + !dynamic, range);
        RUN_FORALL(VertexAttributes::values.resize(im.vertexId() + 1));
      }

      // Store the vertex and its attributes at the correct position in the vertex arrays.
      assert(outEdges[im.vertexId()].first() == -1);
      outEdges[im.vertexId()].first() = 0;
      RUN_FORALL(VertexAttributes::values[im.vertexId()] =
          im.template getValue<VertexAttributes>());
      ++vertexCount;
    }
    assert(numVertices() == vertexCount);
    assert(im.numVertices() == 0 || numVertices() == im.numVertices());

    // Read the edges, one after another. The edges do not have to be in any particular order.
    std::vector<int> edgeTails;
    edgeTails.reserve(im.numEdges());
    bool edgesSorted = true; // Indicates if the edges are already sorted by tail ID.
    int prevTailId = -1; // The tail ID of the previous edge.
    while (im.nextEdge()) {
      assert(im.edgeTail() >= 0); assert(im.edgeTail() < numVertices());
      assert(im.edgeHead() >= 0); assert(im.edgeHead() < numVertices());
      edgeTails.push_back(im.edgeTail());
      edgeHeads.push_back(im.edgeHead());
      RUN_FORALL(EdgeAttributes::values.push_back(im.template getValue<EdgeAttributes>()));
      outEdges[im.edgeTail()].first()++;
      edgesSorted &= prevTailId <= im.edgeTail();
      prevTailId = im.edgeTail();
    }
    assert(im.numEdges() == 0 || edgeTails.size() == im.numEdges());
    edgeCount = edgeTails.size();

    im.close();

    // Before this loop, outEdges[v].first() stores the outdegree of v. After the loop,
    // outEdges[v].first() stores the index of the first edge out of v.
    int firstEdge = 0; // The index of the first edge out of the current/next vertex.
    std::swap(outEdges[0].first(), firstEdge);
    for (int v = 1; v != numVertices(); ++v) {
      std::swap(outEdges[v].first(), firstEdge);
      outEdges[v - dynamic].last() = outEdges[v].first();
      firstEdge += outEdges[v].first();
    }
    outEdges.back().last() = numEdges();

    // Sort the edges by tail ID if they are not already sorted.
    if (!edgesSorted) {
      // Compute a permutation mapping each edge to its correct position.
      Permutation perm(numEdges());
      for (int e = 0; e != numEdges(); ++e)
        perm[e] = outEdges[edgeTails[e]].first()++;
      for (int v = numVertices() - 1; v != 0; --v)
        outEdges[v].first() = outEdges[v - 1].first();
      outEdges[0].first() = 0;

      // Apply the permutation to each edge array.
      permuteEdges(perm);
    }

    assert(validate());
  }

  // Writes a graph to disk. Different exporters support different file formats.
  template <typename ExporterT = DefaultExporter>
  void exportTo(const std::string& /*filename*/, ExporterT /*ex*/ = ExporterT()) const {}

  // Checks if the graph is consistent.
  bool validate() const {
    boost::dynamic_bitset<> visited(edgeHeads.size());
    int numEdges = 0;
    for (int u = 0; u != numVertices(); ++u) {
      // Assert that a valid range of edges is associated with each vertex.
      assert(firstEdge(u) >= 0);
      assert(lastEdge(u) <= edgeHeads.size());
      assert(firstEdge(u) <= lastEdge(u));

      // Assert that all edges incident on u point to a valid vertex, and that the edge ranges do
      // not overlap.
      for (int e = firstEdge(u); e != lastEdge(u); ++e) {
        assert(edgeHeads[e] >= 0); assert(edgeHeads[e] < numVertices());
        assert(!visited[e]);
        visited[e] = true;
        ++numEdges;
      }
    }

    // Assert that all unvisited edges are invalid.
    visited.flip();
    for (int e = visited.find_first(); e != boost::dynamic_bitset<>::npos; e = visited.find_next(e))
      assert(edgeHeads[e] == INVALID_EDGE);

    assert(edgeCount == numEdges);
    return true;
  }

 private:
  using OutEdgeRange =
      std::conditional_t<dynamic, impl::DynamicOutEdgeRange, impl::StaticOutEdgeRange>;

  // If a graph is dynamic, its edge arrays may contain "holes". To indicate that an edge is a hole
  // and not an actual edge, we store this value as its head.
  static constexpr int INVALID_EDGE = -1;

  // Initializes an empty graph.
  void init() {
    edgeCount = 0;

    // If the graph is static, add a sentinel to the vertex array.
    if (!dynamic) {
      OutEdgeRange range;
      range.first() = 0;
      outEdges.push_back(range);
    }
  }

  // Auxiliary functions for the converting constructor. Set the vector of out-edge ranges
  // associated with the destination graph. If the source and the destination graph are both static
  // or are both dynamic, the vector of out-edge ranges is simply copied or, if possible, moved. If
  // they are not, there is an overload for the special case where a static graph should be
  // converted to a dynamic graph or vice versa.

  template <
      typename SourceT,
      typename = std::enable_if_t<
          std::is_same<OutEdgeRange, typename std::remove_reference_t<SourceT>::value_type>::value>>
  void setOutEdges(SourceT&& srcOutEdges) {
    outEdges = std::forward<SourceT>(srcOutEdges);
  }

  template <
      typename SourceOutEdgeRangeT,
      typename = std::enable_if_t<!std::is_same<OutEdgeRange, SourceOutEdgeRangeT>::value>>
  void setOutEdges(const std::vector<SourceOutEdgeRangeT>& srcOutEdges) {
    outEdges.resize(srcOutEdges.size() + !dynamic - dynamic);
    if (numVertices() != 0) {
      outEdges.back().last() = srcOutEdges.back().last();
      for (int v = numVertices() - 1; v != 0; --v) {
        outEdges[v].first() = srcOutEdges[v].first();
        outEdges[v - 1].last() = srcOutEdges[v].first();
      }
      outEdges[0].first() = 0;
    } else if (!dynamic) {
      outEdges[0].first() = 0;
    }
  }

  // Auxiliary functions for the converting constructor. Set the specified vertex or edge attribute
  // Attr associated with the destination graph. If Attr is also associated with the source graph,
  // the attribute's values are copied or, if possible, moved. If not, the attribute's values are
  // set to the attribute's default value. The parameter size has to be the size of the vertex/edge
  // attribute vectors.

  template <
      typename Attr, typename SourceT,
      typename = std::enable_if_t<std::is_base_of<Attr, std::remove_reference_t<SourceT>>::value>>
  void setAttribute(SourceT&& src, const int /*size*/) {
    Attr::values = std::forward<SourceT>(src).Attr::values;
  }

  template <
      typename Attr, typename SourceT,
      typename = std::enable_if_t<!std::is_base_of<Attr, SourceT>::value>>
  void setAttribute(const SourceT& /*src*/, const int size) {
    Attr::values.resize(size, use(Attr::DEFAULT_VALUE));
  }

  // Reorders the edges according to the specified permutation.
  // CAUTION: IT IS THE RESPONSIBILITY OF THE USER TO GUARANTEE THAT CONSISTENCY IS MAINTAINED.
  void permuteEdges(const Permutation& perm) {
    assert(perm.validate());
    perm.applyTo(edgeHeads);
    RUN_FORALL(perm.applyTo(EdgeAttributes::values));
  }

  std::vector<OutEdgeRange> outEdges; // The ranges of outgoing edges of the vertices.
  std::vector<int> edgeHeads;         // The head vertices of the edges.

  int edgeCount; // The number of edges in the graph.
};

// Write a textual representation to the specified output stream.
template <typename VertexAttrs, typename EdgeAttrs, bool dynamic>
inline std::ostream& operator<<(std::ostream& os, const Graph<VertexAttrs, EdgeAttrs, dynamic>& g) {
  os << "#Vertices=" << g.numVertices() << " #Edges=" << g.numEdges() << std::endl;
  for (int u = 0; u != g.numVertices(); ++u) {
    os << u << ":";
    for (int e = g.firstEdge(u); e != g.lastEdge(u); ++e)
      os << " " << g.edgeHead(e);
    os << std::endl;
  }
  return os;
}

// Alias templates for static and dynamic graphs.
template <typename VertexAttributes = VertexAttrs<>, typename EdgeAttributes = EdgeAttrs<>>
using StaticGraph = Graph<VertexAttributes, EdgeAttributes, false>;
template <typename VertexAttributes = VertexAttrs<>, typename EdgeAttributes = EdgeAttrs<>>
using DynamicGraph = Graph<VertexAttributes, EdgeAttributes, true>;

// Iteration macros for conveniently looping through vertices or edges of a graph.
#define FORALL_VERTICES(G, u) for (int u = 0; u != G.numVertices(); ++u)
#define FORALL_EDGES(G, e) for (int e = 0; e <= G.maxEdgeIndex(); ++e)
#define FORALL_VALID_EDGES(G, u, e) for (int u = 0; u != G.numVertices(); ++u) \
    for (int e = G.firstEdge(u); e != G.lastEdge(u); ++e)
#define FORALL_INCIDENT_EDGES(G, u, e) for (int e = G.firstEdge(u); e != G.lastEdge(u); ++e)
