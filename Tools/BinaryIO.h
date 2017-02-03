#pragma once

#include <cassert>
#include <fstream>
#include <vector>

#include <boost/dynamic_bitset.hpp>

// Reads a self-contained object from a binary file.
template <typename T>
inline void read(std::ifstream& in, T& obj) {
  in.read(reinterpret_cast<char*>(&obj), sizeof(T));
  assert(in.good());
}

// Reads a vector of self-contained objects from a binary file.
template <typename T, typename AllocT>
inline void read(std::ifstream& in, std::vector<T, AllocT>& vec) {
  in.read(reinterpret_cast<char*>(vec.data()), vec.size() * sizeof(T));
  assert(in.good());
}

// Reads a bit-vector from a binary file.
template <typename BlockT, typename AllocT>
inline void read(std::ifstream& in, boost::dynamic_bitset<BlockT, AllocT>& vec) {
  std::vector<BlockT> blocks(vec.num_blocks());
  read(in, blocks);
  boost::from_block_range(blocks.begin(), blocks.end(), vec);
}

// Writes a self-contained object to a binary file.
template <typename T>
inline void write(std::ofstream& out, const T& obj) {
  out.write(reinterpret_cast<const char*>(&obj), sizeof(T));
  assert(out.good());
}

// Writes a vector of self-contained objects to a binary file.
template <typename T, typename AllocT>
inline void write(std::ofstream& out, const std::vector<T, AllocT>& vec) {
  out.write(reinterpret_cast<const char*>(vec.data()), vec.size() * sizeof(T));
  assert(out.good());
}

// Writes a bit-vector to a binary file.
template <typename BlockT, typename AllocT>
inline void write(std::ofstream& out, const boost::dynamic_bitset<BlockT, AllocT>& vec) {
  std::vector<BlockT> blocks(vec.num_blocks());
  boost::to_block_range(vec, blocks.begin());
  write(out, blocks);
}
