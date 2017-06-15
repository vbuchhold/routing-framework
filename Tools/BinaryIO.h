#pragma once

#include <cassert>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

#include <boost/dynamic_bitset.hpp>

namespace bio {

// Returns the number of bytes the specified self-contained object occupies on disk.
template <typename T>
inline int size(const T& /*obj*/) {
  return sizeof(T);
}

// Returns the number of bytes the specified string occupies on disk.
inline int size(const std::string& str) {
  return str.size() + 1;
}

// Returns the number of bytes the specified C-style string occupies on disk.
inline int size(const char* const str) {
  return std::strlen(str) + 1;
}

// Returns the number of bytes the specified vector of self-contained objects occupies on disk.
template <typename T, typename AllocT>
inline int size(const std::vector<T, AllocT>& vec) {
  return sizeof(int) + vec.size() * sizeof(T);
}

// Returns the number of bytes the specified bit-vector occupies on disk.
template <typename BlockT, typename AllocT>
inline int size(const boost::dynamic_bitset<BlockT, AllocT>& vec) {
  return 2 * sizeof(int) + vec.num_blocks() * sizeof(BlockT);
}

// Reads a self-contained object from a binary file.
template <typename T>
inline void read(std::ifstream& in, T& obj) {
  in.read(reinterpret_cast<char*>(&obj), sizeof(T));
  assert(in.good());
}

// Reads a string from a binary file.
inline void read(std::ifstream& in, std::string& str) {
  getline(in, str, '\0');
  assert(in.good());
}

// Reads a vector of self-contained objects from a binary file.
template <typename T, typename AllocT>
inline void read(std::ifstream& in, std::vector<T, AllocT>& vec) {
  int size;
  read(in, size);
  vec.resize(size);
  in.read(reinterpret_cast<char*>(vec.data()), size * sizeof(T));
  assert(in.good());
}

// Reads a bit-vector from a binary file.
template <typename BlockT, typename AllocT>
inline void read(std::ifstream& in, boost::dynamic_bitset<BlockT, AllocT>& vec) {
  int size;
  read(in, size);
  vec.resize(size);
  std::vector<BlockT> blocks;
  read(in, blocks);
  boost::from_block_range(blocks.begin(), blocks.end(), vec);
}

// Writes a self-contained object to a binary file.
template <typename T>
inline void write(std::ofstream& out, const T& obj) {
  out.write(reinterpret_cast<const char*>(&obj), sizeof(T));
  assert(out.good());
}

// Writes a string to a binary file.
inline void write(std::ofstream& out, const std::string& str) {
  out.write(str.data(), str.size() + 1);
  assert(out.good());
}

// Writes a C-style string to a binary file.
inline void write(std::ofstream& out, const char* const str) {
  out.write(str, std::strlen(str) + 1);
  assert(out.good());
}

// Writes a vector of self-contained objects to a binary file.
template <typename T, typename AllocT>
inline void write(std::ofstream& out, const std::vector<T, AllocT>& vec) {
  const int size = vec.size();
  write(out, size);
  out.write(reinterpret_cast<const char*>(vec.data()), size * sizeof(T));
  assert(out.good());
}

// Writes a bit-vector to a binary file.
template <typename BlockT, typename AllocT>
inline void write(std::ofstream& out, const boost::dynamic_bitset<BlockT, AllocT>& vec) {
  const int size = vec.size();
  write(out, size);
  std::vector<BlockT> blocks(vec.num_blocks());
  boost::to_block_range(vec, blocks.begin());
  write(out, blocks);
}

}
