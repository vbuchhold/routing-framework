#pragma once

#include <algorithm>
#include <cassert>
#include <fstream>
#include <random>
#include <utility>
#include <vector>

#include <boost/dynamic_bitset.hpp>

#include "Tools/BinaryIO.h"

// A permutation is a rearrangement of members of a sequence into a new sequence.
class Permutation {
 public:
  // Constructs an empty permutation.
  Permutation() = default;

  // Constructs an uninitialized permutation of the specified size.
  explicit Permutation(const int size) : permutation(size) {}

  // Constructs a permutation backed by the specified vector (which must represent a permutation).
  explicit Permutation(const std::vector<int>& perm) : permutation(perm) { assert(validate()); }
  explicit Permutation(std::vector<int>&& perm) : permutation(std::move(perm)) {
    assert(validate());
  }

  // Returns a random permutation of the specified size.
  static Permutation getRandomPermutation(const int size, std::default_random_engine& eng) {
    Permutation perm(size);
    for (int i = 0; i != size; ++i)
      perm[i] = i;
    std::shuffle(perm.permutation.begin(), perm.permutation.end(), eng);
    return perm;
  }

  // Returns the size of the permutation.
  int size() const {
    return permutation.size();
  }

  // Returns the location of the element at oldPos in the new sequence.
  int operator[](const int oldPos) const {
    assert(oldPos >= 0); assert(oldPos < permutation.size());
    return permutation[oldPos];
  }

  // Returns the location of the element at oldPos in the new sequence.
  int& operator[](const int oldPos) {
    assert(oldPos >= 0); assert(oldPos < permutation.size());
    return permutation[oldPos];
  }

  // Inverts the permutation.
  void invert() {
    Permutation temp = getInversePermutation();
    *this = std::move(temp);
  }

  // Returns the inverse permutation that maps new locations to old locations.
  Permutation getInversePermutation() const {
    assert(validate());
    Permutation perm(permutation.size());
    for (int oldPos = 0; oldPos != permutation.size(); ++oldPos)
      perm[permutation[oldPos]] = oldPos;
    return perm;
  }

  // Reorders the elements in the specified container according to the permutation.
  template <typename ContT>
  void applyTo(ContT& cont) const {
    assert(cont.size() == permutation.size());
    ContT temp(cont.size());
    for (int i = 0; i != cont.size(); ++i)
      temp[permutation[i]] = std::move(cont[i]);
    cont.swap(temp);
  }

  // Checks if the permutation indeed represents a valid permutation.
  bool validate() const {
    const int size = permutation.size();
    boost::dynamic_bitset<> isContained(size);
    for (int i = 0; i != size; ++i) {
      const int elem = permutation[i];
      assert(elem >= 0); assert(elem < size);
      assert(!isContained[elem]);
      isContained[elem] = true;
    }
    return true;
  }

  // Reads a permutation from a binary file.
  void readFrom(std::ifstream& in) {
    bio::read(in, permutation);
  }

  // Write a permutation to a binary file.
  void writeTo(std::ofstream& out) const {
    bio::write(out, permutation);
  }

 private:
  std::vector<int> permutation; // permutation[i] represents the new location of the element at i.
};
