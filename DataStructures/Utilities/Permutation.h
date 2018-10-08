#pragma once

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <numeric>
#include <vector>

#include "DataStructures/Containers/BitVector.h"
#include "Tools/BinaryIO.h"

// A permutation is a rearrangement of members of a sequence into a new sequence.
class Permutation {
 public:
  // A const iterator type referring to a new location.
  using ConstIterator = std::vector<int32_t>::const_iterator;

  // Constructs an empty permutation.
  Permutation() = default;

  // Constructs an uninitialized permutation of the specified size.
  explicit Permutation(const int size) : permutation(size) {}

  // Constructs a permutation equal to the one specified as an iterator range.
  template <typename InputIteratorT>
  Permutation(InputIteratorT first, InputIteratorT last) : permutation(first, last) {
    assert(validate());
  }

  // Constructs a permutation from the specified binary file.
  explicit Permutation(std::ifstream& in) {
    readFrom(in);
  }

  // Returns a random permutation of the specified size.
  template <typename RandomNumberGeneratorT>
  static Permutation getRandomPermutation(const int size, RandomNumberGeneratorT&& rand) {
    Permutation perm(size);
    std::iota(perm.permutation.begin(), perm.permutation.end(), 0);
    std::shuffle(perm.permutation.begin(), perm.permutation.end(), rand);
    return perm;
  }

  // Replaces this permutation with the one specified as an iterator range.
  template <typename InputIteratorT>
  void assign(InputIteratorT first, InputIteratorT last) {
    permutation.assign(first, last);
    assert(validate());
  }

  // Returns an iterator referring to the new location of the element at location 0.
  ConstIterator begin() const noexcept {
    return permutation.begin();
  }

  // Returns an iterator which is the past-the-end value for this permutation.
  ConstIterator end() const noexcept {
    return permutation.end();
  }

  // Returns the size of this permutation.
  int size() const noexcept {
    return permutation.size();
  }

  // Returns the new location of the element at oldPos.
  int& operator[](const int oldPos) {
    assert(oldPos >= 0); assert(oldPos < permutation.size());
    return permutation[oldPos];
  }

  // Returns the new location of the element at oldPos.
  int operator[](const int oldPos) const {
    assert(oldPos >= 0); assert(oldPos < permutation.size());
    return permutation[oldPos];
  }

  // Returns true if each element is mapped to the same location in the two permutations.
  friend bool operator==(const Permutation& lhs, const Permutation& rhs) {
    return lhs.permutation == rhs.permutation;
  }

  // Inverts the permutation.
  void invert() {
    *this = getInversePermutation();
  }

  // Returns the inverse permutation that maps new locations to old locations.
  Permutation getInversePermutation() const {
    assert(validate());
    Permutation perm(size());
    for (int oldPos = 0; oldPos < permutation.size(); ++oldPos)
      perm[permutation[oldPos]] = oldPos;
    return perm;
  }

  // Reorders the elements in the specified container according to this permutation.
  template <typename ContT>
  void applyTo(ContT& cont) const {
    assert(validate());
    assert(cont.size() == size());
    ContT temp(cont.size());
    for (int i = 0; i < cont.size(); ++i)
      temp[permutation[i]] = std::move(cont[i]);
    cont.swap(temp);
  }

  // Checks if this permutation indeed represents a valid permutation.
  bool validate() const {
    const auto size = permutation.size();
    BitVector isContained(size);
    for (int i = 0; i < size; ++i) {
      assert(!isContained[permutation[i]]);
      isContained[permutation[i]] = true;
    }
    return true;
  }

  // Reads the permutation from the specified binary file.
  void readFrom(std::ifstream& in) {
    bio::read(in, permutation);
  }

  // Writes the permutation to the specified binary file.
  void writeTo(std::ofstream& out) const {
    bio::write(out, permutation);
  }

 private:
  std::vector<int32_t> permutation; // permutation[i] is the new location of the element at i.
};
