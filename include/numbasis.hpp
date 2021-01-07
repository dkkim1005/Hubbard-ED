// Copyright (c) 2020-2021 Dongkyu Kim (dkkim1005@gmail.com)

#pragma once

#include <iostream>
#include <vector>
#include <assert.h>
#include "bittool.hpp"

namespace {
// combination nCr=n!/(n-r)!/r!
inline unsigned int combination(const unsigned int n, const unsigned int r)
{
  unsigned int r0 = ((r >= n/2) ? r : n-r);
  unsigned long int num = 1ul;
  for (unsigned int i=(r0+1); i<=n; ++i)
    num *= i;
  for (unsigned int i=1; i<=(n-r0); ++i)
    num /= i;
  return num;
}
} // end namespace

// paring of two Fock states of Fermion particles: e.g. spin up and down
typedef std::pair<unsigned int, unsigned int> Bitpair;

namespace Fermion {
// the number of sites and particles are fixed.
class Fockstate
{
public:
  // fix the total number of particles
  Fockstate(const unsigned int nsites, const unsigned int nparticles);
  // fix the number of particles for each spin state
  Fockstate(const unsigned int nsites, const unsigned int nspinsup, const unsigned int nspinsdw);
  Fockstate(const Fockstate & rhs) = delete;
  Fockstate & operator=(const Fockstate & rhs) = delete;
  // size of the Hilbert space 
  unsigned int size() const { return data_.size(); }
  // return the state at the index 'idx'
  unsigned int operator[](const unsigned int idx) const { return data_[idx]; }
  // return the state with the Bitpair type at the index 'idx'
  Bitpair operator()(const unsigned int idx) const;
  // print the sate at the index 'idx'
  void print(const unsigned int idx) const;
  // return the index for the state 's'
  unsigned int search(const unsigned int & s) const;
  // return the index for the state 's' in terms of the Bitpair type
  unsigned int search(const Bitpair & spair) const;
  unsigned get_nsites() const { return knsites; }
  unsigned get_nparticles() const { return knparticles; }
private:
  unsigned int binary_search_(const unsigned int s) const;
  std::vector<unsigned int> data_;
  const unsigned int knsites, knparticles;
};

// pairing of a Fock state and phase that can be two values, 1 or -1.
typedef std::pair<unsigned int, int> BitInt;

// object for creation and annihilation operators
class op
{
public:
  explicit op(const unsigned int nsites);
  // return a Fock state, where the particle is created at the index 'idx'
  BitInt add(unsigned int s, unsigned int idx) const;
  // return a Fock state, where the particle is annihilated at the index 'idx'
  BitInt del(unsigned int s, unsigned int idx) const;
private:
  const unsigned int knsites;
  std::vector<int> lookup_;
};
} // end namespace Fermion
