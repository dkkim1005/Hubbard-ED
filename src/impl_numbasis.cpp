#include "../include/numbasis.hpp"

namespace Fermion {
Fockstate::Fockstate(const unsigned int nsites, const unsigned int nparticles):
  data_(::combination(2u*nsites, nparticles)),
  knsites(nsites),
  knparticles(nparticles)
{
  unsigned int s0 = 0u;
  for (int i=0; i<knparticles; ++i)
    s0 = (s0 << 1)+1;
  for (int i=0; i<data_.size(); ++i)
  {
    data_[i] = s0;
    s0 = bittool::next_set_of_n_elements(s0);
  }
}

Bitpair Fockstate::operator()(const unsigned int idx) const
{
  return Bitpair(data_[idx] >> knsites, data_[idx]-((data_[idx] >> knsites) << knsites));
}

void Fockstate::print(const unsigned int idx) const
{
  std::cout << "|";
  for (unsigned int i=0u; i<knsites; i++)
    std::cout << bittool::bitget(2u*knsites-1u-i, data_[idx]);
  std::cout << "> (x) |";
  for (unsigned int i=knsites; i<2u*knsites; i++)
    std::cout << bittool::bitget(2u*knsites-1u-i, data_[idx]);
  std::cout << ">" << std::endl;
}

unsigned int Fockstate::search(const unsigned int & s) const
{
  assert(bittool::bitcount(s) == knparticles);
  return this->binary_search_(s);
}

unsigned int Fockstate::search(const Bitpair & spair) const
{
  //std::cout << (bittool::bitcount(spair.first)+bittool::bitcount(spair.second)) << std::endl;
  assert((bittool::bitcount(spair.first)+bittool::bitcount(spair.second)) == knparticles);
  unsigned int s = (spair.first << knsites) + spair.second;
  return this->binary_search_(s);
}

unsigned int Fockstate::binary_search_(const unsigned int s) const
{
  unsigned int idx = (data_.size()-1u)/2u, idx0 = 0u, idx1 = data_.size()-1u;
  while (1)
  {
    if (data_[idx] < s)
    {
      idx0 = idx+1;
      idx = (idx0+idx1)/2u;
    }
    else if (data_[idx] > s)
    {
      idx1 = idx-1;
      idx = (idx0+idx1)/2u;
    }
    else if (data_[idx] == s)
      break;
  }
  return idx;
}

op::op(const unsigned int nsites):
  knsites(nsites),
  lookup_(nsites)
{
  for (int i=0; i<lookup_.size(); ++i)
    lookup_[i] = ((i%2==0) ? 1 : -1);
}

BitInt op::add(unsigned int s, unsigned int idx) const
{
  const bool isEmpty = (1u-bittool::bitget(idx, s));
  const int phase = lookup_[bittool::bitcount(s >> (idx+1u))];
  return BitInt(bittool::bitflip(idx, s)*isEmpty, phase);
}

BitInt op::del(unsigned int s, unsigned int idx) const
{
  const bool isFull = bittool::bitget(idx, s);
  const int phase = lookup_[bittool::bitcount(s >> (idx+1u))];
  return BitInt(bittool::bitflip(idx, s)*isFull, phase);
}
} // end namespace Fermion
