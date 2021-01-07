// Copyright (c) 2020-2021 Dongkyu Kim (dkkim1005@gmail.com)

#include "../include/Hamiltonian.hpp"

void construct_hubbard_hamiltonian(BaseMatrixType & mat,
  const Fermion::Fockstate & nbasis, const BaseLattice & NN,
  const HubbardParams & params, const double bias)
{
  const Fermion::op ops(nbasis.get_nsites());
  assert(mat.width() == nbasis.size());
  assert(NN.get_nsites() == params.V.size());
  for (unsigned int idx_r = 0u; idx_r<nbasis.size(); ++idx_r)
  {
    const auto & state = nbasis(idx_r);
    // hopping
    for (unsigned int i=0; i<NN.get_nsites(); ++i)
    {
      // \sum_{j}\sum_{i} c^+_{j,up} c_{i,up} |n>
      for (const int & nnidx : NN.get_nnsite(i))
      {
        if (bittool::bitget(i, state.first) != 1)
          break;
        const auto & state0_up = ops.del(state.first, i);
        if (bittool::bitget(nnidx, state0_up.first) != 0)
          continue;
        const auto & state1_up = ops.add(state0_up.first, nnidx);
        const int phase = (state0_up.second*state1_up.second);
        const unsigned int idx_l = nbasis.search(Bitpair(state1_up.first, state.second));
        mat.add(idx_l, idx_r, -params.t*phase);
      }
      // \sum_{j}\sum_{i} c^+_{j,dw} c_{i,dw} |n>
      for (const int & nnidx : NN.get_nnsite(i))
      {
        if (bittool::bitget(i, state.second) != 1)
          break;
        const auto & state0_dw = ops.del(state.second, i);
        if (bittool::bitget(nnidx, state0_dw.first) != 0)
          continue;
        const auto & state1_dw = ops.add(state0_dw.first, nnidx);
        const int phase = (state0_dw.second*state1_dw.second);
        const unsigned int idx_l = nbasis.search(Bitpair(state.first, state1_dw.first));
        mat.add(idx_l, idx_r, -params.t*phase);
      }
    }
    // on-site interaction
    double diag = 0;
    for (unsigned int i=0; i<NN.get_nsites(); ++i)
      diag += params.U*bittool::bitget(i, state.first)*bittool::bitget(i, state.second);
    mat.add(idx_r, idx_r, diag);
    // Zeeman field
    diag = 0;
    for (unsigned int i=0; i<NN.get_nsites(); ++i)
      diag += params.V[i]*(bittool::bitget(i, state.first)+bittool::bitget(i, state.second));
    mat.add(idx_r, idx_r, diag);
    // mat <-- mat + bias*I
    mat.add(idx_r, idx_r, bias);
  }
}

double meas_polarization(const Fermion::Fockstate & nbasis, const std::vector<double> & eigvec)
{
  assert(nbasis.size() == eigvec.size());
  double dnavg = 0;
  for (unsigned int i=0u; i<nbasis.size(); ++i)
  {
    const auto & state = nbasis(i);
    // dn = n_dw - n_up
    double dn = 0;
    for (unsigned int j=0u; j<nbasis.get_nsites(); ++j)
      dn += (int)bittool::bitget(j, state.second)-(int)bittool::bitget(j, state.first);
    dnavg += dn*std::pow(eigvec[i], 2);
  }
  return dnavg/2.0;
}

std::vector<double> meas_density(const Fermion::Fockstate & nbasis, const std::vector<double> & eigvec)
{
  assert(nbasis.size() == eigvec.size());
  std::vector<double> density(nbasis.get_nsites(), 0.0);
  for (unsigned int idx=0u; idx<nbasis.size(); ++idx)
  {
    const auto & state = nbasis(idx);
    for (unsigned int i=0u; i<nbasis.get_nsites(); ++i)
      density[i] += (bittool::bitget(i, state.first)+bittool::bitget(i, state.second))*std::pow(eigvec[idx], 2);
  }
  return density;
}
