#include "../include/Hamiltonian.hpp"

HubbardHamiltonian::HubbardHamiltonian(const unsigned int nsites, const unsigned int nparticles):
  nbasis_(nsites, nparticles),
  op_(nsites),
  knsites(nsites),
  knparticles(nparticles)
{}

void HubbardHamiltonian::construct_matrix(BaseMatrixType & mat, const BaseLattice & NN,
  const double t, const double U, const double mu1, const double mu2)
{
  assert(mat.width() == nbasis_.size());
  for (unsigned int idx_r = 0u; idx_r<nbasis_.size(); ++idx_r)
  {
    const auto & state = nbasis_(idx_r);
    // hopping
    for (unsigned int i=0; i<NN.get_nsites(); ++i)
    {
      // \sum_{j}\sum_{i} c^+_{j,up} c_{i,up} |n>
      for (unsigned int j=0; j<NN.get_nnnsites(); ++j)
      {
        if (bittool::bitget(i, state.first) != 1)
          break;
        const auto & state0_up = op_.del(state.first, i);
	if (bittool::bitget(NN.get_nnsite(i, j), state0_up.first) != 0)
          continue;
        const auto & state1_up = op_.add(state0_up.first, NN.get_nnsite(i, j));
        const int phase = (state0_up.second*state1_up.second);
        const unsigned int idx_l = nbasis_.search(Bitpair(state1_up.first, state.second));
        mat.add(idx_l, idx_r, -t*phase);
      }
      // \sum_{j}\sum_{i} c^+_{j,dw} c_{i,dw} |n>
      for (unsigned int j=0; j<NN.get_nnnsites(); ++j)
      {
        if (bittool::bitget(i, state.second) != 1)
          break;
        const auto & state0_dw = op_.del(state.second, i);
        if (bittool::bitget(NN.get_nnsite(i, j), state0_dw.first) != 0)
          continue;
        const auto & state1_dw = op_.add(state0_dw.first, NN.get_nnsite(i, j));
        const int phase = (state0_dw.second*state1_dw.second);
        const unsigned int idx_l = nbasis_.search(Bitpair(state.first, state1_dw.first));
        mat.add(idx_l, idx_r, -t*phase);
      }
    }
    // on-site interaction
    double diag = 0;
    for (unsigned int i=0; i<NN.get_nsites(); ++i)
      diag += U*bittool::bitget(i, state.first)*bittool::bitget(i, state.second);
    mat.add(idx_r, idx_r, diag);
    // Zeeman field
    diag = 0;
    for (unsigned int i=0; i<NN.get_nsites(); ++i)
      diag += (mu1*bittool::bitget(i, state.first))+(mu2*bittool::bitget(i, state.second));
    mat.add(idx_r, idx_r, diag);
  }
}

double HubbardHamiltonian::meas_polarization(const double * eigvec) const
{
  double dnavg = 0;
  for (unsigned int i=0u; i<nbasis_.size(); ++i)
  {
    const auto & state = nbasis_(i);
    // dn = n_dw - n_up
    double dn = 0;
    for (unsigned int j=0u; j<knsites; ++j)
      dn += (int)bittool::bitget(j, state.second)-(int)bittool::bitget(j, state.first);
    dnavg += eigvec[i]*eigvec[i]*dn;
  }
  return dnavg/knparticles;
}
