// Copyright (c) 2020-2021 Dongkyu Kim (dkkim1005@gmail.com)

#include <iostream>
#include "./include/Hamiltonian.hpp"
#include "./include/Arnoldi_solver.hpp"
#include "./include/argparse.hpp"

int main(int argc, char * argv[])
{
  std::vector<pair_t> options, defaults;
  // env; explanation of env
  options.push_back(pair_t("lx", "width of x-direction"));
  options.push_back(pair_t("ly", "width of y-direction"));
  options.push_back(pair_t("np", "# of particle numbers"));
  options.push_back(pair_t("t", "hopping element"));
  options.push_back(pair_t("U", "onsite interaction"));
  // env; default value
  defaults.push_back(pair_t("t", "1"));
  // parser for arg list
  argsparse parser(argc, argv, options, defaults);
  const unsigned int lx = parser.find<unsigned int>("lx"),
    ly = parser.find<unsigned int>("ly");
  const unsigned int nsites = lx*ly,
    nparticles = parser.find<unsigned int>("np");
  const double t = parser.find<double>("t"),
    U = parser.find<double>("U");

  // number basis (Fock space of fixed # of particles)
  Fermion::Fockstate nbasis(nsites, nparticles);
  const unsigned int HilbertSize = nbasis.size();
  // matrix container to deploy Hubbard Hamiltonian
  EigenSparseMatrix smatrix(HilbertSize);
  // lattice information
  SquareLattice lattice(lx, ly);

  std::cout << "# construct sparse matrix... " << std::flush;
  // parameter set of Hubbard model
  HubbardParams params;
  params.t = t;
  params.U = U;
  params.V = std::vector<double>(nsites, 0.0);
  construct_hubbard_hamiltonian(smatrix, nbasis, lattice, params);
  std::cout << "done." << std::endl << std::flush;

  auto & rawData = smatrix.get_data();
  rawData.prune(1e-15);
  rawData.makeCompressed();
  Eigen::VectorXd evalues;
  Eigen::MatrixXd evectors;
  const int nev = 1;
  std::cout << "# solve an eigenvalue problem... " << std::flush;
  EIGEN_SOLVER::EIGEN::SPECTRA::SymEigsSolver(smatrix.get_data(), evalues, evectors, nev, 30);
  std::cout << "done." << std::endl << std::flush;
  std::cout << "energy: " << evalues << std::endl;

  return 0;
}
