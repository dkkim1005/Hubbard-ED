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
  options.push_back(pair_t("mu1", "chemical potential (spin up)"));
  options.push_back(pair_t("mu2", "chemical potential (spin down)"));
  // env; default value
  defaults.push_back(pair_t("t", "1"));
  // parser for arg list
  argsparse parser(argc, argv, options, defaults);
  const unsigned int lx = parser.find<unsigned int>("lx"),
    ly = parser.find<unsigned int>("ly");
  const unsigned int nsites = lx*ly,
    nparticles = parser.find<unsigned int>("np");
  const unsigned int HilbertSize = ::combination(2*nsites, nparticles);
  const double t = parser.find<double>("t"),
    U = parser.find<double>("U"),
    mu1 = parser.find<double>("mu1"),
    mu2 = parser.find<double>("mu2");
  EigenSparseMatrix smatrix(HilbertSize);
  SquareLattice lattice(lx, ly);
  HubbardHamiltonian Hamiltonian(nsites, nparticles);
  std::cout << "# construct sparse matrix... " << std::flush;
  Hamiltonian.construct_matrix(smatrix, lattice, t, U, mu1, mu2);
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
  std::cout << "polarization: " << Hamiltonian.meas_polarization(&evectors(0)) << std::endl;
  return 0;
}
