// Copyright (c) 2020-2021 Dongkyu Kim (dkkim1005@gmail.com)

#include <iostream>
#include <memory>
#include <exception>
#include <fstream>
#include "./include/Hamiltonian.hpp"
#include "./include/argparse.hpp"
#include "./include/magma_eigen_solver.cuh"

Fermion::Fockstate * basis_factory(const int L, const std::vector<int> & np);
std::vector<double> generate_harmonic_potential(const int L, const double V);

int main(int argc, char * argv[])
{
  std::vector<pair_t> options, defaults;
  // env; explanation of env
  options.push_back(pair_t("L", "system size"));
  options.push_back(pair_t("np", "# of particle numbers; one elem -> fix total #, two elem -> fix # for each spin"));
  options.push_back(pair_t("t", "hopping element"));
  options.push_back(pair_t("U", "onsite interaction"));
  options.push_back(pair_t("V", "strength of the harmonic potential"));
  options.push_back(pair_t("bias", "mat <-- mat+bias*I"));
  options.push_back(pair_t("pbc", "use periodic boundary condition (true : 1 or false : 0)"));
  // env; default value
  defaults.push_back(pair_t("t", "1"));
  defaults.push_back(pair_t("bias", "0"));
  // parser for arg list
  argsparse parser(argc, argv, options, defaults);
  const unsigned int L = parser.find<unsigned int>("L");
  const unsigned int nsites = L;
  const auto np = parser.mfind<int>("np");
  const double t = parser.find<double>("t"),
    U = parser.find<double>("U"),
    bias = parser.find<double>("bias");
  const bool usePBC = parser.find<bool>("pbc");

  // number basis (Fock space of fixed # of particles)
  std::unique_ptr<Fermion::Fockstate> basis_ptr(basis_factory(L, np));
  const unsigned int HilbertSize = basis_ptr->size();
  const std::vector<double> V = generate_harmonic_potential(L, parser.find<double>("V"));
  // matrix container to deploy Hubbard Hamiltonian
  EigenSparseMatrix smatrix(HilbertSize);
  // lattice information
  ChainLattice lattice(nsites, usePBC);

  std::cout << "# construct sparse matrix... (hilbert size: " << HilbertSize << ") " << std::flush;
  // parameter set of Hubbard model
  HubbardParams params;
  params.t = t;
  params.U = U;
  params.V = V;
  construct_hubbard_hamiltonian(smatrix, *basis_ptr, lattice, params, bias);
  std::cout << "done." << std::endl << std::flush;

  // Eigen sparse --> MAGMA sparse
  auto & rawData = smatrix.get_data();
  rawData.prune(1e-15);
  rawData.makeCompressed();
  int * col = rawData.innerIndexPtr(),
    * row = rawData.outerIndexPtr();
  double * val = rawData.valuePtr();
  magma_init();
  magma_dopts opts;
  magma_queue_t queue;
  magma_queue_create(0, &queue);
  MAGMA_DLOBPCG_WRAPPER solver(HilbertSize);
  solver.construct_CSR_format(row, col, val, queue);

  // ********** MAGMA-GPU eigen value solver **********

  // options for DLOBPCG solver
  opts.solver_par.solver = Magma_LOBPCG;
  opts.solver_par.num_eigenvalues = 1; // number of eigenvalues you want to compute
  opts.solver_par.ev_length = HilbertSize;
  opts.solver_par.maxiter = 1000; // max number of iterations
  opts.solver_par.verbose = true;
  opts.solver_par.rtol = 1e-7;
  opts.solver_par.restart = 8;
  // options for preconditioner
  opts.precond_par.solver = Magma_JACOBI;
  opts.precond_par.levels = 0; // ILU(0) - no fill-in
  opts.precond_par.trisolver = Magma_CUSOLVE; //exact triangular solves

  magma_deigensolverinfo_init(&opts.solver_par, queue);

  solver.run(opts, queue);

  std::cout << "--- lowest energy:" << std::endl;
  for (int i=0; i<opts.solver_par.num_eigenvalues; ++i)
    std::cout << (opts.solver_par.eigenvalues[i]-bias) << " ";
  std::cout << std::endl;

  std::vector<double> eigenvectors(HilbertSize);
  if (cudaSuccess != cudaMemcpy(eigenvectors.data(), opts.solver_par.eigenvectors,
      sizeof(double)*eigenvectors.size(), cudaMemcpyDeviceToHost))
    std::cerr << "check cudaMemcpy!" << (__LINE__-2) << std::endl;

  solver.free(queue);

  magma_queue_destroy(queue);
  magma_finalize();

  // **************************************************

  // plot density profile
  const auto density = meas_density(*basis_ptr, eigenvectors);
  const std::string fname = "density-L" + parser.find<>("L") + "U" + parser.find<>("U")
    + "V" + parser.find<>("V") + "P" + parser.find<>("pbc") + "N" + parser.find<>("np") + ".dat";
  std::ofstream wfile(fname);
  for (int i=0; i<L; ++i)
    wfile << (i+1) << " " << density[i] << std::endl;

  return 0;
}

Fermion::Fockstate * basis_factory(const int L, const std::vector<int> & np)
{
  if (np.size() == 1)
    return new Fermion::Fockstate(L, np[0]);
  else if (np.size() == 2)
    return new Fermion::Fockstate(L, np[0], np[1]);
  else
    throw std::invalid_argument("Check the size of'np' container.");
}

std::vector<double> generate_harmonic_potential(const int L, const double V)
{
  std::vector<double> tmp(L);
  for (int i=0; i<L; ++i)
    tmp[i] = V*std::pow(i-(L-1.0)/2.0, 2);
  return tmp;
}
