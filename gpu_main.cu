#include <iostream>
#include "./include/Hamiltonian.hpp"
#include "./include/argparse.hpp"
#include "./include/magma_eigen_solver.cuh"

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

  int * col = rawData.innerIndexPtr(),
    * row = rawData.outerIndexPtr();
  double * val = rawData.valuePtr();

  magma_init();
  magma_dopts opts;
  magma_queue_t queue;
  magma_queue_create(0, &queue);

  MAGMA_DLOBPCG_WRAPPER solver(HilbertSize);
  solver.construct_CSR_format(row, col, val, queue);

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
    std::cout << opts.solver_par.eigenvalues[i] << " ";
  std::cout << std::endl;

  solver.free(queue);

  magma_queue_destroy(queue);
  magma_finalize();

  return 0;
}
