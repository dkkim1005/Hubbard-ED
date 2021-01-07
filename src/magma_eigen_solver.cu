// Copyright (c) 2020-2021 Dongkyu Kim (dkkim1005@gmail.com)

#include "../include/magma_eigen_solver.cuh"

MAGMA_DLOBPCG_WRAPPER::MAGMA_DLOBPCG_WRAPPER(const int m):
  m_(m),
  A_h_{Magma_CSR},
  A_d_{Magma_CSR},
  b_h_{Magma_CSR},
  b_d_{Magma_CSR},
  IsAllocated(false) {}

void MAGMA_DLOBPCG_WRAPPER::construct_CSR_format(int * row, int * col, double * val, const magma_queue_t & queue)
{
  if (IsAllocated)
  {
    std::cerr << "plz excute 'free' method first." << std::flush << std::endl;
    return;
  }
  // convert to the csr format
  magma_dcsrset(m_, m_, row, col, val, &A_h_, queue);
  // initialize a vector from a random starting point
  std::vector<double> rhs(m_);
  unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::mt19937 gen(seed);
  std::uniform_real_distribution<double> dist(-1.0, 1.0);
  double norm = 0;
  for (auto & item : rhs)
  {
    item = dist(gen);
    norm += item*item;
  }
  norm = std::sqrt(norm);
  for (auto & item : rhs)
    item /= norm;
  magma_dvset(m_, 1, rhs.data(), &b_h_, queue);
  // cpu --> gpu
  magma_dmtransfer(A_h_, &A_d_, Magma_CPU, Magma_DEV, queue);
  magma_dmtransfer(b_h_, &b_d_, Magma_CPU, Magma_DEV, queue);
  IsAllocated = true;
}

void MAGMA_DLOBPCG_WRAPPER::run(magma_dopts & opts, const magma_queue_t & queue)
{
  magma_d_precondsetup(A_d_, b_d_, &opts.solver_par, &opts.precond_par, queue);
  magma_dlobpcg(A_d_, &opts.solver_par, &opts.precond_par, queue);
}

void MAGMA_DLOBPCG_WRAPPER::free(const magma_queue_t & queue)
{
  if (IsAllocated == true)
  {
    magma_dmfree(&A_h_, queue);
    magma_dmfree(&b_h_, queue);
    magma_dmfree(&A_d_, queue);
    magma_dmfree(&b_d_, queue);
  }
  IsAllocated = false;
}
