// Copyright (c) 2020-2021 Dongkyu Kim (dkkim1005@gmail.com)

#pragma once

#include <iostream>
#include <vector>
#include <chrono>
#include <random>
#include <cmath>
#include <magma_v2.h>
#include <magmasparse.h>

#define CHECK_ERROR(SUCCESS_MACRO, STATUS) do {\
  if (SUCCESS_MACRO != (STATUS)) {\
    std::cerr << "# ERROR --- FILE:" << __FILE__ << ", LINE:" << __LINE__ << std::endl;\
    exit(1);\
  }\
} while (false)

// CSR format is employed here.
class MAGMA_DLOBPCG_WRAPPER
{
/***********************************************************
  // parameter settings for DLOBPCG solver and preconditioner

  magma_dopts opts;

  // options for DLOBPCG solver
  opts.solver_par.solver = Magma_LOBPCG;
  opts.solver_par.num_eigenvalues = 1; // number of eigenvalues you want to compute
  opts.solver_par.ev_length = m; // length of an eigenvector 
  opts.solver_par.maxiter = 1000; // max number of iterations
  opts.solver_par.verbose = true;
  opts.solver_par.rtol = 1e-7;
  opts.solver_par.restart = 8;

  // options for preconditioner
  opts.precond_par.solver = Magma_JACOBI; // note that the diagonal part of a matrix should be non-zero.
  opts.precond_par.levels = 0; // ILU(0) - no fill-in
  opts.precond_par.trisolver = Magma_CUSOLVE; //exact triangular solves

  magma_deigensolverinfo_init(opts.solver_par, queue);
***********************************************************/
public:
  explicit MAGMA_DLOBPCG_WRAPPER(const int m);
  void construct_CSR_format(int * row, int * col, double * val, const magma_queue_t & queue);
  void run(magma_dopts & opts, const magma_queue_t & queue);
  void free(const magma_queue_t & queue);
private:
  const int m_;
  magma_d_matrix A_h_, A_d_, b_h_, b_d_;
  bool IsAllocated;
};
