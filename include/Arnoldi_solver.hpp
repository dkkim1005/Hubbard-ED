// Copyright (c) 2020-2021 Dongkyu Kim (dkkim1005@gmail.com)

#pragma once

#include <Eigen/Sparse>
#include <SymEigsSolver.h>
#include <MatOp/SparseSymMatProd.h>
#include <iostream>
#include <cstdlib>

namespace EIGEN_SOLVER {
namespace EIGEN {
namespace SPECTRA {
typedef Eigen::SparseMatrix<double, Eigen::ColMajor, int> SparseDoubleInt;
bool SymEigsSolver(const SparseDoubleInt& A, Eigen::VectorXd& evalues, Eigen::MatrixXd& evectors,
  const int nev, const int ncv, const bool info = false)
{
  /* ncv must satisfy nev < ncv <= n, n is the size of matrix
     This parameter must satisfy nev<ncv<=n, and is advised to take ncv>=2nev.*/
  // Construct matrix operation object using the wrapper class SparseGenMatProd
  Spectra::SparseSymMatProd<double> op(A);
  // Construct eigen solver object, requesting the smallest eigenvalues
  Spectra::SymEigsSolver<double, Spectra::SMALLEST_ALGE, Spectra::SparseSymMatProd<double> > eigs(&op, nev, ncv);
  eigs.init();
  int nconv = eigs.compute();
  bool isconverged = true;
  if(info)
    std::cout << "      The # of converged eigenvalues: " << nconv << " (nev: " << nev << ")\n";
  if(eigs.info() != Spectra::SUCCESSFUL)
  {
    std::cout << "   !Error: EIGEN_SOLVER::EIGEN::SPECTRA::SymEigsSolver\n"
              << "           Solver fails to gain results...\n";
    isconverged = false;
  }
  evalues = eigs.eigenvalues();
  evectors = eigs.eigenvectors(nev);
  return isconverged;
}
} // end namespace SPECTRA 
} // end namespace EIGEN 
} // end namespace EIGEN_SOLVER 
