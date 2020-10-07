#pragma once

#include <Eigen/Sparse>
#include "numbasis.hpp"

class BaseMatrixType
{
public:
  explicit BaseMatrixType(const unsigned int n): kn(n) {}
  unsigned int width() const { return kn; }
  virtual void add(const unsigned int idx_l, const unsigned int idx_r, const double val) = 0;
protected:
  const unsigned int kn;
};

class EigenSparseMatrix : public BaseMatrixType
{
  typedef Eigen::SparseMatrix<double, Eigen::ColMajor, int> SpMat;
public:
  explicit EigenSparseMatrix(const unsigned int n):
    BaseMatrixType(n),
    smat_(n, n) {}
  virtual ~EigenSparseMatrix() {}
  virtual void add(const unsigned int idx_l, const unsigned int idx_r, const double val)
  {
    smat_.coeffRef(idx_l, idx_r) = smat_.coeffRef(idx_l, idx_r) + val;
  }
  SpMat & get_data() { return smat_; }
private:
  SpMat smat_;
};


class BaseLattice
{
public:
  BaseLattice(const unsigned int nsites, const unsigned int nnsites):
    knsites(nsites),
    knnsites(nnsites),
    data_(nsites, std::vector<int>(nnsites)) {}
  unsigned int get_nnsite(const unsigned int idx, const unsigned nnidx) const
  {
    return data_[idx][nnidx];
  }
  unsigned int get_nsites() const { return knsites; }
  unsigned int get_nnnsites() const { return knnsites; }
protected:
  const unsigned int knsites, knnsites;
  std::vector<std::vector<int> > data_;
};


class ChainLattice : public BaseLattice
{
public:
  ChainLattice(const unsigned int nsites):
    BaseLattice(nsites, 2)
  {
    for (int i=0; i<knsites; ++i)
    {
      data_[i][0] = ((i == 0) ? knsites-1 : i-1);
      data_[i][1] = ((i == knsites-1) ? 0 : i+1);
    }
  }
};

class SquareLattice : public BaseLattice
{
public:
  SquareLattice(const unsigned int nx, const unsigned int ny):
    BaseLattice(nx*ny, 4)
  {
    for (int i=0; i<ny; ++i)
    {
      for (int j=0; j<nx; ++j)
      {
        // up
        data_[i*nx+j][0] = ((i == 0) ? (ny-1)*nx+j : (i-1)*nx+j);
	// down
        data_[i*nx+j][1] = ((i == ny-1) ? j : (i+1)*nx+j);
	// left
        data_[i*nx+j][2] = ((j == 0) ? i*nx+nx-1 : i*nx+j-1);
	// right
        data_[i*nx+j][3] = ((j == nx-1) ? i*nx : i*nx+j+1);
      }
    }
  }
};


// H = -t\sum_{ijs}(c^+_{i,s} c_{j,s}) + U\sum_{i}n_{i,up}n_{i,dw} + h\sum_{i}(c^+_{i,up}c_{i,up} - c^+_{i,dw}c_{i,dw})
class HubbardHamiltonian
{
public:
  explicit HubbardHamiltonian(const unsigned int nsites,
    const unsigned int nparticles);
  void construct_matrix(BaseMatrixType & mat, const BaseLattice & NN,
    const double t, const double U, const double mu1, const double mu2);
  double meas_polarization(const double * eigvec) const;
private:
  Fermion::Fockstate nbasis_;
  Fermion::op op_;
  const unsigned knsites, knparticles;
};
