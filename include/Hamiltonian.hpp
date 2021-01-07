// Copyright (c) 2020-2021 Dongkyu Kim (dkkim1005@gmail.com)

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
  explicit BaseLattice(const unsigned int nsites):
    knsites(nsites),
    data_(nsites) {}
  // nearest neighbor sites for lattice site index 'idx'
  const std::vector<int> & get_nnsite(const unsigned int idx) const { return data_[idx]; }
  // return total # of lattice sites
  unsigned int get_nsites() const { return knsites; }
protected:
  const unsigned int knsites;
  // look-up table of nearest neighbor sites
  std::vector<std::vector<int> > data_;
};

class ChainLattice : public BaseLattice
{
public:
  ChainLattice(const unsigned int nsites, const bool usePBC):
    BaseLattice(nsites)
  {
    for (auto & item : data_)
      item.resize(2);
    for (int i=0; i<knsites; ++i)
    {
      data_[i][0] = ((i == 0) ? knsites-1 : i-1);
      data_[i][1] = ((i == knsites-1) ? 0 : i+1);
    }
    // open boundary condition
    if (!usePBC)
    {
      data_[0].resize(1);
      data_[knsites-1].resize(1);
      data_[0][0] = 1;
      data_[knsites-1][0] = knsites-2;
    }
  }
};

class SquareLattice : public BaseLattice
{
public:
  SquareLattice(const unsigned int nx, const unsigned int ny):
    BaseLattice(nx*ny)
  {
    for (auto & item : data_)
      item.resize(4);
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


// H = -t\sum_{ijs}(c^+_{i,s} c_{j,s}) + U\sum_{i}n_{i,up}n_{i,dw} + \sum_{i}V_{i}(c^+_{i,up}c_{i,up} + c^+_{i,dw}c_{i,dw})
struct HubbardParams
{
  double t, U;
  std::vector<double> V;
};

// return a matrix form of Hubbard model for given info
void construct_hubbard_hamiltonian(BaseMatrixType & mat, const Fermion::Fockstate & nbasis,
  const BaseLattice & NN, const HubbardParams & params, const double bias = 0);

// meas: S_{z} = \sum_{i} (c^+_{i,up}c_{i,up} - c^+_{i,dw}c_{i,dw})/2
double meas_polarization(const Fermion::Fockstate & nbasis, const std::vector<double> & eigvec);

// meas: n_{i} = c^+_{i,up}c_{i,up} + c^+_{i,dw}c_{i,dw}
std::vector<double> meas_density(const Fermion::Fockstate & nbasis, const std::vector<double> & eigvec);
