// 
//  Copyright (C) 2007, 2021  Smithsonian Astrophysical Observatory
//
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with this program; if not, write to the Free Software Foundation, Inc.,
//  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
//
#ifndef myArray_hh
#define myArray_hh

#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace sherpa {

#ifndef NDEBUG

#define MYCMPSIZES(arg) this->check_size(arg)
#define STRINGIFY(x) #x
#define PRINTARRAY4(os, arg, num, suffix)                                      \
  arg.myprint(os, num, STRINGIFY(arg), suffix);
#define PRINTARRAY(arg) PRINTARRAY4(std::cout, arg, 6, "\n")

#else

// production code should run silently!
#define MYCMPSIZES(arg)
#define PRINTARRAY4(os, arg, suffix)
#define PRINTARRAY(arg)
// production code should run silently!

#endif

// template<typename T>
// struct Cmp {
//   virtual bool operator()(const T& lhs, const T& rhs) const {
//     return true;
//   }
// };

////////////////////////////////// Array1D //////////////////////////////////
template <typename T> class Array1D {

  friend std::ostream &operator<<(std::ostream &s, const Array1D<T> &A) {
    return A.print(s);
  }

public:
  Array1D() : vec() {}

  Array1D(int arg, T init = 0) : vec(arg, init) {}

  Array1D(T *a, T *b) : vec(a, b) {}

  Array1D(const Array1D &arg) { vec = arg.get_vec(); }

  virtual bool operator<(const Array1D<T> &rhs) const {
    int num = vec.size() - 1;
    return vec[num] < rhs[num];
  }

  const Array1D &operator=(const Array1D &rhs) {
    if (&rhs != this)
      vec = rhs.get_vec();
    return *this;
  }

  T *operator()() const { return &vec[0]; }

  operator T *() { return &vec[0]; }

  T &operator[](int arg) {
#ifndef NDEBUG
    return vec.at(arg);
#endif
    return vec[arg];
  }
  const T &operator[](int arg) const {
#ifndef NDEBUG
    return vec.at(arg);
#endif
    return vec[arg];
  }

  // std::vector<T>::iterator begin() { return vec.begin(); }//
  // std::vector<T>::const_iterator begin() const { return vec.begin(); }
  void resize(std::vector<int>::size_type n) { vec.resize(n); }

  std::vector<int>::size_type size() const { return vec.size(); }

  const std::vector<T> &get_vec() const { return vec; }

  virtual std::ostream &print(std::ostream &s) const {
    for (int ii = 0; ii < size(); ++ii)
      s << vec[ii] << ' ';
    return s;
  }

  virtual void sort() { std::sort(vec.begin(), vec.end()); }

  Array1D<T> &operator+=(const Array1D<T> &rhs) {
    MYCMPSIZES(rhs);
    for (int ii = 0; ii < rhs.size(); ++ii)
      vec[ii] += rhs[ii];
    return *this;
  }

  Array1D<T> &operator-=(const Array1D<T> &rhs) {
    MYCMPSIZES(rhs);
    for (int ii = 0; ii < rhs.size(); ++ii)
      vec[ii] -= rhs[ii];
    return *this;
  }

  Array1D<T> &operator*=(const Array1D<T> &rhs) {
    MYCMPSIZES(rhs);
    for (int ii = 0; ii < rhs.size(); ++ii)
      vec[ii] *= rhs[ii];
    return *this;
  }

  Array1D<T> &operator/=(const Array1D<T> &rhs) {
    MYCMPSIZES(rhs);
    MYCMPSIZES(rhs);
    for (int ii = 0; ii < rhs.size(); ++ii)
      vec[ii] /= rhs[ii];
    return *this;
  }

  Array1D<T> &operator+=(T rhs) {
    for (int ii = 0; ii < this->size(); ++ii)
      vec[ii] += rhs;
    return *this;
  }

  Array1D<T> &operator-=(T rhs) {
    for (int ii = 0; ii < this->size(); ++ii)
      vec[ii] -= rhs;
    return *this;
  }

  Array1D<T> &operator*=(T rhs) {
    for (int ii = 0; ii < this->size(); ++ii)
      vec[ii] *= rhs;
    return *this;
  }

  Array1D<T> &operator/=(T rhs) {
    for (int ii = 0; ii < this->size(); ++ii)
      vec[ii] /= rhs;
    return *this;
  }

  ////////////////////////// create extra tmp? //////////////////////////////
  Array1D<T> &operator*(T rhs) {
    for (int ii = 0; ii < this->size(); ++ii)
      vec[ii] *= rhs;
    return *this;
  }

  Array1D<T> &operator/(T rhs) {
    if (T(0) != rhs)
      for (int ii = 0; ii < this->size(); ++ii)
        vec[ii] /= rhs;
    return *this;
  }

  Array1D<T> &operator+(T rhs) {
    for (int ii = 0; ii < this->size(); ++ii)
      vec[ii] += rhs;
    return *this;
  }

  Array1D<T> &operator-(T rhs) {
    for (int ii = 0; ii < this->size(); ++ii)
      vec[ii] -= rhs;
    return *this;
  }
  ////////////////////////// create extra tmp? //////////////////////////////

protected:
  std::vector<T> vec;

  void check_size(const Array1D<T> &arg) const {
    if (vec.size() != arg.size())
      throw std::runtime_error("Operation attempted on different size arrays");
  }

}; // class Array1D
////////////////////////////////// Array1D //////////////////////////////////

// template<class T>
// std::ostream& operator << (std::ostream& s, const Array1D<T>& A) { return
// A.print(s); }

////////////////////////////////// Array2D //////////////////////////////////
//
// A simple 2d array class written for NelderMead/MultiDirSearch.
// Note Array2D as written, using vector of vector, does not guaranteed
// contiguous data. Check out one of the followings if contiguous
// in memory is necessary:
//
//    http://math.nist.gov/tnt/index.html
//    http://www.oonumerics.org/blitz/
//    http://boost.org/libs/multi_array/doc/user.html
//
//
template <typename S, typename T> class Array2D {

  friend std::ostream &operator<<(std::ostream &os, const Array2D<S, T> &a) {
    return a.print(os);
  }

public:
  virtual ~Array2D() {}

  Array2D(int r = 0, int c = 0) : nrow(r), ncol(c), array(r, S(c)) {}

  S &operator[](int arg) { return array[arg]; }

  const S &operator[](int arg) const { return array[arg]; }
  int ncols() const { return ncol; }

  int nrows() const { return nrow; }

  std::ostream &print(std::ostream &os) const {
    for (int ii = 0; ii < nrow; ++ii) {
      os << array[ii][0];
      for (int jj = 1; jj < ncol; ++jj)
        os << ' ' << array[ii][jj];
      if (nrow - 1 != ii)
        os << '\n';
    }
    return os;
  }

private:
  int nrow, ncol;
  Array1D<S> array;

  Array2D &operator=(Array2D const &); // declare but, purposely, not define
  Array2D(Array2D const &);            // declare but, purposely, not define

}; // class Array2D
////////////////////////////////// Array2D //////////////////////////////////

} // namespace

#endif

#ifdef testArray1D
#include "StopWatch.hh"

template <typename T> bool are_same(const T &tmp1, const T &tmp2) {
  int num = tmp1.size();
  for (int ii = 0; ii < num; ++ii)
    if (tmp1[ii] != tmp2[ii]) {
      std::cerr << "@ ii = " << ii << ' ' << tmp1[ii] << " != " << tmp2[ii]
                << '\n';
      return false;
    }
  return true;
}

int compare(const void *a, const void *b) {
  int ai = *((int *)a);
  int bi = *((int *)b);
  if (ai < bi)
    return -1;
  else if (ai > bi)
    return 1;
  else
    return 0;
}

template <typename T> void tst_qsort(int num, T *ptr) {
  qsort(ptr, num, sizeof(T), &compare);
}

template <typename T>
void tst_qsort(int niter, int num, sherpa::Array1D<T> &array) {
  StopWatch stopwatch("testing qsort");
  for (int ii = 0; ii < niter; ++ii) {
    for (int jj = 0; jj < num; ++jj)
      array[jj] = rand();
    tst_qsort(num, &array[0]);
    tst_qsort(num, &array[0]); // check already sorted
  }
}

template <typename T>
void tst_sort(int niter, int num, sherpa::Array1D<T> &array) {
  StopWatch stopwatch("testing std::sort");
  for (int ii = 0; ii < niter; ++ii) {
    for (int jj = 0; jj < num; ++jj)
      array[jj] = rand();
    array.sort();
    array.sort(); // check already sorted
  }
}

void tst_Array1D(int num) {

  StopWatch("tst_Array1D");
  {
    sherpa::Array1D<int> tmp1(num, 1), tmp2(num, 2), tmp3(num, 3);
    tmp1 += tmp2;
    if (false == are_same(tmp1, tmp3))
      std::cerr << "tst_Array1D: += Array1D failed\n";
  }
  {
    sherpa::Array1D<int> tmp1(num, 3), tmp2(num, 2), tmp3(num, 1);
    tmp1 -= tmp2;
    if (false == are_same(tmp1, tmp3))
      std::cerr << "tst_Array1D: -= Array1D failed\n";
  }
  {
    sherpa::Array1D<int> tmp1(num, 3), tmp2(num, 2), tmp3(num, 6);
    tmp1 *= tmp2;
    if (false == are_same(tmp1, tmp3))
      std::cerr << "tst_Array1D: *= Array1D failed\n";
  }
  {
    sherpa::Array1D<int> tmp1(num, 6), tmp2(num, 2), tmp3(num, 3);
    tmp1 /= tmp2;
    if (false == are_same(tmp1, tmp3))
      std::cerr << "tst_Array1D: /= Array1D failed\n";
  }
  {
    int rhs = 3;
    sherpa::Array1D<int> tmp1(num, 3), tmp2(num, 6);
    tmp1 += rhs;
    if (false == are_same(tmp1, tmp2))
      std::cerr << "tst_Array1D: += val failed\n";
  }
  {
    int rhs = 3;
    sherpa::Array1D<int> tmp1(num, 3), tmp2(num, 0);
    tmp1 -= rhs;
    if (false == are_same(tmp1, tmp2))
      std::cerr << "tst_Array1D: -= val failed\n";
  }
  {
    int rhs = 3;
    sherpa::Array1D<int> tmp1(num, 3), tmp2(num, 9);
    tmp1 *= rhs;
    if (false == are_same(tmp1, tmp2))
      std::cerr << "tst_Array1D: *= val failed\n";
  }
  {
    int rhs = 3;
    sherpa::Array1D<int> tmp1(num, 9), tmp2(num, 3);
    tmp1 /= rhs;
    if (false == are_same(tmp1, tmp2))
      std::cerr << "tst_Array1D: /= val failed\n";
  }
  {
    int rhs = 3;
    sherpa::Array1D<int> tmp1(num, 2), tmp2(num, 2), tmp3(num, 6);
    tmp1 = tmp2 * rhs;
    if (false == are_same(tmp1, tmp3))
      std::cerr << "tst_Array1D: * val failed\n";
  }
  {
    int rhs = 3;
    sherpa::Array1D<int> tmp1(num), tmp2(num, 6), tmp3(num, 2);
    tmp1 = tmp2 / rhs;
    if (false == are_same(tmp1, tmp3))
      std::cerr << "tst_Array1D: / val failed\n";
  }
  {
    int rhs = 3;
    sherpa::Array1D<int> tmp1(num), tmp2(num, 2), tmp3(num, 5);
    tmp1 = tmp2 + rhs;
    if (false == are_same(tmp1, tmp3))
      std::cerr << "tst_Array1D: + val failed\n";
  }
  {
    int rhs = 3;
    sherpa::Array1D<int> tmp1(num), tmp2(num, 2), tmp3(num, -1);
    tmp1 = tmp2 - rhs;
    if (false == are_same(tmp1, tmp3))
      std::cerr << "tst_Array1D: - val failed\n";
  }
}

int main(int argc, char *argv[]) {

  int seed = 1357;
  int niter = 10;
  int num = 1000000;

  if (argc == 2)
    num = atoi(argv[1]);
  if (argc == 3) {
    num = atoi(argv[1]);
    niter = atoi(argv[2]);
  }

  sherpa::Array1D<int> tmp1(num), tmp2(num);

  std::cout << "num = " << num << "\tniter = " << niter << '\n';
  srand(seed); // make sure we are always using the same set of numbers
  tst_qsort(niter, num, tmp1);
  srand(seed); // make sure we are always using the same set of numbers
  tst_sort(niter, num, tmp2);
  if (false == are_same(tmp1, tmp2))
    std::cerr << "Failed sort\n";

  tst_Array1D(num);

  return 0;
}

#endif

#ifdef testArray2D

#include "StopWatch.hh"

template <typename T> void initLaplacian(T &uu, int r, int c) {

  for (int ii = 0; ii < r; ++ii) {
    for (int jj = 0; jj < c; ++jj)
      uu[ii][jj] = ii * ii + jj * jj;
  }
}

template <typename S>
void timeme(S &uu, S &laplacian, int r, int c, const char *header) {

  StopWatch stopwatch(header);

  initLaplacian(uu, r, c);

  for (int kk = 0; kk < 10; kk++) {
    for (int ii = 1; ii < r - 1; ii++) {
      for (int jj = 1; jj < c - 1; jj++) {
        laplacian[ii][jj] = -uu[ii - 1][jj] - uu[ii][jj - 1] +
                            4.0 * uu[ii][jj] - uu[ii][jj + 1] - uu[ii + 1][jj];
      }
    }
  }

  double sum = 0.0;
  for (int ii = 1; ii < r - 1; ++ii)
    for (int jj = 1; jj < c - 1; jj++)
      sum += laplacian[ii][jj];

  std::cout << sum << '\t';

}

template <typename T> T ***alloc3d(int nx, int ny, int nz) {
  T ***ptr = new int **[nx];

  for (int x = 0; x < nx; ++x) {
    ptr[x] = new T *[ny];

    for (int y = 0; y < ny; ++y)
      ptr[x][y] = new T[nz];
  }

  return ptr;
}

template <typename T> void del3d(T ***ptr, int nx, int ny) {
  for (int x = 0; x < nx; ++x) {
    for (int y = 0; y < ny; ++y)
      delete[] ptr[x][y], ptr[x][y] = 0;
    delete[] ptr[x], ptr[x] = 0;
  }
  delete[] ptr;
  ptr = 0;
}

template <typename T> T **alloc2d(int r, int c) {
  T **ptr = new T *[r];
  for (int ii = 0; ii < r; ++ii)
    ptr[ii] = new T[c];
  return ptr;
}

template <typename T> void del2d(T **ptr, int r) {

  for (int ii = 0; ii < r; ++ii)
    delete[] ptr[ii];
  delete[] ptr;
}

template <typename T> void timeclassic(int r, int c) {

  T **uu = alloc2d<T>(r, c);
  T **laplacian = alloc2d<T>(r, c);

  timeme(uu, laplacian, r, c, "classic");

  del2d(laplacian, r);
  del2d(uu, r);
}

template <typename T> void timebracket(int r, int c) {

  sherpa::Array2D<sherpa::Array1D<T>, T> uu(r, c), laplacian(r, c);

  for (int ii = 0; ii < r; ++ii) {
    for (int jj = 0; jj < c; ++jj)
      uu[ii][jj] = ii * ii + jj * jj;
  }

  timeme(uu, laplacian, r, c, "sherpa::Array2D[][]");
}

int main(int argc, char *argv[]) {

  const int num = 4000;

  int row = num;
  int col = num;
  if (argc == 3) {
    row = atoi(argv[1]);
    col = atoi(argv[2]);
    std::cout << "row = " << row << "\tcol = " << col << '\n';
  }

  timeclassic<double>(row, col);
  timebracket<double>(row, col);

  return 0;
}

#endif

/*
// To test Array1D:
//
cp myArray.hh tmp.cc; g++ -Wall -ansi -O3 -DtestArray1D -DNDEBUG tmp.cc; rm -f tmp.cc; a.out
//
// num = 1000000   niter = 10      alreadysorted = 1
// qsort took 2.68 secs
// sort took 1 secs
//
// C++ stl sort is faster then C lib qsort!

// To test Array2D:
//
cp myArray.hh tmp.cc; g++ -Wall -ansi -pedantic -O3 -DtestArray2D -DNDEBUG tmp.cc; rm -f tmp.cc; a.out
*/
