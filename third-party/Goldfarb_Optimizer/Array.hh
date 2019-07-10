// $Id: Array.hh 249 2008-11-20 09:58:23Z schaerf $
// This file is part of EasyLocalpp: a C++ Object-Oriented framework
// aimed at easing the development of Local Search algorithms.
// Copyright (C) 2001--2008 Andrea Schaerf, Luca Di Gaspero.
//
// This software may be modified and distributed under the terms
// of the MIT license.  See the LICENSE file for details.

#if !defined(_ARRAY_HH)
#define _ARRAY_HH

#include <set>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>

enum MType { DIAG };
namespace GolDIdnani{

  template <typename T>
  class GVect
  {
  public:
    GVect();
    GVect(const unsigned int n);
    GVect(const T& a, const unsigned int n); //initialize to constant value
    GVect(const T* a, const unsigned int n); // Initialize to array
    GVect(const GVect &rhs); // copy constructor
    ~GVect(); // destructor

    inline void set(const T* a, const unsigned int n);
    GVect<T> extract(const std::set<unsigned int>& indexes) const;
    inline T& operator[](const unsigned int& i); //i-th element
    inline const T& operator[](const unsigned int& i) const;

    inline unsigned int size() const;
    inline void resize(const unsigned int n);
    inline void resize(const T& a, const unsigned int n);

    GVect<T>& operator=(const GVect<T>& rhs); //assignment
    GVect<T>& operator=(const T& a); //assign a to every element
    inline GVect<T>& operator+=(const GVect<T>& rhs);
    inline GVect<T>& operator-=(const GVect<T>& rhs);
    inline GVect<T>& operator*=(const GVect<T>& rhs);
    inline GVect<T>& operator/=(const GVect<T>& rhs);
    inline GVect<T>& operator^=(const GVect<T>& rhs);
    inline GVect<T>& operator+=(const T& a);
    inline GVect<T>& operator-=(const T& a);
    inline GVect<T>& operator*=(const T& a);
    inline GVect<T>& operator/=(const T& a);
    inline GVect<T>& operator^=(const T& a);
  private:
    unsigned int row; // size of array. upper index is n-1
    T* v; // storage for data
  };
  template <typename T>
  inline GVect<T> sqrt(const GVect<T>& v)
  {
    GVect<T> tmp(v.size());
    for (unsigned int i = 0; i < v.size(); i++)
      tmp[i] = sqrt(v[i]);

    return tmp;
  }

  template <typename T>
  GVect<T>::GVect()
    : row(0), v(0)
  {}

  template <typename T>
  GVect<T>::GVect(const unsigned int n)
    : v(new T[n])
  {
    this->row = n;
  }

  template <typename T>
  GVect<T>::GVect(const T& a, const unsigned int n)
    : v(new T[n])
  {
    this->row = n;
    for (unsigned int i = 0; i < n; i++)
      v[i] = a;
  }

  template <typename T>
  GVect<T>::GVect(const T* a, const unsigned int n)
    : v(new T[n])
  {
    this->row = n;
    for (unsigned int i = 0; i < n; i++)
      v[i] = *a++;
  }

  template <typename T>
  GVect<T>::GVect(const GVect<T>& rhs)
    : v(new T[rhs.n])
  {
    this->row = rhs.n;
    for (unsigned int	i = 0; i < row; i++)
      v[i] = rhs[i];
  }

  template <typename T>
  GVect<T>::~GVect()
  {
    if (v != 0)
      delete[] (v);
  }

  template <typename T>
  void GVect<T>::resize(const unsigned int n)
  {
    if (n == this->row)
      return;
    if (v != 0)
      delete[] (v);
    v = new T[n];
    this->row = n;
  }

  template <typename T>
  void GVect<T>::resize(const T& a, const unsigned int n)
  {
    resize(n);
    for (unsigned int i = 0; i < n; i++)
      v[i] = a;
  }


  template <typename T>
  inline GVect<T>& GVect<T>::operator=(const GVect<T>& rhs)
  // postcondition: normal assignment via copying has been performed;
  // if vector and rhs were different sizes, vector
  // has been resized to match the size of rhs
  {
    if (this != &rhs)
      {
        resize(rhs.n);
        for (unsigned int i = 0; i < row; i++)
          v[i] = rhs[i];
      }
    return *this;
  }

  template <typename T>
  inline GVect<T> & GVect<T>::operator=(const T& a) //assign a to every element
  {
    for (unsigned int i = 0; i < row; i++)
      v[i] = a;
    return *this;
  }

  template <typename T>
  inline T & GVect<T>::operator[](const unsigned int& i) //subscripting
  {
    return v[i];
  }

  template <typename T>
  inline const T& GVect<T>::operator[](const unsigned int& i) const //subscripting
  {
    return v[i];
  }

  template <typename T>
  inline unsigned int GVect<T>::size() const
  {
    return row;
  }

  template <typename T>
  inline void GVect<T>::set(const T* a, unsigned int n)
  {
    resize(n);
    for (unsigned int i = 0; i < n; i++)
      v[i] = a[i];
  }

  template <typename T>
  inline GVect<T> GVect<T>::extract(const std::set<unsigned int>& indexes) const
  {
    GVect<T> tmp(indexes.size());
    unsigned int i = 0;

    for (std::set<unsigned int>::const_iterator el = indexes.begin(); el != indexes.end(); el++)
      {
        if (*el >= row)
          throw std::logic_error("Error extracting subvector: the indexes are out of vector bounds");
        tmp[i++] = v[*el];
      }

    return tmp;
  }

  template <typename T>
  inline GVect<T>& GVect<T>::operator+=(const GVect<T>& rhs)
  {
    if (this->size() != rhs.size())
      throw std::logic_error("Operator+=: vectors have different sizes");
    for (unsigned int i = 0; i < row; i++)
      v[i] += rhs[i];

    return *this;
  }


  template <typename T>
  inline GVect<T>& GVect<T>::operator+=(const T& a)
  {
    for (unsigned int i = 0; i < row; i++)
      v[i] += a;

    return *this;
  }

  template <typename T>
  inline GVect<T> operator+(const GVect<T>& rhs)
  {
    return rhs;
  }

  template <typename T>
  inline GVect<T> operator+(const GVect<T>& lhs, const GVect<T>& rhs)
  {
    if (lhs.size() != rhs.size())
      throw std::logic_error("Operator+: vectors have different sizes");
    GVect<T> tmp(lhs.size());
    for (unsigned int i = 0; i < lhs.size(); i++)
      tmp[i] = lhs[i] + rhs[i];

    return tmp;
  }

  template <typename T>
  inline GVect<T> operator+(const GVect<T>& lhs, const T& a)
  {
    GVect<T> tmp(lhs.size());
    for (unsigned int i = 0; i < lhs.size(); i++)
      tmp[i] = lhs[i] + a;

    return tmp;
  }

  template <typename T>
  inline GVect<T> operator+(const T& a, const GVect<T>& rhs)
  {
    GVect<T> tmp(rhs.size());
    for (unsigned int i = 0; i < rhs.size(); i++)
      tmp[i] = a + rhs[i];

    return tmp;
  }

  template <typename T>
  inline GVect<T>& GVect<T>::operator-=(const GVect<T>& rhs)
  {
    if (this->size() != rhs.size())
      throw std::logic_error("Operator-=: vectors have different sizes");
    for (unsigned int i = 0; i < row; i++)
      v[i] -= rhs[i];

    return *this;
  }


  template <typename T>
  inline GVect<T>& GVect<T>::operator-=(const T& a)
  {
    for (unsigned int i = 0; i < row; i++)
      v[i] -= a;

    return *this;
  }

  template <typename T>
  inline GVect<T> operator-(const GVect<T>& rhs)
  {
    return (T)(-1) * rhs;
  }

  template <typename T>
  inline GVect<T> operator-(const GVect<T>& lhs, const GVect<T>& rhs)
  {
    if (lhs.size() != rhs.size())
      throw std::logic_error("Operator-: vectors have different sizes");
    GVect<T> tmp(lhs.size());
    for (unsigned int i = 0; i < lhs.size(); i++)
      tmp[i] = lhs[i] - rhs[i];

    return tmp;
  }

  template <typename T>
  inline GVect<T> operator-(const GVect<T>& lhs, const T& a)
  {
    GVect<T> tmp(lhs.size());
    for (unsigned int i = 0; i < lhs.size(); i++)
      tmp[i] = lhs[i] - a;

    return tmp;
  }

  template <typename T>
  inline GVect<T> operator-(const T& a, const GVect<T>& rhs)
  {
    GVect<T> tmp(rhs.size());
    for (unsigned int i = 0; i < rhs.size(); i++)
      tmp[i] = a - rhs[i];

    return tmp;
  }

  template <typename T>
  inline GVect<T>& GVect<T>::operator*=(const GVect<T>& rhs)
  {
    if (this->size() != rhs.size())
      throw std::logic_error("Operator*=: vectors have different sizes");
    for (unsigned int i = 0; i < row; i++)
      v[i] *= rhs[i];

    return *this;
  }


  template <typename T>
  inline GVect<T>& GVect<T>::operator*=(const T& a)
  {
    for (unsigned int i = 0; i < row; i++)
      v[i] *= a;

    return *this;
  }

  template <typename T>
  inline GVect<T> operator*(const GVect<T>& lhs, const GVect<T>& rhs)
  {
    if (lhs.size() != rhs.size())
      throw std::logic_error("Operator*: vectors have different sizes");
    GVect<T> tmp(lhs.size());
    for (unsigned int i = 0; i < lhs.size(); i++)
      tmp[i] = lhs[i] * rhs[i];

    return tmp;
  }

  template <typename T>
  inline GVect<T> operator*(const GVect<T>& lhs, const T& a)
  {
    GVect<T> tmp(lhs.size());
    for (unsigned int i = 0; i < lhs.size(); i++)
      tmp[i] = lhs[i] * a;

    return tmp;
  }

  template <typename T>
  inline GVect<T> operator*(const T& a, const GVect<T>& rhs)
  {
    GVect<T> tmp(rhs.size());
    for (unsigned int i = 0; i < rhs.size(); i++)
      tmp[i] = a * rhs[i];

    return tmp;
  }

  template <typename T>
  inline GVect<T>& GVect<T>::operator/=(const GVect<T>& rhs)
  {
    if (this->size() != rhs.size())
      throw std::logic_error("Operator/=: vectors have different sizes");
    for (unsigned int i = 0; i < row; i++)
      v[i] /= rhs[i];

    return *this;
  }


  template <typename T>
  inline GVect<T>& GVect<T>::operator/=(const T& a)
  {
    for (unsigned int i = 0; i < row; i++)
      v[i] /= a;

    return *this;
  }

  template <typename T>
  inline GVect<T> operator/(const GVect<T>& lhs, const GVect<T>& rhs)
  {
    if (lhs.size() != rhs.size())
      throw std::logic_error("Operator/: vectors have different sizes");
    GVect<T> tmp(lhs.size());
    for (unsigned int i = 0; i < lhs.size(); i++)
      tmp[i] = lhs[i] / rhs[i];

    return tmp;
  }

  template <typename T>
  inline GVect<T> operator/(const GVect<T>& lhs, const T& a)
  {
    GVect<T> tmp(lhs.size());
    for (unsigned int i = 0; i < lhs.size(); i++)
      tmp[i] = lhs[i] / a;

    return tmp;
  }

  template <typename T>
  inline GVect<T> operator/(const T& a, const GVect<T>& rhs)
  {
    GVect<T> tmp(rhs.size());
    for (unsigned int i = 0; i < rhs.size(); i++)
      tmp[i] = a / rhs[i];

    return tmp;
  }

  template <typename T>
  inline GVect<T> operator^(const GVect<T>& lhs, const GVect<T>& rhs)
  {
    if (lhs.size() != rhs.size())
      throw std::logic_error("Operator^: vectors have different sizes");
    GVect<T> tmp(lhs.size());
    for (unsigned int i = 0; i < lhs.size(); i++)
      tmp[i] = pow(lhs[i], rhs[i]);

    return tmp;
  }

  template <typename T>
  inline GVect<T> operator^(const GVect<T>& lhs, const T& a)
  {
    GVect<T> tmp(lhs.size());
    for (unsigned int i = 0; i < lhs.size(); i++)
      tmp[i] = pow(lhs[i], a);

    return tmp;
  }

  template <typename T>
  inline GVect<T> operator^(const T& a, const GVect<T>& rhs)
  {
    GVect<T> tmp(rhs.size());
    for (unsigned int i = 0; i < rhs.size(); i++)
      tmp[i] = pow(a, rhs[i]);

    return tmp;
  }

  template <typename T>
  inline GVect<T>& GVect<T>::operator^=(const GVect<T>& rhs)
  {
    if (this->size() != rhs.size())
      throw std::logic_error("Operator^=: vectors have different sizes");
    for (unsigned int i = 0; i < row; i++)
      v[i] = pow(v[i], rhs[i]);

    return *this;
  }

  template <typename T>
  inline GVect<T>& GVect<T>::operator^=(const T& a)
  {
    for (unsigned int i = 0; i < row; i++)
      v[i] = pow(v[i], a);

    return *this;
  }

  template <typename T>
  inline bool operator==(const GVect<T>& v, const GVect<T>& w)
  {
    if (v.size() != w.size())
      throw std::logic_error("GVects of different size are not confrontable");
    for (unsigned i = 0; i < v.size(); i++)
      if (v[i] != w[i])
        return false;
    return true;
  }

  template <typename T>
  inline bool operator!=(const GVect<T>& v, const GVect<T>& w)
  {
    if (v.size() != w.size())
      throw std::logic_error("GVects of different size are not confrontable");
    for (unsigned i = 0; i < v.size(); i++)
      if (v[i] != w[i])
        return true;
    return false;
  }

  template <typename T>
  inline bool operator<(const GVect<T>& v, const GVect<T>& w)
  {
    if (v.size() != w.size())
      throw std::logic_error("GVects of different size are not confrontable");
    for (unsigned i = 0; i < v.size(); i++)
      if (v[i] >= w[i])
        return false;
    return true;
  }

  template <typename T>
  inline bool operator<=(const GVect<T>& v, const GVect<T>& w)
  {
    if (v.size() != w.size())
      throw std::logic_error("GVects of different size are not confrontable");
    for (unsigned i = 0; i < v.size(); i++)
      if (v[i] > w[i])
        return false;
    return true;
  }

  template <typename T>
  inline bool operator>(const GVect<T>& v, const GVect<T>& w)
  {
    if (v.size() != w.size())
      throw std::logic_error("GVects of different size are not confrontable");
    for (unsigned i = 0; i < v.size(); i++)
      if (v[i] <= w[i])
        return false;
    return true;
  }

  template <typename T>
  inline bool operator>=(const GVect<T>& v, const GVect<T>& w)
  {
    if (v.size() != w.size())
      throw std::logic_error("GVects of different size are not confrontable");
    for (unsigned i = 0; i < v.size(); i++)
      if (v[i] < w[i])
        return false;
    return true;
  }

  /**
     Input/Output
  */
  template <typename T>
  inline std::ostream& operator<<(std::ostream& os, const GVect<T>& v)
  {
    os << std::endl << v.size() << std::endl;
    for (unsigned int i = 0; i < v.size() - 1; i++)
      os << std::setw(20) << std::setprecision(16) << v[i] << ", ";
    os << std::setw(20) << std::setprecision(16) << v[v.size() - 1] << std::endl;

    return os;
  }

  template <typename T>
  std::istream& operator>>(std::istream& is, GVect<T>& v)
  {
    int elements;
    char comma;
    is >> elements;
    v.resize(elements);
    for (unsigned int i = 0; i < elements; i++)
      is >> v[i] >> comma;

    return is;
  }

  /**
     Index utilities
  */

  std::set<unsigned int> seq(unsigned int s, unsigned int e);

  std::set<unsigned int> singleton(unsigned int i);

  template <typename T>
  class CanonicalBaseGVect : public GVect<T>
  {
  public:
    CanonicalBaseGVect(unsigned int i, unsigned int n);
    inline void reset(unsigned int i);
  private:
    unsigned int e;
  };

  template <typename T>
  CanonicalBaseGVect<T>::CanonicalBaseGVect(unsigned int i, unsigned int n)
    : GVect<T>((T)0, n), e(i)
  { (*this)[e] = (T)1; }

  template <typename T>
  inline void CanonicalBaseGVect<T>::reset(unsigned int i)
  {
    (*this)[e] = (T)0;
    e = i;
    (*this)[e] = (T)1;
  }

#include <stdexcept>

  template <typename T>
  inline T sum(const GVect<T>& v)
  {
    T tmp = (T)0;
    for (unsigned int i = 0; i < v.size(); i++)
      tmp += v[i];

    return tmp;
  }

  template <typename T>
  inline T prod(const GVect<T>& v)
  {
    T tmp = (T)1;
    for (unsigned int i = 0; i < v.size(); i++)
      tmp *= v[i];

    return tmp;
  }

  template <typename T>
  inline T mean(const GVect<T>& v)
  {
    T sum = (T)0;
    for (unsigned int i = 0; i < v.size(); i++)
      sum += v[i];
    return sum / v.size();
  }

  template <typename T>
  inline T median(const GVect<T>& v)
  {
    GVect<T> tmp = sort(v);
    if (v.size() % 2 == 1) // it is an odd-sized vector
      return tmp[v.size() / 2];
    else
      return 0.5 * (tmp[v.size() / 2 - 1] + tmp[v.size() / 2]);
  }

  template <typename T>
  inline T stdev(const GVect<T>& v, bool sample_correction = false)
  {
    return sqrt(var(v, sample_correction));
  }

  template <typename T>
  inline T var(const GVect<T>& v, bool sample_correction = false)
  {
    T sum = (T)0, ssum = (T)0;
    unsigned int n = v.size();
    for (unsigned int i = 0; i < n; i++)
      {
        sum += v[i];
        ssum += (v[i] * v[i]);
      }
    if (!sample_correction)
      return (ssum / n) - (sum / n) * (sum / n);
    else
      return n * ((ssum / n) - (sum / n) * (sum / n)) / (n - 1);
  }

  template <typename T>
  inline T max(const GVect<T>& v)
  {
    T value = v[0];
    for (unsigned int i = 1; i < v.size(); i++)
      value = std::max(v[i], value);

    return value;
  }

  template <typename T>
  inline T min(const GVect<T>& v)
  {
    T value = v[0];
    for (unsigned int i = 1; i < v.size(); i++)
      value = std::min(v[i], value);

    return value;
  }

  template <typename T>
  inline unsigned int index_max(const GVect<T>& v)
  {
    unsigned int max = 0;
    for (unsigned int i = 1; i < v.size(); i++)
      if (v[i] > v[max])
        max = i;

    return max;
  }

  template <typename T>
  inline unsigned int index_min(const GVect<T>& v)
  {
    unsigned int min = 0;
    for (unsigned int i = 1; i < v.size(); i++)
      if (v[i] < v[min])
        min = i;

    return min;
  }


  template <typename T>
  inline T dot_prod(const GVect<T>& a, const GVect<T>& b)
  {
    T sum = (T)0;
    if (a.size() != b.size())
      throw std::logic_error("Dotprod error: the vectors are not the same size");
    for (unsigned int i = 0; i < a.size(); i++)
      sum += a[i] * b[i];

    return sum;
  }

  /**
     Single element mathematical functions
  */

  template <typename T>
  inline GVect<T> exp(const GVect<T>& v)
  {
    GVect<T> tmp(v.size());
    for (unsigned int i = 0; i < v.size(); i++)
      tmp[i] = exp(v[i]);

    return tmp;
  }

  template <typename T>
  inline GVect<T> log(const GVect<T>& v)
  {
    GVect<T> tmp(v.size());
    for (unsigned int i = 0; i < v.size(); i++)
      tmp[i] = log(v[i]);

    return tmp;
  }


  template <typename T>
  inline GVect<T> pow(const GVect<T>& v, double a)
  {
    GVect<T> tmp(v.size());
    for (unsigned int i = 0; i < v.size(); i++)
      tmp[i] = pow(v[i], a);

    return tmp;
  }

  template <typename T>
  inline GVect<T> abs(const GVect<T>& v)
  {
    GVect<T> tmp(v.size());
    for (unsigned int i = 0; i < v.size(); i++)
      tmp[i] = (T)fabs(v[i]);

    return tmp;
  }

  template <typename T>
  inline GVect<T> sign(const GVect<T>& v)
  {
    GVect<T> tmp(v.size());
    for (unsigned int i = 0; i < v.size(); i++)
      tmp[i] = v[i] > 0 ? +1 : v[i] == 0 ? 0 : -1;

    return tmp;
  }

  template <typename T>
  inline unsigned int partition(GVect<T>& v, unsigned int begin, unsigned int end)
  {
    unsigned int i = begin + 1, j = begin + 1;
    T pivot = v[begin];
    while (j <= end)
      {
        if (v[j] < pivot) {
          std::swap(v[i], v[j]);
          i++;
        }
        j++;
      }
    v[begin] = v[i - 1];
    v[i - 1] = pivot;
    return i - 2;
  }


  template <typename T>
  inline void quicksort(GVect<T>& v, unsigned int begin, unsigned int end)
  {
    if (end > begin)
      {
        unsigned int index = partition(v, begin, end);
        quicksort(v, begin, index);
        quicksort(v, index + 2, end);
      }
  }

  template <typename T>
  inline GVect<T> sort(const GVect<T>& v)
  {
    GVect<T> tmp(v);

    quicksort<T>(tmp, 0, tmp.size() - 1);

    return tmp;
  }

  template <typename T>
  inline GVect<double> rank(const GVect<T>& v)
  {
    GVect<T> tmp(v);
    GVect<double> tmp_rank(0.0, v.size());

    for (unsigned int i = 0; i < tmp.size(); i++)
      {
        unsigned int smaller = 0, equal = 0;
        for (unsigned int j = 0; j < tmp.size(); j++)
          if (i == j)
            continue;
          else
            if (tmp[j] < tmp[i])
              smaller++;
            else if (tmp[j] == tmp[i])
              equal++;
        tmp_rank[i] = smaller + 1;
        if (equal > 0)
          {
            for (unsigned int j = 1; j <= equal; j++)
              tmp_rank[i] += smaller + 1 + j;
            tmp_rank[i] /= (double)(equal + 1);
          }
      }

    return tmp_rank;
  }

  //enum MType { DIAG };

  template <typename T>
  class GMatr
  {
  public:
    GMatr(); // Default constructor
    GMatr(const unsigned int n, const unsigned int m); // Construct a n x m matrix
    GMatr(const T& a, const unsigned int n, const unsigned int m); // Initialize the content to constant a
    GMatr(MType t, const T& a, const T& o, const unsigned int n, const unsigned int m);
    GMatr(MType t, const GVect<T>& v, const T& o, const unsigned int n, const unsigned int m);
    GMatr(const T* a, const unsigned int n, const unsigned int m); // Initialize to array
    GMatr(const GMatr<T>& rhs); // Copy constructor
    ~GMatr(); // destructor

    inline T* operator[](const unsigned int& i) { return v[i]; } // Subscripting: row i
    inline const T* operator[](const unsigned int& i) const { return v[i]; }; // const subsctipting

    inline void resize(const unsigned int n, const unsigned int m);
    inline void resize(const T& a, const unsigned int n, const unsigned int m);


    inline GVect<T> extractRow(const unsigned int i) const;
    inline GVect<T> extractColumn(const unsigned int j) const;
    inline GVect<T> extractDiag() const;
    inline GMatr<T> extractRows(const std::set<unsigned int>& indexes) const;
    inline GMatr<T> extractColumns(const std::set<unsigned int>& indexes) const;
    inline GMatr<T> extract(const std::set<unsigned int>& r_indexes, const std::set<unsigned int>& c_indexes) const;

    inline void set(const T* a, unsigned int n, unsigned int m);
    inline void set(const std::set<unsigned int>& r_indexes, const std::set<unsigned int>& c_indexes, const GMatr<T>& m);
    inline void setRow(const unsigned int index, const GVect<T>& v);
    inline void setRow(const unsigned int index, const GMatr<T>& v);
    inline void setRows(const std::set<unsigned int>& indexes, const GMatr<T>& m);
    inline void setColumn(const unsigned int index, const GVect<T>& v);
    inline void setColumn(const unsigned int index, const GMatr<T>& v);
    inline void setColumns(const std::set<unsigned int>& indexes, const GMatr<T>& m);


    inline unsigned int nrows() const { return row; } // number of rows
    inline unsigned int ncols() const { return col; } // number of columns

    inline GMatr<T>& operator=(const GMatr<T>& rhs); // Assignment operator
    inline GMatr<T>& operator=(const T& a); // Assign to every element value a
    inline GMatr<T>& operator+=(const GMatr<T>& rhs);
    inline GMatr<T>& operator-=(const GMatr<T>& rhs);
    inline GMatr<T>& operator*=(const GMatr<T>& rhs);
    inline GMatr<T>& operator/=(const GMatr<T>& rhs);
    inline GMatr<T>& operator^=(const GMatr<T>& rhs);
    inline GMatr<T>& operator+=(const T& a);
    inline GMatr<T>& operator-=(const T& a);
    inline GMatr<T>& operator*=(const T& a);
    inline GMatr<T>& operator/=(const T& a);
    inline GMatr<T>& operator^=(const T& a);
    inline operator GVect<T>();
  private:
    unsigned int row; // number of rows
    unsigned int col; // number of columns
    T **v; // storage for data
  };

  template <typename T>
  GMatr<T>::GMatr()
    : row(0), col(0), v(0)
  {}

  template <typename T>
  GMatr<T>::GMatr(unsigned int n, unsigned int m)
    : v(new T*[n])
  {
    unsigned int i;
    this->row = n; this->col = m;
    v[0] = new T[m * n];
    for (i = 1; i < n; i++)
      v[i] = v[i - 1] + m;
  }

  template <typename T>
  GMatr<T>::GMatr(const T& a, unsigned int n, unsigned int m)
    : v(new T*[n])
  {
    unsigned int i, j;
    this->row = n; this->col = m;
    v[0] = new T[m * n];
    for (i = 1; i < n; i++)
      v[i] = v[i - 1] + m;
    for (i = 0; i < n; i++)
      for (j = 0; j < m; j++)
        v[i][j] = a;
  }

  template <class T>
  GMatr<T>::GMatr(const T* a, unsigned int n, unsigned int m)
    : v(new T*[n])
  {
    unsigned int i, j;
    this->row = n; this->col = m;
    v[0] = new T[m * n];
    for (i = 1; i < n; i++)
      v[i] = v[i - 1] + m;
    for (i = 0; i < n; i++)
      for (j = 0; j < m; j++)
        v[i][j] = *a++;
  }

  template <class T>
  GMatr<T>::GMatr(MType t, const T& a, const T& o, unsigned int n, unsigned int m)
    : v(new T*[n])
  {
    unsigned int i, j;
    this->row = n; this->col = m;
    v[0] = new T[m * n];
    for (i = 1; i < n; i++)
      v[i] = v[i - 1] + m;
    switch (t)
      {
      case DIAG:
        for (i = 0; i < n; i++)
          for (j = 0; j < m; j++)
            if (i != j)
              v[i][j] = o;
            else
              v[i][j] = a;
        break;
      default:
        throw std::logic_error("GMatr type not supported");
      }
  }

  template <class T>
  GMatr<T>::GMatr(MType t, const GVect<T>& a, const T& o, unsigned int n, unsigned int m)
    : v(new T*[n])
  {
    unsigned int i, j;
    this->row = n; this->col = m;
    v[0] = new T[m * n];
    for (i = 1; i < n; i++)
      v[i] = v[i - 1] + m;
    switch (t)
      {
      case DIAG:
        for (i = 0; i < n; i++)
          for (j = 0; j < m; j++)
            if (i != j)
              v[i][j] = o;
            else
              v[i][j] = a[i];
        break;
      default:
        throw std::logic_error("GMatr type not supported");
      }
  }

  template <typename T>
  GMatr<T>::GMatr(const GMatr<T>& rhs)
    : v(new T*[rhs.n])
  {
    unsigned int i, j;
    row = rhs.row; col = rhs.col;
    v[0] = new T[col * row];
    for (i = 1; i < row; i++)
      v[i] = v[i - 1] + col;
    for (i = 0; i < row; i++)
      for (j = 0; j < col; j++)
        v[i][j] = rhs[i][j];
  }

  template <typename T>
  GMatr<T>::~GMatr()
  {
    if (v != 0) {
      delete[] (v[0]);
      delete[] (v);
    }
  }

  template <typename T>
  inline GMatr<T>& GMatr<T>::operator=(const GMatr<T> &rhs)
  // postcondition: normal assignment via copying has been performed;
  // if matrix and rhs were different sizes, matrix
  // has been resized to match the size of rhs
  {
    unsigned int i, j;
    if (this != &rhs)
      {
        resize(rhs.n, rhs.m);
        for (i = 0; i < row; i++)
          for (j = 0; j < col; j++)
            v[i][j] = rhs[i][j];
      }
    return *this;
  }

  template <typename T>
  inline GMatr<T>& GMatr<T>::operator=(const T& a) // assign a to every element
  {
    unsigned int i, j;
    for (i = 0; i < row; i++)
      for (j = 0; j < col; j++)
        v[i][j] = a;
    return *this;
  }


  template <typename T>
  inline void GMatr<T>::resize(const unsigned int n, const unsigned int m)
  {
    unsigned int i;
    if (n == this->row && m == this->col)
      return;
    if (v != 0)
      {
        delete[] (v[0]);
        delete[] (v);
      }
    this->row = n; this->col = m;
    v = new T*[n];
    v[0] = new T[m * n];
    for (i = 1; i < n; i++)
      v[i] = v[i - 1] + m;
  }

  template <typename T>
  inline void GMatr<T>::resize(const T& a, const unsigned int n, const unsigned int m)
  {
    unsigned int i, j;
    resize(n, m);
    for (i = 0; i < n; i++)
      for (j = 0; j < m; j++)
        v[i][j] = a;
  }



  template <typename T>
  inline GVect<T> GMatr<T>::extractRow(const unsigned int i) const
  {
    if (i >= row)
      throw std::logic_error("Error in extractRow: trying to extract a row out of matrix bounds");
    GVect<T> tmp(v[i], col);

    return tmp;
  }

  template <typename T>
  inline GVect<T> GMatr<T>::extractColumn(const unsigned int j) const
  {
    unsigned int i;
    if (j >= col)
      throw std::logic_error("Error in extractRow: trying to extract a row out of matrix bounds");
    GVect<T> tmp(row);

    for (i = 0; i < row; i++)
      tmp[i] = v[i][j];

    return tmp;
  }

  template <typename T>
  inline GVect<T> GMatr<T>::extractDiag() const
  {
    unsigned int d = std::min(row, col), i;

    GVect<T> tmp(d);

    for (i = 0; i < d; i++)
      tmp[i] = v[i][i];

    return tmp;

  }

  template <typename T>
  inline GMatr<T> GMatr<T>::extractRows(const std::set<unsigned int>& indexes) const
  {
    GMatr<T> tmp(indexes.size(), col);
    unsigned int i = 0, j;

    for (std::set<unsigned int>::const_iterator el = indexes.begin(); el != indexes.end(); el++)
      {
        for (j = 0; j < col; j++)
          {
            if (*el >= row)
              throw std::logic_error("Error extracting rows: the indexes are out of matrix bounds");
            tmp[i][j] = v[*el][j];
          }
        i++;
      }

    return tmp;
  }

  template <typename T>
  inline GMatr<T> GMatr<T>::extractColumns(const std::set<unsigned int>& indexes) const
  {
    GMatr<T> tmp(row, indexes.size());
    unsigned int i, j = 0;

    for (std::set<unsigned int>::const_iterator el = indexes.begin(); el != indexes.end(); el++)
      {
        for (i = 0; i < row; i++)
          {
            if (*el >= col)
              throw std::logic_error("Error extracting columns: the indexes are out of matrix bounds");
            tmp[i][j] = v[i][*el];
          }
        j++;
      }

    return tmp;
  }

  template <typename T>
  inline GMatr<T> GMatr<T>::extract(const std::set<unsigned int>& r_indexes, const std::set<unsigned int>& c_indexes) const
  {
    GMatr<T> tmp(r_indexes.size(), c_indexes.size());
    unsigned int i = 0, j;

    for (std::set<unsigned int>::const_iterator r_el = r_indexes.begin(); r_el != r_indexes.end(); r_el++)
      {
        if (*r_el >= row)
          throw std::logic_error("Error extracting submatrix: the indexes are out of matrix bounds");
        j = 0;
        for (std::set<unsigned int>::const_iterator c_el = c_indexes.begin(); c_el != c_indexes.end(); c_el++)
          {
            if (*c_el >= col)
              throw std::logic_error("Error extracting rows: the indexes are out of matrix bounds");
            tmp[i][j] = v[*r_el][*c_el];
            j++;
          }
        i++;
      }

    return tmp;
  }

  template <typename T>
  inline void GMatr<T>::setRow(unsigned int i, const GVect<T>& a)
  {
    if (i >= row)
      throw std::logic_error("Error in setRow: trying to set a row out of matrix bounds");
    if (this->col != a.size())
      throw std::logic_error("Error setting matrix row: ranges are not compatible");
    for (unsigned int j = 0; j < ncols(); j++)
      v[i][j] = a[j];
  }

  template <typename T>
  inline void GMatr<T>::setRow(unsigned int i, const GMatr<T>& a)
  {
    if (i >= row)
      throw std::logic_error("Error in setRow: trying to set a row out of matrix bounds");
    if (this->col != a.ncols())
      throw std::logic_error("Error setting matrix column: ranges are not compatible");
    if (a.nrows() != 1)
      throw std::logic_error("Error setting matrix column with a non-row matrix");
    for (unsigned int j = 0; j < ncols(); j++)
      v[i][j] = a[0][j];
  }

  template <typename T>
  inline void GMatr<T>::setRows(const std::set<unsigned int>& indexes, const GMatr<T>& m)
  {
    unsigned int i = 0;

    if (indexes.size() != m.nrows() || this->col != m.ncols())
      throw std::logic_error("Error setting matrix rows: ranges are not compatible");
    for (std::set<unsigned int>::const_iterator el = indexes.begin(); el != indexes.end(); el++)
      {
        for (unsigned int j = 0; j < ncols(); j++)
          {
            if (*el >= row)
              throw std::logic_error("Error in setRows: trying to set a row out of matrix bounds");
            v[*el][j] = m[i][j];
          }
        i++;
      }
  }

  template <typename T>
  inline void GMatr<T>::setColumn(unsigned int j, const GVect<T>& a)
  {
    if (j >= col)
      throw std::logic_error("Error in setColumn: trying to set a column out of matrix bounds");
    if (this->row != a.size())
      throw std::logic_error("Error setting matrix column: ranges are not compatible");
    for (unsigned int i = 0; i < nrows(); i++)
      v[i][j] = a[i];
  }

  template <typename T>
  inline void GMatr<T>::setColumn(unsigned int j, const GMatr<T>& a)
  {
    if (j >= col)
      throw std::logic_error("Error in setColumn: trying to set a column out of matrix bounds");
    if (this->row != a.nrows())
      throw std::logic_error("Error setting matrix column: ranges are not compatible");
    if (a.ncols() != 1)
      throw std::logic_error("Error setting matrix column with a non-column matrix");
    for (unsigned int i = 0; i < nrows(); i++)
      v[i][j] = a[i][0];
  }


  template <typename T>
  inline void GMatr<T>::setColumns(const std::set<unsigned int>& indexes, const GMatr<T>& a)
  {
    unsigned int j = 0;

    if (indexes.size() != a.ncols() || this->row != a.nrows())
      throw std::logic_error("Error setting matrix columns: ranges are not compatible");
    for (std::set<unsigned int>::const_iterator el = indexes.begin(); el != indexes.end(); el++)
      {
        for (unsigned int i = 0; i < nrows(); i++)
          {
            if (*el >= col)
              throw std::logic_error("Error in setColumns: trying to set a column out of matrix bounds");
            v[i][*el] = a[i][j];
          }
        j++;
      }
  }

  template <typename T>
  inline void GMatr<T>::set(const std::set<unsigned int>& r_indexes, const std::set<unsigned int>& c_indexes, const GMatr<T>& a)
  {
    unsigned int i = 0, j;
    if (c_indexes.size() != a.ncols() || r_indexes.size() != a.nrows())
      throw std::logic_error("Error setting matrix elements: ranges are not compatible");

    for (std::set<unsigned int>::const_iterator r_el = r_indexes.begin(); r_el != r_indexes.end(); r_el++)
      {
        if (*r_el >= row)
          throw std::logic_error("Error in set: trying to set a row out of matrix bounds");
        j = 0;
        for (std::set<unsigned int>::const_iterator c_el = c_indexes.begin(); c_el != c_indexes.end(); c_el++)
          {
            if (*c_el >= col)
              throw std::logic_error("Error in set: trying to set a column out of matrix bounds");
            v[*r_el][*c_el] = a[i][j];
            j++;
          }
        i++;
      }
  }

  template <typename T>
  inline void GMatr<T>::set(const T* a, unsigned int n, unsigned int m)
  {
    if (this->row != n || this->col != m)
      resize(n, m);
    unsigned int k = 0;
    for (unsigned int i = 0; i < n; i++)
      for (unsigned int j = 0; j < m; j++)
        v[i][j] = a[k++];
  }


  template <typename T>
  GMatr<T> operator+(const GMatr<T>& rhs)
  {
    return rhs;
  }

  template <typename T>
  GMatr<T> operator+(const GMatr<T>& lhs, const GMatr<T>& rhs)
  {
    if (lhs.ncols() != rhs.ncols() || lhs.nrows() != rhs.nrows())
      throw std::logic_error("Operator+: matrices have different sizes");
    GMatr<T> tmp(lhs.nrows(), lhs.ncols());
    for (unsigned int i = 0; i < lhs.nrows(); i++)
      for (unsigned int j = 0; j < lhs.ncols(); j++)
        tmp[i][j] = lhs[i][j] + rhs[i][j];

    return tmp;
  }

  template <typename T>
  GMatr<T> operator+(const GMatr<T>& lhs, const T& a)
  {
    GMatr<T> tmp(lhs.nrows(), lhs.ncols());
    for (unsigned int i = 0; i < lhs.nrows(); i++)
      for (unsigned int j = 0; j < lhs.ncols(); j++)
        tmp[i][j] = lhs[i][j] + a;

    return tmp;
  }

  template <typename T>
  GMatr<T> operator+(const T& a, const GMatr<T>& rhs)
  {
    GMatr<T> tmp(rhs.nrows(), rhs.ncols());
    for (unsigned int i = 0; i < rhs.nrows(); i++)
      for (unsigned int j = 0; j < rhs.ncols(); j++)
        tmp[i][j] = a + rhs[i][j];

    return tmp;
  }

  template <typename T>
  inline GMatr<T>& GMatr<T>::operator+=(const GMatr<T>& rhs)
  {
    if (col != rhs.ncols() || row != rhs.nrows())
      throw std::logic_error("Operator+=: matrices have different sizes");
    for (unsigned int i = 0; i < row; i++)
      for (unsigned int j = 0; j < col; j++)
        v[i][j] += rhs[i][j];

    return *this;
  }

  template <typename T>
  inline GMatr<T>& GMatr<T>::operator+=(const T& a)
  {
    for (unsigned int i = 0; i < row; i++)
      for (unsigned int j = 0; j < col; j++)
        v[i][j] += a;

    return *this;
  }

  template <typename T>
  GMatr<T> operator-(const GMatr<T>& rhs)
  {
    return (T)(-1) * rhs;
  }

  template <typename T>
  GMatr<T> operator-(const GMatr<T>& lhs, const GMatr<T>& rhs)
  {
    if (lhs.ncols() != rhs.ncols() || lhs.nrows() != rhs.nrows())
      throw std::logic_error("Operator-: matrices have different sizes");
    GMatr<T> tmp(lhs.nrows(), lhs.ncols());
    for (unsigned int i = 0; i < lhs.nrows(); i++)
      for (unsigned int j = 0; j < lhs.ncols(); j++)
        tmp[i][j] = lhs[i][j] - rhs[i][j];

    return tmp;
  }

  template <typename T>
  GMatr<T> operator-(const GMatr<T>& lhs, const T& a)
  {
    GMatr<T> tmp(lhs.nrows(), lhs.ncols());
    for (unsigned int i = 0; i < lhs.nrows(); i++)
      for (unsigned int j = 0; j < lhs.ncols(); j++)
        tmp[i][j] = lhs[i][j] - a;

    return tmp;
  }

  template <typename T>
  GMatr<T> operator-(const T& a, const GMatr<T>& rhs)
  {
    GMatr<T> tmp(rhs.nrows(), rhs.ncols());
    for (unsigned int i = 0; i < rhs.nrows(); i++)
      for (unsigned int j = 0; j < rhs.ncols(); j++)
        tmp[i][j] = a - rhs[i][j];

    return tmp;
  }

  template <typename T>
  inline GMatr<T>& GMatr<T>::operator-=(const GMatr<T>& rhs)
  {
    if (col != rhs.ncols() || row != rhs.nrows())
      throw std::logic_error("Operator-=: matrices have different sizes");
    for (unsigned int i = 0; i < row; i++)
      for (unsigned int j = 0; j < col; j++)
        v[i][j] -= rhs[i][j];

    return *this;
  }

  template <typename T>
  inline GMatr<T>& GMatr<T>::operator-=(const T& a)
  {
    for (unsigned int i = 0; i < row; i++)
      for (unsigned int j = 0; j < col; j++)
        v[i][j] -= a;

    return *this;
  }

  template <typename T>
  GMatr<T> operator*(const GMatr<T>& lhs, const GMatr<T>& rhs)
  {
    if (lhs.ncols() != rhs.ncols() || lhs.nrows() != rhs.nrows())
      throw std::logic_error("Operator*: matrices have different sizes");
    GMatr<T> tmp(lhs.nrows(), lhs.ncols());
    for (unsigned int i = 0; i < lhs.nrows(); i++)
      for (unsigned int j = 0; j < lhs.ncols(); j++)
        tmp[i][j] = lhs[i][j] * rhs[i][j];

    return tmp;
  }

  template <typename T>
  GMatr<T> operator*(const GMatr<T>& lhs, const T& a)
  {
    GMatr<T> tmp(lhs.nrows(), lhs.ncols());
    for (unsigned int i = 0; i < lhs.nrows(); i++)
      for (unsigned int j = 0; j < lhs.ncols(); j++)
        tmp[i][j] = lhs[i][j] * a;

    return tmp;
  }

  template <typename T>
  GMatr<T> operator*(const T& a, const GMatr<T>& rhs)
  {
    GMatr<T> tmp(rhs.nrows(), rhs.ncols());
    for (unsigned int i = 0; i < rhs.nrows(); i++)
      for (unsigned int j = 0; j < rhs.ncols(); j++)
        tmp[i][j] = a * rhs[i][j];

    return tmp;
  }

  template <typename T>
  inline GMatr<T>& GMatr<T>::operator*=(const GMatr<T>& rhs)
  {
    if (col != rhs.ncols() || row != rhs.nrows())
      throw std::logic_error("Operator*=: matrices have different sizes");
    for (unsigned int i = 0; i < row; i++)
      for (unsigned int j = 0; j < col; j++)
        v[i][j] *= rhs[i][j];

    return *this;
  }

  template <typename T>
  inline GMatr<T>& GMatr<T>::operator*=(const T& a)
  {
    for (unsigned int i = 0; i < row; i++)
      for (unsigned int j = 0; j < col; j++)
        v[i][j] *= a;

    return *this;
  }

  template <typename T>
  GMatr<T> operator/(const GMatr<T>& lhs, const GMatr<T>& rhs)
  {
    if (lhs.ncols() != rhs.ncols() || lhs.nrows() != rhs.nrows())
      throw std::logic_error("Operator+: matrices have different sizes");
    GMatr<T> tmp(lhs.nrows(), lhs.ncols());
    for (unsigned int i = 0; i < lhs.nrows(); i++)
      for (unsigned int j = 0; j < lhs.ncols(); j++)
        tmp[i][j] = lhs[i][j] / rhs[i][j];

    return tmp;
  }

  template <typename T>
  GMatr<T> operator/(const GMatr<T>& lhs, const T& a)
  {
    GMatr<T> tmp(lhs.nrows(), lhs.ncols());
    for (unsigned int i = 0; i < lhs.nrows(); i++)
      for (unsigned int j = 0; j < lhs.ncols(); j++)
        tmp[i][j] = lhs[i][j] / a;

    return tmp;
  }

  template <typename T>
  GMatr<T> operator/(const T& a, const GMatr<T>& rhs)
  {
    GMatr<T> tmp(rhs.nrows(), rhs.ncols());
    for (unsigned int i = 0; i < rhs.nrows(); i++)
      for (unsigned int j = 0; j < rhs.ncols(); j++)
        tmp[i][j] = a / rhs[i][j];

    return tmp;
  }

  template <typename T>
  inline GMatr<T>& GMatr<T>::operator/=(const GMatr<T>& rhs)
  {
    if (col != rhs.ncols() || row != rhs.nrows())
      throw std::logic_error("Operator+=: matrices have different sizes");
    for (unsigned int i = 0; i < row; i++)
      for (unsigned int j = 0; j < col; j++)
        v[i][j] /= rhs[i][j];

    return *this;
  }

  template <typename T>
  inline GMatr<T>& GMatr<T>::operator/=(const T& a)
  {
    for (unsigned int i = 0; i < row; i++)
      for (unsigned int j = 0; j < col; j++)
        v[i][j] /= a;

    return *this;
  }

  template <typename T>
  GMatr<T> operator^(const GMatr<T>& lhs, const T& a)
  {
    GMatr<T> tmp(lhs.nrows(), lhs.ncols());
    for (unsigned int i = 0; i < lhs.nrows(); i++)
      for (unsigned int j = 0; j < lhs.ncols(); j++)
        tmp[i][j] = pow(lhs[i][j], a);

    return tmp;
  }

  template <typename T>
  inline GMatr<T>& GMatr<T>::operator^=(const GMatr<T>& rhs)
  {
    if (col != rhs.ncols() || row != rhs.nrows())
      throw std::logic_error("Operator^=: matrices have different sizes");
    for (unsigned int i = 0; i < row; i++)
      for (unsigned int j = 0; j < col; j++)
        v[i][j] = pow(v[i][j], rhs[i][j]);

    return *this;
  }


  template <typename T>
  inline GMatr<T>& GMatr<T>::operator^=(const T& a)
  {
    for (unsigned int i = 0; i < row; i++)
      for (unsigned int j = 0; j < col; j++)
        v[i][j] = pow(v[i][j], a);

    return *this;
  }

  template <typename T>
  inline GMatr<T>::operator GVect<T>()
  {
    if (row > 1 && col > 1)
      throw std::logic_error("Error matrix cast to vector: trying to cast a multi-dimensional matrix");
    if (row == 1)
      return extractRow(0);
    else
      return extractColumn(0);
  }

  template <typename T>
  inline bool operator==(const GMatr<T>& a, const GMatr<T>& b)
  {
    if (a.nrows() != b.nrows() || a.ncols() != b.ncols())
      throw std::logic_error("Matrices of different size are not confrontable");
    for (unsigned i = 0; i < a.nrows(); i++)
      for (unsigned j = 0; j < a.ncols(); j++)
        if (a[i][j] != b[i][j])
          return false;
    return true;
  }

  template <typename T>
  inline bool operator!=(const GMatr<T>& a, const GMatr<T>& b)
  {
    if (a.nrows() != b.nrows() || a.ncols() != b.ncols())
      throw std::logic_error("Matrices of different size are not confrontable");
    for (unsigned i = 0; i < a.nrows(); i++)
      for (unsigned j = 0; j < a.ncols(); j++)
        if (a[i][j] != b[i][j])
          return true;
    return false;
  }



  /**
     Input/Output
  */
  template <typename T>
  std::ostream& operator<<(std::ostream& os, const GMatr<T>& m)
  {
    os << std::endl << m.nrows() << " " << m.ncols() << std::endl;
    for (unsigned int i = 0; i < m.nrows(); i++)
      {
        for (unsigned int j = 0; j < m.ncols() - 1; j++)
          os << std::setw(20) << std::setprecision(16) << m[i][j] << ", ";
        os << std::setw(20) << std::setprecision(16) << m[i][m.ncols() - 1] << std::endl;
      }

    return os;
  }

  template <typename T>
  std::istream& operator>>(std::istream& is, GMatr<T>& m)
  {
    int rows, cols;
    char comma;
    is >> rows >> cols;
    m.resize(rows, cols);
    for (unsigned int i = 0; i < rows; i++)
      for (unsigned int j = 0; j < cols; j++)
        is >> m[i][j] >> comma;

    return is;
  }

  template <typename T>
  T sign(const T& v)
  {
    if (v >= (T)0.0)
      return (T)1.0;
    else
      return (T)-1.0;
  }

  template <typename T>
  T dist(const T& a, const T& b)
  {
    T abs_a = (T)fabs(a), abs_b = (T)fabs(b);
    if (abs_a > abs_b)
      return abs_a * sqrt((T)1.0 + (abs_b / abs_a) * (abs_b / abs_a));
    else
      return (abs_b == (T)0.0 ? (T)0.0 : abs_b * sqrt((T)1.0 + (abs_a / abs_b) * (abs_a / abs_b)));
  }

  template <typename T>
  void svd(const GMatr<T>& A, GMatr<T>& U, GVect<T>& W, GMatr<T>& V)
  {
    int m = A.nrows(), n = A.ncols(), i, j, k, l, jj, nm;
    const unsigned int max_its = 30;
    bool flag;
    GVect<T> rv1(n);
    U = A;
    W.resize(n);
    V.resize(n, n);
    T anorm, c, f, g, h, s, scale, x, y, z;
    g = scale = anorm = (T)0.0;

    // Householder reduction to bidiagonal form
    for (i = 0; i < n; i++)
      {
        l = i + 1;
        rv1[i] = scale * g;
        g = s = scale = (T)0.0;
        if (i < m)
          {
            for (k = i; k < m; k++)
              scale += fabs(U[k][i]);
            if (scale != (T)0.0)
              {
                for (k = i; k < m; k++)
                  {
                    U[k][i] /= scale;
                    s += U[k][i] * U[k][i];
                  }
                f = U[i][i];
                g = -sign(f) * sqrt(s);
                h = f * g - s;
                U[i][i] = f - g;
                for (j = l; j < n; j++)
                  {
                    s = (T)0.0;
                    for (k = i; k < m; k++)
                      s += U[k][i] * U[k][j];
                    f = s / h;
                    for (k = i; k < m; k++)
                      U[k][j] += f * U[k][i];
                  }
                for (k = i; k < m; k++)
                  U[k][i] *= scale;
              }
          }
        W[i] = scale * g;
        g = s = scale = (T)0.0;
        if (i < m && i != n - 1)
          {
            for (k = l; k < n; k++)
              scale += fabs(U[i][k]);
            if (scale != (T)0.0)
              {
                for (k = l; k < n; k++)
                  {
                    U[i][k] /= scale;
                    s += U[i][k] * U[i][k];
                  }
                f = U[i][l];
                g = -sign(f) * sqrt(s);
                h = f * g - s;
                U[i][l] = f - g;
                for (k = l; k <n; k++)
                  rv1[k] = U[i][k] / h;
                for (j = l; j < m; j++)
                  {
                    s = (T)0.0;
                    for (k = l; k < n; k++)
                      s += U[j][k] * U[i][k];
                    for (k = l; k < n; k++)
                      U[j][k] += s * rv1[k];
                  }
                for (k = l; k < n; k++)
                  U[i][k] *= scale;
              }
          }
        anorm = std::max(anorm, fabs(W[i]) + fabs(rv1[i]));
      }
    // Accumulation of right-hand transformations
    for (i = n - 1; i >= 0; i--)
      {
        if (i < n - 1)
          {
            if (g != (T)0.0)
              {
                for (j = l; j < n; j++)
                  V[j][i] = (U[i][j] / U[i][l]) / g;
                for (j = l; j < n; j++)
                  {
                    s = (T)0.0;
                    for (k = l; k < n; k++)
                      s += U[i][k] * V[k][j];
                    for (k = l; k < n; k++)
                      V[k][j] += s * V[k][i];
                  }
              }
            for (j = l; j < n; j++)
              V[i][j] = V[j][i] = (T)0.0;
          }
        V[i][i] = (T)1.0;
        g = rv1[i];
        l = i;
      }
    // Accumulation of left-hand transformations
    for (i = std::min(m, n) - 1; i >= 0; i--)
      {
        l = i + 1;
        g = W[i];
        for (j = l; j < n; j++)
          U[i][j] = (T)0.0;
        if (g != (T)0.0)
          {
            g = (T)1.0 / g;
            for (j = l; j < n; j++)
              {
                s = (T)0.0;
                for (k = l; k < m; k++)
                  s += U[k][i] * U[k][j];
                f = (s / U[i][i]) * g;
                for (k = i; k < m; k++)
                  U[k][j] += f * U[k][i];
              }
            for (j = i; j < m; j++)
              U[j][i] *= g;
          }
        else
          for (j = i; j < m; j++)
            U[j][i] = (T)0.0;
        U[i][i]++;
      }
    // Diagonalization of the bidiagonal form: loop over singular values, and over allowed iterations.
    for (k = n - 1; k >= 0; k--)
      {
        for (unsigned int its = 0; its < max_its; its++)
          {
            flag = true;
            for (l = k; l >= 0; l--) // FIXME: in NR it was l >= 1 but there subscripts start from one
              { // Test for splitting
                nm = l - 1; // Note that rV[0] is always zero
                if ((T)(fabs(rv1[l]) + anorm) == anorm)
                  {
                    flag = false;
                    break;
                  }
                if ((T)(fabs(W[nm]) + anorm) == anorm)
                  break;
              }
            if (flag)
              {
                // Cancellation of rv1[l], if l > 0 FIXME: it was l > 1 in NR
                c = (T)0.0;
                s = (T)1.0;
                for (i = l; i <= k; i++)
                  {
                    f = s * rv1[i];
                    rv1[i] *= c;
                    if ((T)(fabs(f) + anorm) == anorm)
                      break;
                    g = W[i];
                    h = dist(f, g);
                    W[i] = h;
                    h = (T)1.0 / h;
                    c = g * h;
                    s = -f * h;
                    for (j = 0; j < m; j++)
                      {
                        y = U[j][nm];
                        z = U[j][i];
                        U[j][nm] = y * c + z * s;
                        U[j][i] = z * c - y * s;
                      }
                  }
              }
            z = W[k];
            if (l == k)
              {  // Convergence
                if (z < (T)0.0)
                  { // Singular value is made nonnegative
                    W[k] = -z;
                    for (j = 0; j < n; j++)
                      V[j][k] = -V[j][k];
                  }
                break;
              }
            if (its == max_its)
              throw std::logic_error("Error svd: no convergence in the maximum number of iterations");
            x = W[l];
            nm = k - 1;
            y = W[nm];
            g = rv1[nm];
            h = rv1[k];
            f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
            g = dist(f, (T)1.0);
            f = ((x - z) * (x + z) + h * ((y / (f + sign(f)*fabs(g))) - h)) / x;
            c = s = (T)1.0; // Next QR transformation
            for (j = l; j <= nm; j++)
              {
                i = j + 1;
                g = rv1[i];
                y = W[i];
                h = s * g;
                g *= c;
                z = dist(f, h);
                rv1[j] = z;
                c = f / z;
                s = h / z;
                f = x * c + g * s;
                g = g * c - x * s;
                h = y * s;
                y *= c;
                for (jj = 0; jj < n; jj++)
                  {
                    x = V[jj][j];
                    z = V[jj][i];
                    V[jj][j] = x * c + z * s;
                    V[jj][i] = z * c - x * s;
                  }
                z = dist(f, h);
                W[j] = z;
                if (z != 0) // Rotation can be arbitrary if z = 0
                  {
                    z = (T)1.0 / z;
                    c = f * z;
                    s = h * z;
                  }
                f = c * g + s * y;
                x = c * y - s * g;
                for (jj = 0; jj < m; jj++)
                  {
                    y = U[jj][j];
                    z = U[jj][i];
                    U[jj][j] = y * c + z * s;
                    U[jj][i] = z * c - y * s;
                  }
              }
            rv1[l] = (T)0.0;
            rv1[k] = f;
            W[k] = x;
          }
      }
  }

  template <typename T>
  GMatr<T> pinv(const GMatr<T>& A)
  {
    GMatr<T> U, V, x, tmp(A.ncols(), A.nrows());
    GVect<T> W;
    CanonicalBaseGVect<T> e(0, A.nrows());
    svd(A, U, W, V);
    for (unsigned int i = 0; i < A.nrows(); i++)
      {
        e.reset(i);
        tmp.setColumn(i, dot_prod(dot_prod(dot_prod(V, GMatr<double>(DIAG, 1.0 / W, 0.0, W.size(), W.size())), t(U)), e));
      }

    return tmp;
  }

  template <typename T>
  int lu(const GMatr<T>& A, GMatr<T>& LU, GVect<unsigned int>& index)
  {
    if (A.ncols() != A.nrows())
      throw std::logic_error("Error in LU decomposition: matrix must be squared");
    int i, p, j, k, n = A.ncols(), ex;
    T val, tmp;
    GVect<T> d(n);
    LU = A;
    index.resize(n);

    ex = 1;
    for (i = 0; i < n; i++)
      {
        index[i] = i;
        val = (T)0.0;
        for (j = 0; j < n; j++)
          val = std::max(val, (T)fabs(LU[i][j]));
        if (val == (T)0.0)
          std::logic_error("Error in LU decomposition: matrix was singular");
        d[i] = val;
      }

    for (k = 0; k < n - 1; k++)
      {
        p = k;
        val = fabs(LU[k][k]) / d[k];
        for (i = k + 1; i < n; i++)
          {
            tmp = fabs(LU[i][k]) / d[i];
            if (tmp > val)
              {
                val = tmp;
                p = i;
              }
          }
        if (val == (T)0.0)
          std::logic_error("Error in LU decomposition: matrix was singular");
        if (p > k)
          {
            ex = -ex;
            std::swap(index[k], index[p]);
            std::swap(d[k], d[p]);
            for (j = 0; j < n; j++)
              std::swap(LU[k][j], LU[p][j]);
          }

        for (i = k + 1; i < n; i++)
          {
            LU[i][k] /= LU[k][k];
            for (j = k + 1; j < n; j++)
              LU[i][j] -= LU[i][k] * LU[k][j];
          }
      }
    if (LU[n - 1][n - 1] == (T)0.0)
      std::logic_error("Error in LU decomposition: matrix was singular");

    return ex;
  }

  template <typename T>
  GVect<T> lu_solve_gold(const GMatr<T>& LU, const GVect<T>& b, GVect<unsigned int>& index)
  {
    if (LU.ncols() != LU.nrows())
      throw std::logic_error("Error in LU solve: LU matrix should be squared");
    unsigned int n = LU.ncols();
    if (b.size() != n)
      throw std::logic_error("Error in LU solve: b vector must be of the same dimensions of LU matrix");
    GVect<T> x((T)0.0, n);
    int i, j, p;
    T sum;

    p = index[0];
    x[0] = b[p];

    for (i = 1; i < n; i++)
      {
        sum = (T)0.0;
        for (j = 0; j < i; j++)
          sum += LU[i][j] * x[j];
        p = index[i];
        x[i] = b[p] - sum;
      }
    x[n - 1] /= LU[n - 1][n - 1];
    for (i = n - 2; i >= 0; i--)
      {
        sum = (T)0.0;
        for (j = i + 1; j < n; j++)
          sum += LU[i][j] * x[j];
        x[i] = (x[i] - sum) / LU[i][i];
      }
    return x;
  }

  template <typename T>
  void lu_solve_gold(const GMatr<T>& LU, GVect<T>& x, const GVect<T>& b, GVect<unsigned int>& index)
  {
    x = lu_solve_gold(LU, b, index);
  }

  template <typename T>
  GMatr<T> lu_inverse(const GMatr<T>& A)
  {
    if (A.ncols() != A.nrows())
      throw std::logic_error("Error in LU invert: matrix must be squared");
    unsigned int n = A.ncols();
    GMatr<T> A1(n, n), LU;
    GVect<unsigned int> index;

    lu(A, LU, index);
    CanonicalBaseGVect<T> e(0, n);
    for (unsigned i = 0; i < n; i++)
      {
        e.reset(i);
        A1.setColumn(i, lu_solve_gold(LU, e, index));
      }

    return A1;
  }

  template <typename T>
  T lu_det(const GMatr<T>& A)
  {
    if (A.ncols() != A.nrows())
      throw std::logic_error("Error in LU determinant: matrix must be squared");
    unsigned int d;
    GMatr<T> LU;
    GVect<unsigned int> index;

    d = lu(A, LU, index);

    return d * prod(LU.extractDiag());
  }

  template <typename T>
  void cholesky(const GMatr<T> A, GMatr<T>& LL)
  {
    if (A.ncols() != A.nrows())
      throw std::logic_error("Error in Cholesky decomposition: matrix must be squared");
    int i, j, k, n = A.ncols();
    double sum;
    LL = A;

    for (i = 0; i < n; i++)
      {
        for (j = i; j < n; j++)
          {
            sum = LL[i][j];
            for (k = i - 1; k >= 0; k--)
              sum -= LL[i][k] * LL[j][k];
            if (i == j)
              {
                if (sum <= 0.0)
                  throw std::logic_error("Error in Cholesky decomposition: matrix is not postive definite");
                LL[i][i] = std::sqrt(sum);
              }
            else
              LL[j][i] = sum / LL[i][i];
          }
        for (k = i + 1; k < n; k++)
          LL[i][k] = LL[k][i];
      }
  }

  template <typename T>
  GMatr<T> cholesky(const GMatr<T> A)
  {
    GMatr<T> LL;
    cholesky(A, LL);

    return LL;
  }

  template <typename T>
  GVect<T> cholesky_solve(const GMatr<T>& LL, const GVect<T>& b)
  {
    if (LL.ncols() != LL.nrows())
      throw std::logic_error("Error in Cholesky solve: matrix must be squared");
    unsigned int n = LL.ncols();
    if (b.size() != n)
      throw std::logic_error("Error in Cholesky decomposition: b vector must be of the same dimensions of LU matrix");
    GVect<T> x, y;

    /* Solve L * y = b */
    forward_elimination(LL, y, b);
    /* Solve L^T * x = y */
    backward_elimination(LL, x, y);

    return x;
  }

  template <typename T>
  void cholesky_solve(const GMatr<T>& LL, GVect<T>& x, const GVect<T>& b)
  {
    x = cholesky_solve(LL, b);
  }

  template <typename T>
  void forward_elimination(const GMatr<T>& L, GVect<T>& y, const GVect<T> b)
  {
    if (L.ncols() != L.nrows())
      throw std::logic_error("Error in Forward elimination: matrix must be squared (lower triangular)");
    if (b.size() != L.nrows())
      throw std::logic_error("Error in Forward elimination: b vector must be of the same dimensions of L matrix");
    int i, j, n = b.size();
    y.resize(n);

    y[0] = b[0] / L[0][0];
    for (i = 1; i < n; i++)
      {
        y[i] = b[i];
        for (j = 0; j < i; j++)
          y[i] -= L[i][j] * y[j];
        y[i] = y[i] / L[i][i];
      }
  }

  template <typename T>
  GVect<T> forward_elimination(const GMatr<T>& L, const GVect<T> b)
  {
    GVect<T> y;
    forward_elimination(L, y, b);

    return y;
  }

  template <typename T>
  void backward_elimination(const GMatr<T>& U, GVect<T>& x, const GVect<T>& y)
  {
    if (U.ncols() != U.nrows())
      throw std::logic_error("Error in Backward elimination: matrix must be squared (upper triangular)");
    if (y.size() != U.nrows())
      throw std::logic_error("Error in Backward elimination: b vector must be of the same dimensions of U matrix");
    int i, j, n = y.size();
    x.resize(n);

    x[n - 1] = y[n - 1] / U[n - 1][n - 1];
    for (i = n - 2; i >= 0; i--)
      {
        x[i] = y[i];
        for (j = i + 1; j < n; j++)
          x[i] -= U[i][j] * x[j];
        x[i] = x[i] / U[i][i];
      }
  }

  template <typename T>
  GVect<T> backward_elimination(const GMatr<T>& U, const GVect<T> y)
  {
    GVect<T> x;
    forward_elimination(U, x, y);

    return x;
  }

  /* Setting default linear systems machinery */

#define det lu_det
//#define inverse lu_inverse
#define solve_gold lu_solve_gold

  /* Random */

  template <typename T>
  void random(GMatr<T>& m)
  {
    for (unsigned int i = 0; i < m.nrows(); i++)
      for (unsigned int j = 0; j < m.ncols(); j++)
        m[i][j] = (T)(rand() / double(RAND_MAX));
  }

  /**
     Aggregate functions
  */

  template <typename T>
  GVect<T> sum(const GMatr<T>& m)
  {
    GVect<T> tmp((T)0, m.ncols());
    for (unsigned int j = 0; j < m.ncols(); j++)
      for (unsigned int i = 0; i < m.nrows(); i++)
        tmp[j] += m[i][j];
    return tmp;
  }

  template <typename T>
  GVect<T> r_sum(const GMatr<T>& m)
  {
    GVect<T> tmp((T)0, m.nrows());
    for (unsigned int i = 0; i < m.nrows(); i++)
      for (unsigned int j = 0; j < m.ncols(); j++)
        tmp[i] += m[i][j];
    return tmp;
  }

  template <typename T>
  T all_sum(const GMatr<T>& m)
  {
    T tmp = (T)0;
    for (unsigned int i = 0; i < m.nrows(); i++)
      for (unsigned int j = 0; j < m.ncols(); j++)
        tmp += m[i][j];
    return tmp;
  }

  template <typename T>
  GVect<T> prod(const GMatr<T>& m)
  {
    GVect<T> tmp((T)1, m.ncols());
    for (unsigned int j = 0; j < m.ncols(); j++)
      for (unsigned int i = 0; i < m.nrows(); i++)
        tmp[j] *= m[i][j];
    return tmp;
  }

  template <typename T>
  GVect<T> r_prod(const GMatr<T>& m)
  {
    GVect<T> tmp((T)1, m.nrows());
    for (unsigned int i = 0; i < m.nrows(); i++)
      for (unsigned int j = 0; j < m.ncols(); j++)
        tmp[i] *= m[i][j];
    return tmp;
  }

  template <typename T>
  T all_prod(const GMatr<T>& m)
  {
    T tmp = (T)1;
    for (unsigned int i = 0; i < m.nrows(); i++)
      for (unsigned int j = 0; j < m.ncols(); j++)
        tmp *= m[i][j];
    return tmp;
  }

  template <typename T>
  GVect<T> mean(const GMatr<T>& m)
  {
    GVect<T> res((T)0, m.ncols());
    for (unsigned int j = 0; j < m.ncols(); j++)
      {
        for (unsigned int i = 0; i < m.nrows(); i++)
          res[j] += m[i][j];
        res[j] /= m.nrows();
      }

    return res;
  }

  template <typename T>
  GVect<T> r_mean(const GMatr<T>& m)
  {
    GVect<T> res((T)0, m.rows());
    for (unsigned int i = 0; i < m.nrows(); i++)
      {
        for (unsigned int j = 0; j < m.ncols(); j++)
          res[i] += m[i][j];
        res[i] /= m.nrows();
      }

    return res;
  }

  template <typename T>
  T all_mean(const GMatr<T>& m)
  {
    T tmp = (T)0;
    for (unsigned int i = 0; i < m.nrows(); i++)
      for (unsigned int j = 0; j < m.ncols(); j++)
        tmp += m[i][j];
    return tmp / (m.nrows() * m.ncols());
  }

  template <typename T>
  GVect<T> var(const GMatr<T>& m, bool sample_correction = false)
  {
    GVect<T> res((T)0, m.ncols());
    unsigned int n = m.nrows();
    double sum, ssum;
    for (unsigned int j = 0; j < m.ncols(); j++)
      {
        sum = (T)0.0; ssum = (T)0.0;
        for (unsigned int i = 0; i < m.nrows(); i++)
          {
            sum += m[i][j];
            ssum += (m[i][j] * m[i][j]);
          }
        if (!sample_correction)
          res[j] = (ssum / n) - (sum / n) * (sum / n);
        else
          res[j] = n * ((ssum / n) - (sum / n) * (sum / n)) / (n - 1);
      }

    return res;
  }

  template <typename T>
  GVect<T> stdev(const GMatr<T>& m, bool sample_correction = false)
  {
    return sqrt(var(m, sample_correction));
  }

  template <typename T>
  GVect<T> r_var(const GMatr<T>& m, bool sample_correction = false)
  {
    GVect<T> res((T)0, m.nrows());
    double sum, ssum;
    unsigned int n = m.ncols();
    for (unsigned int i = 0; i < m.nrows(); i++)
      {
        sum = 0.0; ssum = 0.0;
        for (unsigned int j = 0; j < m.ncols(); j++)
          {
            sum += m[i][j];
            ssum += (m[i][j] * m[i][j]);
          }
        if (!sample_correction)
          res[i] = (ssum / n) - (sum / n) * (sum / n);
        else
          res[i] = n * ((ssum / n) - (sum / n) * (sum / n)) / (n - 1);
      }

    return res;
  }

  template <typename T>
  GVect<T> r_stdev(const GMatr<T>& m, bool sample_correction = false)
  {
    return sqrt(r_var(m, sample_correction));
  }

  template <typename T>
  GVect<T> max(const GMatr<T>& m)
  {
    GVect<T> res(m.ncols());
    double value;
    for (unsigned int j = 0; j < m.ncols(); j++)
      {
        value = m[0][j];
        for (unsigned int i = 1; i < m.nrows(); i++)
          value = std::max(m[i][j], value);
        res[j] = value;
      }

    return res;
  }

  template <typename T>
  GVect<T> r_max(const GMatr<T>& m)
  {
    GVect<T> res(m.nrows());
    double value;
    for (unsigned int i = 0; i < m.nrows(); i++)
      {
        value = m[i][0];
        for (unsigned int j = 1; j < m.ncols(); j++)
          value = std::max(m[i][j], value);
        res[i] = value;
      }

    return res;
  }

  template <typename T>
  GVect<T> min(const GMatr<T>& m)
  {
    GVect<T> res(m.ncols());
    double value;
    for (unsigned int j = 0; j < m.ncols(); j++)
      {
        value = m[0][j];
        for (unsigned int i = 1; i < m.nrows(); i++)
          value = std::min(m[i][j], value);
        res[j] = value;
      }

    return res;
  }

  template <typename T>
  GVect<T> r_min(const GMatr<T>& m)
  {
    GVect<T> res(m.nrows());
    double value;
    for (unsigned int i = 0; i < m.nrows(); i++)
      {
        value = m[i][0];
        for (unsigned int j = 1; j < m.ncols(); j++)
          value = std::min(m[i][j], value);
        res[i] = value;
      }

    return res;
  }



  /**
     Single element mathematical functions
  */

  template <typename T>
  GMatr<T> exp(const GMatr<T>&m)
  {
    GMatr<T> tmp(m.nrows(), m.ncols());

    for (unsigned int i = 0; i < m.nrows(); i++)
      for (unsigned int j = 0; j < m.ncols(); j++)
        tmp[i][j] = exp(m[i][j]);

    return tmp;
  }

  template <typename T>
  GMatr<T> sqrt(const GMatr<T>&m)
  {
    GMatr<T> tmp(m.nrows(), m.ncols());

    for (unsigned int i = 0; i < m.nrows(); i++)
      for (unsigned int j = 0; j < m.ncols(); j++)
        tmp[i][j] = sqrt(m[i][j]);

    return tmp;
  }

  /**
     GMatr operators
  */

  template <typename T>
  GMatr<T> kron(const GVect<T>& b, const GVect<T>& a)
  {
    GMatr<T> tmp(b.size(), a.size());
    for (unsigned int i = 0; i < b.size(); i++)
      for (unsigned int j = 0; j < a.size(); j++)
        tmp[i][j] = a[j] * b[i];

    return tmp;
  }

  template <typename T>
  GMatr<T> t(const GMatr<T>& a)
  {
    GMatr<T> tmp(a.ncols(), a.nrows());
    for (unsigned int i = 0; i < a.nrows(); i++)
      for (unsigned int j = 0; j < a.ncols(); j++)
        tmp[j][i] = a[i][j];

    return tmp;
  }

  template <typename T>
  GMatr<T> dot_prod(const GMatr<T>& a, const GMatr<T>& b)
  {
    if (a.ncols() != b.nrows())
      throw std::logic_error("Error matrix dot product: dimensions of the matrices are not compatible");
    GMatr<T> tmp(a.nrows(), b.ncols());
    for (unsigned int i = 0; i < tmp.nrows(); i++)
      for (unsigned int j = 0; j < tmp.ncols(); j++)
        {
          tmp[i][j] = (T)0;
          for (unsigned int k = 0; k < a.ncols(); k++)
            tmp[i][j] += a[i][k] * b[k][j];
        }

    return tmp;
  }

  template <typename T>
  GMatr<T> dot_prod(const GMatr<T>& a, const GVect<T>& b)
  {
    if (a.ncols() != b.size())
      throw std::logic_error("Error matrix dot product: dimensions of the matrix and the vector are not compatible");
    GMatr<T> tmp(a.nrows(), 1);
    for (unsigned int i = 0; i < tmp.nrows(); i++)
      {
        tmp[i][0] = (T)0;
        for (unsigned int k = 0; k < a.ncols(); k++)
          tmp[i][0] += a[i][k] * b[k];
      }

    return tmp;
  }

  template <typename T>
  GMatr<T> dot_prod(const GVect<T>& a, const GMatr<T>& b)
  {
    if (a.size() != b.ncols())
      throw std::logic_error("Error matrix dot product: dimensions of the vector and matrix are not compatible");
    GMatr<T> tmp(1, b.ncols());
    for (unsigned int j = 0; j < tmp.ncols(); j++)
      {
        tmp[0][j] = (T)0;
        for (unsigned int k = 0; k < a.size(); k++)
          tmp[0][j] += a[k] * b[k][j];
      }

    return tmp;
  }

  template <typename T>
  inline GMatr<double> rank(const GMatr<T> m)
  {
    GMatr<double> tmp(m.nrows(), m.ncols());
    for (unsigned int j = 0; j < m.ncols(); j++)
      tmp.setColumn(j, rank<T>(m.extractColumn(j)));

    return tmp;
  }

  template <typename T>
  inline GMatr<double> r_rank(const GMatr<T> m)
  {
    GMatr<double> tmp(m.nrows(), m.ncols());
    for (unsigned int i = 0; i < m.nrows(); i++)
      tmp.setRow(i, rank<T>(m.extractRow(i)));

    return tmp;
  }
}

#endif // define _ARRAY_HH_
