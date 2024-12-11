#pragma once
#include <cassert>
#include <array>
#include <chrono>

#include "Rational.h"
#include "Residue.h"

template<size_t N>
bool operator!=(const Residue<N>& a, const Residue<N>& b) {
  return !(a == b);
}

template<size_t M, size_t N, typename Field = Rational>
struct Matrix {
  using Row = std::array<Field, N>;
  using Column = std::array<Field, M>;
  std::array<Row, M> table;
  Matrix() {
    SetOperation([](Field& x) { x = 0; });
  }
  template<typename Func, typename... Args>
  void SetOperationWithIndex(Func function, Args... args) {
    for (size_t i = 0; i < M; ++i) {
      for (size_t j = 0; j < N; ++j) {
        function(table[i][j], i, j, args...);
      }
    }
  }
  template<typename Func, typename... Args>
  void SetOperation(Func function, Args... args) {
    for (size_t i = 0; i < M; ++i) {
      for (size_t j = 0; j < N; ++j) {
        function(table[i][j], args...);
      }
    }
  }

  static Matrix unityMatrix() {
    static_assert(M == N);
    Matrix m;
    for (size_t i = 0; i < M; ++i) {
      m[i, i] = 1;
    }
    return m;
  }
  Matrix<N, M> transposed() const {
    Matrix<N, M, Field> m;
    m.SetOperationWithIndex([*this](Field& x, size_t i, size_t j) {
      x = (*this)[j, i];
    });
    return m;
  }
  Row getRow(size_t i) const {
    return table[i];
  }
  Column getColumn(size_t i) const {
    Column c;
    for (size_t k = 0; k < M; ++k) {
      c[k] = table[k][i];
    }
    return c;
  }
  Field trace() const {
    static_assert(M == N);
    Field sum = 0;
    for (size_t i = 0; i < M; ++i) {
      sum += table[i, i];
    }
    return sum;
  }
  Field det() const {
    return 0;
  }
  Field rank() const {
    return 0;
  }
  Field inverted() const {
    return 0;
  }
  void invert() const {
  }
  Field& operator[](size_t i, size_t j) {
    return table[i][j];
  }
  const Field& operator[](size_t i, size_t j) const {
    return table[i][j];
  }
  Matrix& operator+=(const Matrix& a) {
    SetOperationWithIndex([a](Field& x, size_t i, size_t j) { x += a[i, j]; });
    return *this;
  }
  Matrix& operator-=(const Matrix& a) {
    SetOperationWithIndex([a](Field& x, size_t i, size_t j) { x -= a[i, j]; });
    return *this;
  }
  Matrix& operator-() {
    SetOperation([](Field& x) { x = -x; });
    return *this;
  }
  Matrix& operator*=(Field a) {
    SetOperation([a](Field& x) { x *= a; });
    return *this;
  }
  void Print() const {
    for (size_t i = 0; i < M; ++i) {
      for (size_t j = 0; j < N; ++j) {
        cout << table[i][j] << " ";
      }
      cout << '\n';
    }
  }
};
template<size_t M, typename Field = Rational>
using SquareMatrix = Matrix<M, M, Field>;
template<size_t M, size_t N, typename T, typename Field = Rational>
Matrix<M, N, Field> operator*(const Matrix<M, N, Field>& m, T a) {
  Matrix<M, N, Field> copy = m;
  copy *= a;
  return copy;
}
template<size_t M, size_t N, typename T, typename Field = Rational>
Matrix<M, N, Field> operator*(T a, const Matrix<M, N, Field>& m) {
  Matrix<M, N, Field> copy = m;
  copy *= a;
  return copy;
}
template<size_t M, size_t N, typename Field = Rational>
Matrix<M, N, Field> operator+(const Matrix<M, N, Field> a,
                              const Matrix<M, N, Field> b) {
  Matrix<M, N, Field> copy = a;
  copy += b;
  return copy;
}
template<size_t M, size_t N, typename Field = Rational>
Matrix<M, N, Field> operator-(const Matrix<M, N, Field> a,
                              const Matrix<M, N, Field> b) {
  Matrix<M, N, Field> copy = a;
  copy -= b;
  return copy;
}
// в случае если нет неявной конверсии из Field1 в Field2, то будет CE
template<size_t M, size_t N, size_t K, typename Field = Rational>
Matrix<M, K, Field> operator*(const Matrix<M, N, Field> a,
                              const Matrix<N, K, Field> b) {
  Matrix<M, K, Field> ans;
  for (size_t i = 0; i < M; ++i) {
    for (size_t j = 0; j < K; ++j) {
      for (size_t k = 0; k < N; ++k) {
        ans[i, j] += a[i, k] * b[k, j];
      }
    }
  }
  return ans;
}
