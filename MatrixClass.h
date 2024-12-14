#pragma once
#include <cassert>
#include <array>
#include <chrono>
#include "Rational.h"
#include "Residue.h"
template<size_t M, size_t N, typename Field = Rational>
struct Matrix {
  using Row = std::array<Field, N>;
  using Column = std::array<Field, M>;
  std::array<Row, M> table;
  Matrix() {
    SetOperation([](Field& x) { x = 0; });
  }
  Matrix(const Matrix&) = default;
  Matrix(const std::initializer_list<std::array<Field,N>>& list) {
    auto it = list.begin();
    for (int i = 0; i < M; ++i) {
      table[i] = *it;
      ++it;
    }
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
        Matrix<M, N, Field> temp = *this;
        Field det_scale = 1;
        for (size_t i = 0; i < M; ++i) {
            // Поиск максимального элемента в столбце
            size_t pivot = i;
            for (size_t j = i + 1; j < M; ++j) {
                if (abs(temp[j, i]) > abs(temp[pivot, i])) {
                    pivot = j;
                }
            }
            if (temp[pivot, i] == 0) {
                return 0;
            }
            if (pivot != i) {
                // Меняем строки местами
                for (size_t j = 0; j < N; ++j) {
                    std::swap(temp[i, j], temp[pivot, j]);
                }
                det_scale = -det_scale;
            }
            det_scale *= temp[i, i];
            // Приведение элементов ниже диагонали к нулю
            for (size_t j = i + 1; j < M; ++j) {
                Field factor = temp[j, i] / temp[i, i];
                for (size_t k = i; k < N; ++k) {
                    temp[j, k] -= factor * temp[i, k];
                }
            }
        }
        return det_scale;
    }

  size_t rank() const {
        Matrix<M, N, Field> temp = *this;
        size_t rank = 0;
        for (size_t i = 0; i < M; ++i) {
            // Поиск ненулевого элемента
            size_t pivot = i;
            while (pivot < M && temp[pivot, i] == 0) {
                ++pivot;
            }
            if (pivot == M) {
                continue;
            }
            if (pivot != i) {
                // Меняем строки местами
                for (size_t j = 0; j < N; ++j) {
                    std::swap(temp[i, j], temp[pivot, j]);
                }
            }
            // Нормализация строки
            Field factor = temp[i, i];
            for (size_t j = 0; j < N; ++j) {
                temp[i, j] /= factor;
            }
            // Обнуление элементов в столбце
            for (size_t j = 0; j < M; ++j) {
                if (j != i && temp[j, i] != 0) {
                    Field multiple = temp[j, i];
                    for (size_t k = 0; k < N; ++k) {
                        temp[j, k] -= multiple * temp[i, k];
                    }
                }
            }
            ++rank;
        }
        return rank;
    }

  Matrix inverted() const {
        Matrix temp = *this;
        Matrix inverse;
        // Приведение матрицы к ступенчатому виду
        for (size_t i = 0; i < M; ++i) {
            // Поиск ненулевого элемента
            size_t pivot = i;
            while (pivot < M && temp[pivot, i] == 0) {
                ++pivot;
            }
            if (pivot == M) {
                throw std::runtime_error("Матрица вырожденная и не имеет обратной.");
            }
            if (pivot != i) {
                // Меняем строки местами в обеих матрицах
                for (size_t j = 0; j < N; ++j) {
                    std::swap(temp[i, j], temp[pivot, j]);
                    std::swap(inverse[i, j], inverse[pivot, j]);
                }
            }
            // Нормализация строки
            Field factor = temp[i, i];
            for (size_t j = 0; j < N; ++j) {
                temp[i, j] /= factor;
                inverse[i, j] /= factor;
            }
            // Обнуление остальных элементов в столбце
            for (size_t j = 0; j < M; ++j) {
                if (j != i && temp[j, i] != 0) {
                    Field multiple = temp[j, i];
                    for (size_t k = 0; k < N; ++k) {
                        temp[j, k] -= multiple * temp[i, k];
                        inverse[j, k] -= multiple * inverse[i, k];
                    }
                }
            }
        }
        return inverse;
    }
  void invert() {
    static_assert(M!=N);
    Matrix copy = *this;
    *this = copy.inverted();
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
  Matrix& operator*=(const Matrix& a) {
    static_assert(M==N);
    for (size_t i = 0; i < M; ++i) {
      for (size_t j = 0; j < M; ++j) {
        for (size_t k = 0; k < M; ++k) {
          (*this)[i, j] += (*this)[i, k] * a[k, j];
        }
      }
    }
    return *this;
  }
  void Print() const {
    for (size_t i = 0; i < M; ++i) {
      for (size_t j = 0; j < N; ++j) {
        std::cout << table[i][j] << " ";
      }
      std::cout << '\n';
    }
  }
};
template<size_t M, size_t N, typename T, typename Field = Rational>
bool operator==(const Matrix<M,N,Field>& a, const Matrix<M, N, Field>& b) {
  for (int i = 0; i < M; ++i) {
    for (int j = 0; j < N; ++j) {
      if(a[i,j] != b[i,j]) {
        return false;
      }
    }
  }
  return true;
}
template<size_t M, size_t N, typename T, typename Field = Rational>
bool operator!=(const Matrix<M,N,Field>& a, const Matrix<M, N, Field>& b) {
  return !(a==b);
}

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
