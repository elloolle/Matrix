// #define CATCH_CONFIG_MAIN
//
// #include "catch.hpp"
// #include <array>
// #include <type_traits>
// #include "matrix.h"
// #include "matrix.h"  // check include guards
//
// template <size_t M, size_t N, class T>
// void EqualMatrix(const Matrix<M, N, T> &matrix, const std::array<std::array<T, N>, M> &arr) {
//   for (size_t i = 0u; i < M; ++i) {
//     for (size_t j = 0u; j < N; ++j) {
//       REQUIRE(matrix(i, j) == arr[i][j]);
//     }
//   }
// }
//
// TEST_CASE("Matrix", "[Public]") {
//   static_assert(sizeof(Matrix<17, 2, int>) == sizeof(int) * 34);
//
//   Matrix<1, 2, int> a{1, -1};
//   Matrix<2, 1, int> b{-1, 2};
//   REQUIRE(std::as_const(a).RowsNumber() == 1);
//   REQUIRE(std::as_const(a).ColumnsNumber() == 2);
//   REQUIRE(std::as_const(b).RowsNumber() == 2);
//   REQUIRE(std::as_const(b).ColumnsNumber() == 1);
//
//   a(0, 1) = -2;
//   REQUIRE(std::as_const(a)(0, 1) == -2);
//
//   b.At(1, 0) = 3;
//   REQUIRE(std::as_const(b).At(1, 0) == 3);
//
//   REQUIRE_THROWS_AS(std::as_const(a).At(5, 5), MatrixOutOfRange);
//
//   EqualMatrix(a += a, {2, -4});
//   EqualMatrix(b + b, {-2, 6});
//
//   EqualMatrix(a -= (GetTransposed(b) + GetTransposed(b)), {4, -10});
//   EqualMatrix(b - b, {0, 0});
//
//   EqualMatrix(a * b, {-34});
//   EqualMatrix(b *= Matrix<1, 1, int>{-1}, {1, -3});
//
//   EqualMatrix(a *= -1l, {-4, 10});
//   EqualMatrix(2ul * b, {2, -6});
//
//   EqualMatrix(a / 2ul, {-2, 5});
//   EqualMatrix(b /= -1l, {-1, 3});
//
//   REQUIRE(a == a);
//   REQUIRE(a != GetTransposed(b));
//
//   std::stringstream is{"-5 1\n0 10"};
//
//   Matrix<2, 2, int> matrix{};
//   is >> matrix;
//   EqualMatrix(matrix, std::array<std::array<int, 2>, 2>{-5, 1, 0, 10});
//
//   std::stringstream os;
//   os << std::as_const(matrix);
//   REQUIRE(os.str() == "-5 1\n0 10\n");
// }
//
// #ifdef MATRIX_SQUARE_MATRIX_IMPLEMENTED
//
// TEST_CASE("SquareMatrix", "[Public]") {
//   {
//     Matrix<2, 2, int> matrix{-1, 4, 9, 2};
//     Transpose(matrix);
//     EqualMatrix(matrix, std::array<std::array<int, 2>, 2>{-1, 9, 4, 2});
//   }
//
//   {
//     const Matrix<2, 2, int> matrix{-1, 4, 9, 2};
//     REQUIRE(Trace(matrix) == 1);
//   }
//
//   {
//     const Matrix<2, 2, int> matrix{-1, 4, 9, 2};
//     REQUIRE(Determinant(matrix) == -38);
//   }
//
//   {
//     Matrix<2, 2, float> matrix{-1, 4, 9, 2};
//     Inverse(matrix);
//     REQUIRE(matrix(0, 0) == Approx(-1.f / 19));
//     REQUIRE(matrix(0, 1) == Approx(2.f / 19));
//     REQUIRE(matrix(1, 0) == Approx(9.f / 38));
//     REQUIRE(matrix(1, 1) == Approx(1.f / 38));
//   }
//
//   {
//     Matrix<2, 2, float> matrix{-1, 4, 9, 2};
//     const auto inversed = GetInversed(matrix);
//     REQUIRE(inversed(0, 0) == Approx(-1.f / 19));
//     REQUIRE(inversed(0, 1) == Approx(2.f / 19));
//     REQUIRE(inversed(1, 0) == Approx(9.f / 38));
//     REQUIRE(inversed(1, 1) == Approx(1.f / 38));
//   }
//
//   {
//     Matrix<2, 2, int> matrix{};
//     REQUIRE_THROWS_AS(Inverse(matrix), MatrixIsDegenerateError);
//   }
//
//   {
//     Matrix<2, 2, int> matrix{};
//     REQUIRE_THROWS_AS(GetInversed(matrix), MatrixIsDegenerateError);
//   }
// }
//
// #endif  // MATRIX_SQUARE_MATRIX_IMPLEMENTED