#include <iostream>
using namespace std;
// #include "MatrixClass.h"
#include "matrix.h"
template<size_t I = 0>
void test_impl() {
  if constexpr (I < 8) {
    constexpr std::array<size_t,8> primes = {2,3,5,7,11,13,17,239};
    constexpr size_t prime = primes[I];
    Residue<prime> a = 239;
    std::cout << Residue<prime>::pow(a, prime-1) << ' ';
    test_impl<I + 1>();
  }
}
int main(){
  // array<array<int,2>,3> a = {{3,1}};
  // cout<<a[0][0];
// test_impl<>();
const Matrix<3, 2, Residue<6>> b = {{2, 0},{1, 3},{0, 4}};
}
