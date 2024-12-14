#include <iostream>
using namespace std;
#include "tests.h"
template<size_t I = 0>
void test_impl() {
  if constexpr (I < 9) {
    constexpr std::array<size_t,9> primes = {2,3,5,7,11,13,17,239,99527};
    // constexpr std::array<size_t,1> primes = {99527};
    constexpr size_t prime = primes[I];
    Residue<prime> a = 17095;
    // std::cout << Residue<prime>::pow(a, prime-2)*(a-10) << ' ';
    std::cout << Residue<prime>::pow(a, prime-1) << ' ';
    test_impl<I + 1>();
  }
}
int main(){
  // array<array<int,2>,3> a = {{3,1}};
  // cout<<a[0][0];
// test_impl<>();

  // Residue<99527> a = 17095;
  // std::cout << Residue<prime>::pow(a, prime-1) << ' ';
  // auto m = SquareMatrix<5>::unityMatrix();
  // auto m1 = SquareMatrix<5>::unityMatrix();
  // m==m1;
  // m.Print();
  Matrix<2, 3> b = {{0, 1,3},{0, 3,2}};
  b.GaussMethod();
  b.Print();
  // b.inverted().Print();
  // tester();

}
