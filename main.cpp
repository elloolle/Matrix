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
  // Matrix<2, 3> b = {{0, 1,3},{0, 3,2}};
  // Matrix<4, 4> b = {{10, 30,3,-2},{1, 3,2,-2},{3,4,0,-3},{4,5,2,4}};
  Matrix<3, 4> b = {{0, 1,3,3284029},{1, 3,2,103},{3,4,0,193}};
  // TestRank();
  cout<<Matrix<30,30>::unityMatrix().det();
  // cout<<b.rank();
}
