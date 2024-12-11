#include <type_traits>
#include <cmath>
template<size_t N, size_t D>
struct IsPrime_Helper {
  static constexpr bool value = N % D != 0 && std::conditional_t<
    (D * D > N), IsPrime_Helper<0, 0>, IsPrime_Helper<N, D + 1> >::value;
};
template<>
struct IsPrime_Helper<0, 0> {
  static constexpr bool value = true;
};
template<size_t N>
struct IsPrime {
  static constexpr bool value = IsPrime_Helper<N, 2>::value;
};
template<>
struct IsPrime<1> {
  static constexpr bool value = false;
};
template<size_t N>
constexpr bool IsPrime_v = IsPrime<N>::value;

template<size_t N>
struct Residue {
  int number;
  Residue() = default;
  Residue(int number) {
    if (number >= 0) {
      this->number = number % N;
    } else {
      this->number = -((-number) % N);
      if (this->number != 0) {
        this->number += N;
      }
    }
  }
  operator int() const {
    return number;
  }
  static Residue pow(Residue a, size_t n) {
    if(n == 0) {
      return 1;
    }
    if(n%2 == 0) {
      return pow(a,n/2)*pow(a,n/2);
    }
    return a*pow(a,n-1);
  }
  Residue& operator+=(Residue a) {
    number = (number + a) % N;
    return *this;
  }
  Residue& operator-() {
    if (number > 0) {
      number = -number + N;
    }
    return *this;
  }
  Residue& operator-=(Residue a) {
    return operator+=(-a);
  }
  Residue& operator*=(Residue a) {
    number = static_cast<long long>(number * a) % N;
    return *this;
  }
  Residue& operator/=(Residue a) {
    static_assert(IsPrime_v<N>);
    //a^(p-1) = 1, a^(p-2) = 1/a => a^(p-2) * b
    return pow(*this,N-2)*a;
  }
};
template<size_t N>
Residue<N> operator+(Residue<N> a, Residue<N> b) {
  return a += b;
}

template<size_t N>
Residue<N> operator-(Residue<N> a, Residue<N> b) {
  return a -= b;
}

template<size_t N>
Residue<N> operator*(Residue<N> a, Residue<N> b) {
  return a *= b;
}

template<size_t N>
Residue<N> operator/(Residue<N> a, Residue<N> b) {
  return a /= b;
}

template<size_t N>
bool operator==(Residue<N> a, Residue<N> b) {
  return a.number == b.number;
}
