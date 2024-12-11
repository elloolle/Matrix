#pragma once
#include <algorithm>
#include <cmath>
#include <compare>
#include <iostream>
#include <string>
#include <vector>

class BigInteger;
using std::string;
using std::vector;
using std::cerr;

BigInteger operator*(const BigInteger& a, const BigInteger& b);
std::strong_ordering operator<=>(const BigInteger& a, const BigInteger& b);
bool operator==(const BigInteger& a, const BigInteger& b);
BigInteger operator/(const BigInteger& a, const BigInteger& b);
BigInteger operator+(const BigInteger& a, const BigInteger& b);
BigInteger operator-(const BigInteger& a, const BigInteger& b);

std::ostream& operator<<(std::ostream& stream, const BigInteger& n);

class BigInteger {
  std::vector<int> digits_;
  int sign_ = 0;
  static constexpr int kBase = 1e9;
  static constexpr int kCountDigits = 9;
  friend class Rational;
  void delete_last_nulls() {
    for (int i = static_cast<int>(size()) - 1; i >= 1 && digits_[i] == 0; --i) {
      digits_.pop_back();
    }
  }

  void change_to_zero() {
    digits_.clear();
    sign_ = 0;
  }

 public:
  BigInteger() = default;
  BigInteger(const BigInteger&) = default;
  explicit BigInteger(const vector<int>& n) : digits_(n), sign_(1) {}
  BigInteger& operator=(const BigInteger&) = default;

  BigInteger(long long n) {
    if (n > 0) {
      sign_ = 1;
    } else if (n == 0) {
      sign_ = 0;
      digits_.push_back(0);
      return;
    } else {
      sign_ = -1;
      n = -n;
    }

    while (n > 0) {
      digits_.push_back(n % kBase);
      n /= kBase;
    }
  }
  BigInteger(const BigInteger& n, size_t start, size_t end)
      : digits_(n.digits_.begin() + start, n.digits_.begin() + end), sign_(1) {
    delete_last_nulls();

    if (digits_.size() == 1 && digits_[0] == 0) {
      change_to_zero();
    }
  }

  BigInteger(const std::vector<int>& n, int sign) : digits_(n), sign_(sign) {}

  BigInteger(const string& s) : sign_(1) {
    int k = 0;

    if (s[0] == '-') {
      sign_ = -1;
      k = 1;
    }
    if(s == "-0") {
      sign_ = 0;
      return;
    }
    if(s[0] == '0') {
      sign_ = 0;
      return;
    }
    for (int i = static_cast<int>(s.size()) - 1; i >= k; i -= kCountDigits) {
      int start = i - kCountDigits + 1;
      int end = start + kCountDigits;

      if (start < k) {
        start = k;
        end = i + 1;
      }

      string digit(s.begin() + start, s.begin() + end);
      digits_.push_back(std::stoi(digit));
    }
  }

  BigInteger(const char* s, size_t n) : sign_(1) {
    int k = 0;

    if (s[0] == '-') {
      sign_ = -1;
      k = 1;
    }

    for (int i = static_cast<int>(n) - 1; i >= k; i -= kCountDigits) {
      int start = i - kCountDigits + 1;
      int end = i + 1;

      if (start < k) {
        start = k;
      }

      string digit(s + start, s + end);
      digits_.push_back(std::stoi(digit));
    }
  }

  size_t size() const { return digits_.size(); }

  void change_sign() { sign_ = -sign_; }

  void reserve(size_t n) { digits_.reserve(n); }

  void swap(BigInteger& n) noexcept {
    std::swap(digits_, n.digits_);
    std::swap(sign_, n.sign_);
  }

  static void swap(BigInteger& a, BigInteger& b) {
    std::swap(a.digits_, b.digits_);
    std::swap(a.sign_, b.sign_);
  }

  void reverse() { std::reverse(digits_.begin(), digits_.end()); }

 private:
  BigInteger& plus_of_positive(const std::vector<int>& n) {
    int delta = 0;

    if (size() < n.size()) {
      digits_.resize(n.size(), 0);
    }

    for (size_t i = 0; i < n.size(); ++i) {
      digits_[i] += n[i] + delta;
      delta = digits_[i] / kBase;
      digits_[i] %= kBase;
    }

    size_t i = n.size();
    while (delta != 0 && i < size()) {
      digits_[i] += delta;
      delta = digits_[i] / kBase;
      digits_[i] %= kBase;
      ++i;
    }

    if (delta != 0) {
      digits_.push_back(delta);
    }

    return *this;
  }

  BigInteger& minus_of_positive(const std::vector<int>& n) {
    int delta = 0;

    for (size_t i = 0; i < n.size() || (delta != 0 && i < size()); ++i) {
      digits_[i] -= delta;

      if (i < n.size()) {
        digits_[i] -= n[i];
      }

      if (digits_[i] < 0) {
        digits_[i] += kBase;
        delta = 1;
      } else {
        delta = 0;
      }
    }

    delete_last_nulls();

    if (size() == 1 && digits_[0] == 0) {
      sign_ = 0;
    }

    return *this;
  }

  BigInteger& composite_plus(const std::vector<int>& n, int n_sign) {
    if (sign_ == 0) {
      *this = BigInteger(n, n_sign);
      return *this;
    }

    if (n_sign == 0) {
      return *this;
    }

    if (sign_ == n_sign) {
      plus_of_positive(n);
      return *this;
    }

    if (!compare_less(digits_, n)) {
      minus_of_positive(n);
    } else {
      BigInteger copy{n, n_sign};
      swap(copy);
      minus_of_positive(copy.digits_);
    }

    return *this;
  }

  static bool compare_less(const vector<int>& a, const vector<int>& b) {
    if (a.size() != b.size()) {
      return (a.size() < b.size());
    }

    for (int i = static_cast<int>(a.size()) - 1; i >= 0; --i) {
      if (a[i] != b[i]) {
        return (a[i] < b[i]);
      }
    }

    return false;
  }

  static bool compare_great(const std::vector<int>& a, const vector<int>& b) {
    return compare_less(b, a);
  }

  static BigInteger pow_of_base(int p) {
    BigInteger b;
    b.sign_ = 1;
    b.digits_.resize(p, 0);
    b.digits_.push_back(1);
    return b;
  }

  BigInteger& stupid_multiplicate(const vector<int>& n) {
    size_t old_size = digits_.size();
    BigInteger copy;
    copy.digits_.resize(old_size + n.size(), 0);
    copy.sign_ = sign_;
    swap(copy);

    for (size_t i = 0; i < old_size; ++i) {
      int delta = 0;

      for (size_t j = 0; j < n.size(); ++j) {
        long long new_digit = digits_[i + j] +
                              static_cast<long long>(copy.digits_[i]) * n[j] +
                              delta;
        digits_[i + j] = static_cast<int>(new_digit % kBase);
        delta = (new_digit / kBase);
      }

      size_t pos = i + n.size();

      while (delta > 0) {
        if (pos >= digits_.size()) {
          digits_.push_back(delta);
          break;
        }

        int new_digit = digits_[pos] + static_cast<int>(delta);
        digits_[pos] = static_cast<int>(new_digit % kBase);
        delta = static_cast<int>(new_digit / kBase);
        ++pos;
      }
    }

    delete_last_nulls();
    return *this;
  }

  BigInteger& multiplicate_to_base(int degree) {
    digits_.resize(digits_.size() + degree);
    std::copy_backward(digits_.begin(),
                       digits_.begin() + digits_.size() - degree,
                       digits_.end());
    std::fill(digits_.begin(), digits_.begin() + degree, 0);
    return *this;
  }

  static BigInteger caracuba_multiplication(const BigInteger& A,
                                            const BigInteger& B) {
    if (A.sign_ == 0 || B.sign_ == 0) {
      return BigInteger(0);
    }

    int max_size = std::max(A.size(), B.size());

    if (max_size <= 100) {
      BigInteger copy_a = A;
      copy_a.stupid_multiplicate(B.digits_);
      return copy_a;
    }

    int size_of_half = (max_size + 1) / 2;
    BigInteger a = A;
    BigInteger b = B;
    a.sign_ = 1;
    b.sign_ = 1;
    a.digits_.resize(size_of_half * 2, 0);
    b.digits_.resize(size_of_half * 2, 0);

    BigInteger a1(a, 0, size_of_half);
    BigInteger a2(a, size_of_half, a.size());
    BigInteger b1(b, 0, size_of_half);
    BigInteger b2(b, size_of_half, b.size());
    BigInteger c1 = a1 + a2;
    BigInteger c2 = b1 + b2;
    BigInteger c = c1 * c2;
    a1 *= b1;
    a2 *= b2;
    c -= a1;
    c -= a2;
    c.multiplicate_to_base(size_of_half);
    a2.multiplicate_to_base(size_of_half * 2);
    BigInteger ans = a1 + a2 + c;
    return ans;
  }

  BigInteger multiplicate(const BigInteger& n) {
    return caracuba_multiplication(*this, n);
  }

  BigInteger multiplicate(const vector<int>& n) {
    sign_ = 1;
    return caracuba_multiplication(*this, BigInteger(n));
  }

  // эта функция находит k = a/b,
  // причём известно, что k от 0 до BASE(это достигается в том случае, если
  // a.size()-b.size() = 0 или 1)
  static int find_digit_of_division(const BigInteger& a, const vector<int>& b) {
    int l = 0;
    int r = kBase;
    BigInteger d(b);

    while (r - l > 1) {
      int m = (l + r) / 2;

      if (!compare_great((d * m).digits_, a.digits_)) {
        l = m;
      } else {
        r = m;
      }
    }

    return l;
  }

  BigInteger& divide(const vector<int>& n) {
    if (n.size() > size()) {
      change_to_zero();
      return *this;
    }

    vector<int> result;
    bool flag = false;
    BigInteger suffix(*this, size() - n.size(), size());

    for (int i = static_cast<int>(size()) - n.size() - 1; i >= -1; --i) {
      int digit = find_digit_of_division(suffix, n);

      if (digit != 0) {
        flag = true;
      }

      if (flag) {
        result.push_back(digit);
      }

      suffix -= BigInteger(digit).multiplicate(n);

      if (i == -1) {
        continue;
      }

      suffix *= kBase;
      suffix += digits_[i];
    }

    std::reverse(result.begin(), result.end());
    digits_ = std::move(result);
    delete_last_nulls();

    if (digits_.empty()) {
      sign_ = 0;
    }

    if (digits_.size() == 1 && digits_[0] == 0) {
      sign_ = 0;
    }

    return *this;
  }

 public:
  BigInteger& operator+=(const BigInteger& n) {
    return composite_plus(n.digits_, n.sign_);
  }

  BigInteger& operator-=(const BigInteger& n) {
    return composite_plus(n.digits_, -n.sign_);
  }

  BigInteger operator-() {
    BigInteger copy = *this;
    copy.change_sign();
    return copy;
  }

  BigInteger& operator++() {
    operator+=(1);
    return *this;
  }

  BigInteger& operator--() {
    operator-=(1);
    return *this;
  }

  BigInteger operator++(int) {
    BigInteger copy = *this;
    operator+=(1);
    return copy;
  }

  BigInteger operator--(int) {
    BigInteger copy = *this;
    operator-=(1);
    return copy;
  }

  BigInteger& operator*=(const BigInteger& n) {
    if (n.sign_ == 0) {
      change_to_zero();
      return *this;
    }

    if (sign_ == 0) {
      return *this;
    }

    int answer_sign = n.sign_ * sign_;
    *this = caracuba_multiplication(*this, n);
    this->sign_ = answer_sign;
    return *this;
  }

  BigInteger& operator/=(const BigInteger& n) {
    if (sign_ == 0) {
      return *this;
    }

    if (n.sign_ == 0) {
      throw std::runtime_error("division by zero");
    }

    sign_ *= n.sign_;
    return divide(n.digits_);
  }

  BigInteger& operator%=(const BigInteger& n) {
    return *this -= (*this / n) * n;
  }
  string toString() const {
    string s;

    if (sign_ == -1) {
      s.push_back('-');
    }

    if (sign_ == 0) {
      return "0";
    }

    s += std::to_string(digits_.back());

    for (int i = static_cast<int>(digits_.size()) - 2; i >= 0; --i) {
      string digit = std::to_string(digits_[i]);
      s += string(kCountDigits - digit.size(), '0') + digit;
    }

    return s;
  }

  static std::pair<vector<int>, vector<int>> divide_with_precision(
      const BigInteger& a, const BigInteger& b, int precision) {
    if (b.size() > a.size() && precision == 0) {
      return {{0}, {0}};
    }

    vector<int> integer;
    vector<int> frac;
    BigInteger suffix = 0;
    bool first_null = true;

    for (int i = static_cast<int>(a.size()) - 1; i >= -1; --i) {
      int digit = find_digit_of_division(suffix, b.digits_);

      if (digit != 0) {
        first_null = false;
      }

      if (!first_null) {
        integer.push_back(digit);
      }

      suffix -= BigInteger(digit).multiplicate(b.digits_);
      suffix *= kBase;

      if (i > -1) {
        suffix += a.digits_[i];
      }
    }

    if (integer.empty()) {
      integer.push_back(0);
    }

    for (int i = -2; i >= -1 - precision; --i) {
      int digit = find_digit_of_division(suffix, b.digits_);
      suffix -= BigInteger(digit).multiplicate(b.digits_);
      suffix *= kBase;
      frac.push_back(digit);
    }

    return {integer, frac};
  }

  static BigInteger gcd(BigInteger a, BigInteger b) {
    while (b.sign_ != 0) {
      a %= b;
      swap(a, b);
    }

    a.sign_ = 1;
    return a;
  }

  explicit operator bool() const {
    if(digits_.empty()) {
      return false;
    }
    if(digits_.size()==1 && digits_[0] == 0) {
      return false;
    }
    return true;
  }

  friend std::strong_ordering operator<=>(const BigInteger& a,
                                          const BigInteger& b) {
    if(a.sign_ == 0 && b.sign_ == 0) {
      return std::strong_ordering::equal;
    }
    if (a.sign_ != b.sign_) {
      return a.sign_ <=> b.sign_;
    }

    if (a.size() != b.size()) {
      return (a.sign_ == 1) ? (a.size() <=> b.size()) : (b.size() <=> a.size());
    }

    for (int i = static_cast<int>(a.size()) - 1; i >= 0; --i) {
      if (a.digits_[i] != b.digits_[i]) {
        return (a.sign_ == 1) ? (a.digits_[i] <=> b.digits_[i])
                              : (b.digits_[i] <=> a.digits_[i]);
      }
    }
    return std::strong_ordering::equal;
  }

  friend bool operator==(const BigInteger& a, const BigInteger& b) {
    if (a.sign_ == 0 && b.sign_ == 0) {
      return true;
    }

    if (a.size() != b.size() || a.sign_ != b.sign_) {
      return false;
    }
    return a.digits_ == b.digits_;
  }
  friend BigInteger operator"" _bi(unsigned long long n);

};

bool operator!=(const BigInteger& a, const BigInteger& b) { return !(a == b); }

BigInteger operator+(const BigInteger& a, const BigInteger& b) {
  BigInteger copy = a;
  copy += b;
  return copy;
}

BigInteger operator-(const BigInteger& a, const BigInteger& b) {
  BigInteger copy = a;
  copy -= b;
  return copy;
}

inline BigInteger operator*(const BigInteger& a, const BigInteger& b) {
  BigInteger copy = a;
  copy *= b;
  return copy;
}

inline BigInteger operator/(const BigInteger& a, const BigInteger& b) {
  BigInteger copy = a;
  copy /= b;
  return copy;
}

BigInteger operator%(const BigInteger& a, const BigInteger& b) {
  BigInteger copy = a;
  copy %= b;
  return copy;
}
BigInteger operator"" _bi(unsigned long long n) {
  BigInteger num;
  if (n > 0) {
    num.sign_ = 1;
  } else if (n == 0) {
    num.sign_ = 0;
    num.digits_.push_back(0);
    return num;
  } else {
    num.sign_ = -1;
    n = -n;
  }
  while (n > 0) {
    num.digits_.push_back(n % BigInteger::kBase);
    n /= BigInteger::kBase;
  }
  return num;
}

BigInteger operator"" _bi(const char* s, size_t n) { return {string(s, n)}; }

std::ostream& operator<<(std::ostream& stream, const BigInteger& n) {
  stream << n.toString();
  return stream;
}

std::istream& operator>>(std::istream& stream, BigInteger& n) {
  string s;
  stream >> s;
  n = BigInteger(s);
  return stream;
}
class Rational;
std::ostream& operator<<(std::ostream& stream, const Rational& n);

class Rational {
  BigInteger p;
  BigInteger q;

  void simplify() {
    BigInteger g = BigInteger::gcd(p, q);
    p /= g;
    q /= g;

    if (q.sign_ == -1) {
      q.sign_ = 1;
      p.sign_ *= -1;
    }
  }

 public:
  Rational() : p(0), q(1) {}

  Rational(const Rational&) = default;

  Rational(const BigInteger& n) : p(n), q(1) {}

  Rational(int n) : Rational(static_cast<BigInteger>(n)) {}

  Rational& operator=(const Rational&) = default;

  Rational& operator+=(const Rational& frac) {
    p = p * frac.q + q * frac.p;
    q *= frac.q;
    simplify();
    return *this;
  }

  Rational& operator-=(const Rational& frac) {
    p = p * frac.q - q * frac.p;
    q *= frac.q;
    simplify();
    return *this;
  }

  Rational& operator*=(const Rational& frac) {
    p *= frac.p;
    q *= frac.q;
    simplify();
    return *this;
  }

  Rational& operator/=(const Rational& frac) {
    p *= frac.q;
    q *= frac.p;
    p.sign_ *= frac.p.sign_;
    q.sign_ = 1;
    simplify();
    return *this;
  }

  Rational operator-() {
    Rational copy = *this;
    if (copy.p.sign_ != 0) {
      copy.p.sign_ = -copy.p.sign_;
    }
    return copy;
  }
  string toString() const {
    if (q == BigInteger(1)) {
      return p.toString();
    }

    BigInteger g = BigInteger::gcd(p, q);
    BigInteger p1 = p / g;
    BigInteger q1 = q / g;
    string s = p1.toString();
    s.push_back('/');
    s += q1.toString();
    return s;
  }

  string asDecimal(char ch, int precision) const {
    int precision1 = precision;
    precision /= BigInteger::kCountDigits;
    ++precision;

    std::pair<vector<int>, vector<int>> integer_and_frac =
        BigInteger::divide_with_precision(p, q, precision);
    string s;

    if (p.sign_ == -1) {
      s.push_back('-');
    }

    for (int el : integer_and_frac.first) {
      string digit = std::to_string(el);
      s += digit;
    }

    int count = BigInteger::kCountDigits * integer_and_frac.second.size();
    s.push_back(ch);

    for (int el : integer_and_frac.second) {
      string digit = std::to_string(el);
      s += (string(BigInteger::kCountDigits - digit.size(), '0') + digit);
    }

    while (count > precision1) {
      --count;
      s.pop_back();
    }

    if (precision == 0) {
      return s;
    }

    return s;
  }

  string asDecimal(int precision = 0) const {
    return asDecimal('.', precision);
  }

  explicit operator double() const { return std::stod(asDecimal('.', 100)); }

  friend std::strong_ordering operator<=>(const Rational& a,
                                          const Rational& b) {
    return operator<=>(a.p * b.q, a.q * b.p);
  }
    friend bool operator==(const Rational& a, const Rational& b) {
      if ((a.p != b.p) || (a.q != b.q)) {
        return false;
      }
      return true;
    }
};

bool operator!=(const Rational& a, const Rational& b) { return !(a == b); }

Rational operator+(const Rational& a, const Rational& b) {
  Rational copy = a;
  copy += b;
  return copy;
}

Rational operator-(const Rational& a, const Rational& b) {
  Rational copy = a;
  copy -= b;
  return copy;
}

Rational operator*(const Rational& a, const Rational& b) {
  Rational copy = a;
  copy *= b;
  return copy;
}

Rational operator/(const Rational& a, const Rational& b) {
  Rational copy = a;
  copy /= b;
  return copy;
}

std::ostream& operator<<(std::ostream& stream, const Rational& n) {
  std::string s = n.toString();
  stream << s;
  return stream;
}