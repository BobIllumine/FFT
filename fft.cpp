#include <iostream>
#include <complex>
#include <vector>
#include <math.h>
using namespace std;

typedef complex<double> base;
const double PI = acos(-1);

void fft(vector<base> &p, bool invert) {
  int n = p.size();
  if(n == 1)
    return;
  vector<base> p_e(n / 2), p_o(n / 2);
  for(int i = 0; i * 2 < n; ++i) {
    p_e[i] = p[2 * i];
    p_o[i] = p[2 * i + 1];
  }
  fft(p_e, invert);
  fft(p_o, invert);
  
  double ang = 2 * PI / n * (invert ? -1 : 1);
  base w(1), w_n(cos(ang), sin(ang));
  for(int i = 0; i * 2 < n; ++i) {
    p[i] = p_e[i] + w * p_o[i];
    p[i + n / 2] = p_e[i] - w * p_o[i];
    if(invert) {
      p[i] /= 2;
      p[i + n / 2] /= 2;
    }
    w *= w_n;
  }
}

vector<int> multiply(vector<int> const& a, vector<int> const& b) {
  vector<base> _a(a.begin(), a.end()), _b(b.begin(), b.end());
  int n = 1;
  while(n < a.size() + b.size())
    n <<= 1;
  _a.resize(n);
  _b.resize(n);

  fft(_a, false);
  fft(_b, false);
  for(int i = 0; i < n; ++i)
    _a[i] *= _b[i];
  fft(_a, true);

  vector<int> result(n);
  for(int i = 0; i < n; ++i)
    result[i] = round(_a[i].real());
  return result;
}

int main(int argc, char argv[]) {
  vector<int> a(4), b(4);
  for(auto &e : a)
    cin >> e;
  for(auto &e : b)
    cin >> e;
  vector<int> result = multiply(a, b);
  for(auto e : result)
    cout << e << ' ';
  return 0;
}