#include <cstdio>
#include <vector>
#include <algorithm>
#include <map>
#include <cstring>
#include <iostream>
#include <cassert>
#include <cstring>
using namespace std;

const int mod = 1e6 + 7;
const int maxN = 25;
const int maxFFT = maxN + maxN*(maxN-1)/2 + 2;
 
const int dx[] = {-2, -2, 2, 2, -1, -1, 1, 1};
const int dy[] = {1, -1, 1, -1, 2, -2, 2, -2};

int n, k;

struct Matrix {
  int n;
  int data[maxFFT][maxFFT];
  long long tmp[maxFFT][maxFFT];
  Matrix(int _n) : n(_n) {
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        data[i][j] = tmp[i][j] = 0;
      }
    }
  }
  void square() {
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        tmp[i][j] = 0;
      }
    }
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        for (int k = 0; k < n; k++) {
          tmp[i][k] = tmp[i][k] + 1LL*data[i][j] * data[j][k];
        }
      }
    }
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        data[i][j] = tmp[i][j] % mod;
      }
    }
  }
};
 
int encode(int x, int y) {
  x = min(x, 2*n - 1 - x);
  y = min(y, 2*n - 1 - y);
  if (x > y) {
    swap(x, y);
  }
  return x + y * (y+1) / 2;
}

inline int read_int() {
 register int c = getchar();
 int x;
 x = 0;

 for(; ((c<48 || c>57)); c = getchar());

 for(; c>47 && c<58 ; c = getchar()) {
  x = x * 10 + c - 48;
 }

 return x;
}

long long memo2[maxFFT];
long long memo[maxFFT];

inline int solve(int n, int k) {
  int size = n-1 + n * (n-1) / 2 + 1;
  Matrix fft(size + 1);
  
  for (int i = 0; i <= size; i++) {
    memo[i] = memo2[i] = 0;
  }

  memo2[0] = 1;
  memo2[size] = 1;
  fft.data[size][size] = 1;
  fft.data[size][0] = 1;

  for (int i = 0; i < n; i++) {
    for (int j = 0; j <= i; j++) {
      for (int k = 0; k < 8; k++) {
        int x = i + dx[k];
        int y = j + dy[k];
        if (min(x, y) >= 0 && max(x, y) < 2*n) {
          fft.data[encode(i, j)][encode(x, y)]++;
        }
      }
    }
  }

  while (k > 0) {
    if (k & 1) {
      for (int i = 0; i <= size; i++) {
        memo[i] = memo2[i];
        memo2[i] = 0;
      }
      for (int i = 0; i <= size; i++) {
        for (int j = 0; j <= size; j++) {
          memo2[j] = memo2[j] + 1LL * memo[i] * fft.data[i][j];
        }
      }
      for (int i = 0; i <= size; i++) {
        memo2[i] %= mod;
      }
      if (memo2[0] == 0) return 0;
    }
    fft.square();
    k >>= 1;
  }
  return memo2[0];
}

int main() {
  cin.sync_with_stdio(0);

  int t;
  t = read_int();
  while (t--) {
    n = read_int();
    k = read_int();
    cout << solve(n, k) << "\n";
  }

  return 0;
}
