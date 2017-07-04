#include <cstdio>
#include <algorithm>
#include <iostream>
using namespace std;

const int mod = 1e6 + 7; // 29*34483, 2^64%mod = 11119L
const int maxN = 13;
const int maxFFT = 302;
const int maxBit = 30; 
 
const int dx[] = {-2, -2, 2, 2, -1, -1, 1, 1};
const int dy[] = {1, -1, 1, -1, 2, -2, 2, -2};

long long memo2[maxFFT];
int memo[maxFFT];

struct Matrix {
  int n;
  int** data;
  long long** tmp;
  Matrix(int _n) : n(_n) {
    data = new int*[n];
    tmp = new long long*[n];
    for (int i = 0; i < n; i++) {
      data[i] = new int[n];
      tmp[i] = new long long[n];
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
        if (data[i][j] == 0) continue;
        for (int k = 0; k < n; k++) {
          tmp[i][k] += 1LL*data[i][j]*data[j][k];
        }
      }
    }
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        if (tmp[i][j] >= mod) {
          data[i][j] = tmp[i][j] % mod;
        } else {
          data[i][j] = tmp[i][j];
        }
      }
    }
  }
};
 
inline int encode(int x, int y, int n) {
  x = min(x, 2*n - 1 - x);
  y = min(y, 2*n - 1 - y);
  if (x > y) {
    swap(x, y);
  }
  return x + y * (y+1) / 2;
}

int dp[maxN][31][maxFFT][maxFFT];


void precompute_dp(int ub) {
  for (int n = 2; n < ub; n++) {
    int size = n-1 + n * (n-1) / 2 + 1;
    Matrix fft(size + 1);

    fft.data[size][size] = 1;
    fft.data[size][0] = 1;

    for (int i = 0; i < n; i++) {
      for (int j = 0; j <= i; j++) {
        for (int k = 0; k < 8; k++) {
          int x = i + dx[k];
          int y = j + dy[k];
          if (min(x, y) >= 0 && max(x, y) < 2*n) {
            fft.data[encode(i, j, n)][encode(x, y, n)]++;
          }
        }
      }
    }

    for (int bit = 0; bit < maxBit; bit++) {
      for (int i = 0; i <= size; i++) {
        for (int j = 0; j <= size; j++) {
          dp[n][bit][i][j] = fft.data[i][j];
        }
      }
      fft.square();
    } 
  }
}

inline int solve(int n, int k) {
  int size = n-1 + n * (n-1) / 2 + 1;
  
  for (int i = 0; i <= size; i++) {
    memo[i] = memo2[i] = 0;
  }

  memo2[0] = 1;
  memo2[size] = 1;

  for (int bit = 0; bit < maxBit; bit++) {
    if (k & 1<<bit) {
      for (int i = 0; i <= size; i++) {
        memo[i] = memo2[i];
        memo2[i] = 0;
      }
      for (int i = 0; i <= size; i++) {
        for (int j = 0; j <= size; j++) {
          memo2[j] = memo2[j] + 1LL * memo[i] * dp[n][bit][i][j];
        }
      }
      for (int i = 0; i <= size; i++) {
        if (memo2[i] >= mod) {
          memo2[i] %= memo2[i];
        }
      }
    }
  }
  return memo2[0];
}


int main() {
  cin.sync_with_stdio(0);
  
  clock_t begin = clock();
  
  precompute_dp(maxN);

  int t;
  cin >> t;
  while (t--) {
    int n, k;
    cin >> n >> k;
    cout << solve(n, k) << "\n";
  }
  bool measure_time = false;
  if (measure_time) {
    clock_t end = clock();
    cout << "time = " << double(end - begin) / CLOCKS_PER_SEC << "s" << endl;
  }
  return 0;
}


