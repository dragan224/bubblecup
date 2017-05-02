#include <iostream>
#include <cstdio>
#include <vector>
using namespace std;

const int mod = 1e9 + 7;
const int maxN = 502;

int n;
long long dp[maxN];

int main() {
  cin >> n;

  dp[0] = 1;
  dp[1] = 1;
  for (long long i = 2; i <= n; i++) {
    dp[i] = i*dp[i-1] + (i-1)*dp[i-2];
    dp[i] %= mod;
    //cout << dp[i] << " ";
  }
  cout << dp[n-1] << endl;
  return 0;
}
