#include <iostream>
#include <cstdio>
using namespace std;

const int maxN = 1007;
const int maxK = 207;
const int inf = 1<<30;

int n, k;

int cost[maxN][maxN];
int sum[maxN][maxN];
int dp[maxK][maxN];

void go(int p, int a, int b, int l, int r) {
  int m = (a+b)/2;
  int op = -1;
  for (int i = l; i <= min(m, r); i++) {
    if (dp[p][m] > dp[p-1][i-1] + sum[i][m]) {
      op = i;
      dp[p][m] = dp[p-1][i-1] + sum[i][m];
    }
  }

  if (a <= m-1) {
    go(p, a, m-1, l, op);
  }

  if (m+1 <= b) {
    go(p, m+1, b, op, r);
  }
}

int main() {
  cin.sync_with_stdio(0);
  cin >> n >> k;

  for (int i = 1; i <= n; i++) {
    for (int j = 1; j <= n; j++) {
      cin >> cost[i][j];
    }
  }

  for (int i = 1; i <= n; i++) {
    for (int j = 1; j <= n; j++) {
      cost[i][j] += cost[i-1][j];
    }
  }

  for (int i = 1; i <= n; i++) {
    for (int j = 1; j <= n; j++) {
      cost[i][j] += cost[i][j-1];
    }
  }

  for (int i = 1; i <= n; i++) {
    for (int j = 1; j <= n; j++) {
      sum[i][j] = (cost[j][j] + cost[i-1][i-1] - cost[i-1][j] - cost[j][i-1]) / 2;
    }
  }

  for (int i = 0; i <= k; i++) {
    for (int j = 0; j <= n; j++) {
      dp[i][j] = inf;
    }
  }

  for (int i = 1; i <= n; i++) {
    dp[1][i] = sum[1][i];
  }

  for (int i = 2; i <= k; i++) {
    go(i, 1, n, 1, n);
  }

  cout << dp[k][n] << endl;
}
