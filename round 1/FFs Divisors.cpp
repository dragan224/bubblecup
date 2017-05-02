#include <iostream>
#include <cstdio>
#include <vector>
using namespace std;

const int maxN = 10000007;

int n;
bool is_composite[maxN];
vector<int> primes;

void make_sieve() {
  for (int i = 2; i < maxN; i++) {
    if (is_composite[i]) continue;
    primes.push_back(i);
    for (int j = i*2; j < maxN; j += i) {
      is_composite[j] = true;
    }
  }
}

int main() {
  cin.sync_with_stdio(0);
  make_sieve();
  while (true) {
    cin >> n;
    if (n == 0) break;
    long long ans = 1;
    if (n == 1) {
      ans = 1;
    } else if (!is_composite[n]) {
      ans = 2;
    } else {
      int i = 0;
      vector<int> counts;
      while (n > 1) {
        int cnt = 0;
        while (n % primes[i] == 0) {
          n /= primes[i];
          cnt++;
        }
        if (cnt) {
          counts.push_back(cnt);
        }
        i++;
      }
      long long sum = 1;
      for (int i = 0; i < counts.size(); i++) {
        sum *= (counts[i]+1);
      }
      for (int i = 0; i < counts.size(); i++) {
        ans *= (sum/(counts[i]+1))*(counts[i]*(counts[i]+1))/2 + 1;
      }
    }
    cout << ans << "\n";
  }
}
