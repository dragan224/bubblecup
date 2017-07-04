#include <iostream>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <vector>
#include <cmath>
using namespace std;

#define a first
#define b second

const int maxN = 1000224;


long long d[maxN], p[maxN], t[maxN];
int n;
vector<pair<int, int> > S;

inline double inter(pair<int, int> x, pair<int, int> y){ return 1.* (y.b - x.b) / (x.a - y.a); }

void add(pair<int, int> x) {
  while(S.size() >= 2){
    int sz = S.size();
    if(inter(S[sz - 2], S[sz - 1]) >= inter(S[sz - 2], x)) S.pop_back();
    else break;
  }
  S.push_back(x);
}

inline long long value(int idx, long long x){ return (S[idx].a * x + S[idx].b); }

long long query(long long x){
  int lo = 0, hi = S.size() - 1, mid;
  while(lo < hi){
    mid = (lo + hi)>>1;
    if(value(mid, x) > value(mid + 1, x)) lo = mid + 1;
    else hi = mid;
  }
  return value(lo, x);
}

//dp[i] = p[i] + t[i]*d[i] + dp[j] - t[i]*d[j] 

inline int read_int() {
 register int c = getchar();
 int x;
 x = 0;
 int neg = 0;

 for(; ((c<48 || c>57) && c != '-'); c = getchar());

 if(c=='-') {
  neg = 1;
  c = getchar();
 }

 for(; c>47 && c<58 ; c = getchar()) {
  x = (x<<1) + (x<<3) + c - 48;
 }

 if(neg)
  x = -x;
 return x;
}

int main() {
  cin.sync_with_stdio(0);
  n = read_int();
  S.reserve(n);
  for (int i = 0; i < n; i++) {
    d[i] = read_int();
    p[i] = read_int();
    t[i] = read_int();
  }
  long long cost = p[n-1] + t[n-1]*d[n-1];
  add(make_pair(-d[n-1], cost));
  for (int i = n-2; i >= 0; i--) {
    cost = query(t[i]) + p[i] + t[i]*d[i];
    cost = min(cost, p[i] + t[i]*d[i]);
    add(make_pair(-d[i], cost));
  }
  cout << cost << "\n";
  return 0;
}
