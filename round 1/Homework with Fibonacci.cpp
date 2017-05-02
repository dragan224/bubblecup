#include <iostream>
#include <cstring>
#include <cassert>
#include <algorithm>
#include <unordered_map>
using namespace std;

const int maxN = 100224;
const int mod = 1e9 + 7;
const int maxFib = 10000009;

typedef long long ll;

ll arr[maxN];
int fib[maxFib];
pair<ll, ll> tree[maxN * 4];
ll lazy[maxN*4];

unordered_map<ll,ll> Fib;

void build_tree(int node, int a, int b) {
  if(a > b) return;
    
  if(a == b) { 
    tree[node].first = fib[arr[a]]; 
    tree[node].second = fib[arr[a] - 1];
    return;
  }
  
  build_tree(node*2, a, (a+b)/2); 
  build_tree(node*2+1, 1+(a+b)/2, b); 
  
  tree[node].first = (tree[node*2].first + tree[node*2+1].first) % mod;
  tree[node].second = (tree[node*2].second + tree[node*2+1].second) % mod;
  
}

inline void update(pair<ll, ll>& value, ll d) {
  if (d == 0) return;
  d--; 

  ll hSum = value.first;
  ll lSum = value.second;
  value.first = (fib[d]*lSum + fib[d+1]*hSum) % mod;
  if (d > 0) {
    value.second = (fib[d]*hSum + fib[d-1]*lSum) % mod;
  } else {
    value.second = (fib[d]*hSum) % mod;
  }
}

inline void lazy_update(int a, int b, int node) {
  update(tree[node], lazy[node]);
  
  if (a != b) {
    lazy[node*2]= (lazy[node*2] + lazy[node]) % mod;
    lazy[node*2+1] = (lazy[node*2+1] + lazy[node]) % mod;
  }

  lazy[node] = 0;
}

void update_tree(int node, int a, int b, int i, int j, int value) {

  if (lazy[node] != 0) { 
    lazy_update(a, b, node);
  }
  
  if (a > b || a > j || b < i) return;
     
  if(a >= i && b <= j) { 
    lazy[node] += value;
    lazy_update(a, b, node);  // Nisam siguran da li treba?! (TREBA)
    return;
  }

  update_tree(node*2, a, (a+b)/2, i, j, value); 
  update_tree(1+node*2, 1+(a+b)/2, b, i, j, value); 

  tree[node].first = (tree[node*2].first + tree[node*2+1].first) % mod;
  tree[node].second = (tree[node*2].second + tree[node*2+1].second) % mod;
}

ll query_tree(int node, int a, int b, int i, int j) {
  
  if(a > b || a > j || b < i) return 0;
  
  if(lazy[node] != 0) { 
    lazy_update(a, b, node);
  }

  if(a >= i && b <= j) {
    return tree[node].first;
  }

  ll q1 = query_tree(node*2, a, (a+b)/2, i, j); 
  ll q2 = query_tree(node*2 + 1, 1+(a+b)/2, b, i, j);

  return (q1+q2) % mod;
}

int n, t;
int m;

int main() {
  fib[0] = fib[1] = 1;
  for (int i = 2; i < maxFib; i++) {
    fib[i] = fib[i-1] + fib[i-2];
    if (fib[i] >= mod) {
      fib[i] -= mod;
    }
  }
  //assert(false);

  cin.sync_with_stdio(0);
  cin >> t;
  while (t--) {
    memset(tree, 0, sizeof(tree));
    memset(lazy, 0, sizeof(lazy));

    cin >> n;
    for (int i = 0; i < n; i++) {
      cin >> arr[i];
    }
    build_tree(1, 0, n-1);

    cin >> m;
    while (m--) {
      int type, l, r, d;
      cin >> type >> l >> r;
      --l; --r;
      if (type == 1) {
        cout << query_tree(1, 0, n-1, l, r) % mod << "\n"; 
      } else {
        cin >> d;
        assert(d >= 0);
        if (d > 0) {
          update_tree(1, 0, n-1, l, r, d);
        }
      }
    }
  }
  return 0;
}
