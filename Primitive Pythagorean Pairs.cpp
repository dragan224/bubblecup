#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <cassert>
#include <ctime>
#include <set>
#include <map>
#include <unordered_map>
#include <cstring>
using namespace std;


const int maxN = 1e6 + 7;
const int maxV = 1e6 + 7;
const int mod = 1e9 + 7;

int solve(vector<int>& arr);

inline int mod_inverse(int a) {
  int b = mod - 2;
  int result = 1;
  while (b > 0) {
    if (b & 1) {
      result = 1LL*result*a % mod;
    }
    b >>= 1;
    a = 1LL*a*a % mod;
  }
  return result;
}

struct Triple {
  long long a, b, c;
  Triple(long long _a, long long _b, long long _c) : 
    a(_a),
    b(_b),
    c(_c) {}
  Triple() {}
};

inline void mod_add(int& x, int delta) {
  assert(delta < mod);
  x += delta;
  if (x < 0) x += delta;
  if (x >= mod) x -= mod;
}

inline void mod_mul(int& x, int delta) {
  x = 1LL*x*delta % mod;
}

namespace Pythagora {

  int A[][3] = {
    {1, -2, 2},
    {2, -1, 2},
    {2, -2, 3}
  };
  int B[][3] = {
    {1, 2, 2},
    {2, 1, 2},
    {2, 2, 3}
  };
  int C[][3] = {
    {-1, 2, 2},
    {-2, 1, 2},
    {-2, 2, 3}
  };
  inline Triple multiply(Triple triple, int mat[][3]) {
    Triple result;
    result.a = mat[0][0]*triple.a + mat[0][1]*triple.b + mat[0][2]*triple.c;
    result.b = mat[1][0]*triple.a + mat[1][1]*triple.b + mat[1][2]*triple.c;
    result.c = mat[2][0]*triple.a + mat[2][1]*triple.b + mat[2][2]*triple.c;
    return result;
  }

  void dfs(vector<Triple>& result, Triple triple, int upper_bound) {
    if (triple.a > upper_bound || triple.b > upper_bound) return;
    result.push_back(triple); 
    dfs(result, multiply(triple, A), upper_bound);
    dfs(result, multiply(triple, B), upper_bound);
    dfs(result, multiply(triple, C), upper_bound);
  }

  vector<Triple> generate() {
    vector<Triple> triples;
    triples.reserve(170000);  // About 168k triples up to 1 million.
    dfs(triples, Triple(3, 4, 5), maxV);
    return triples;
  }
} // Pythagora

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

int pow2[maxV];
int visited[maxV];
unordered_map<int, int> cnt;
vector<int> edges[maxV];
int parent[maxV];
int global_counter = 0;
int to_decrease[maxV];

void init(vector<int>& arr) {
  memset(visited, 0, sizeof(visited));
  memset(parent, 0, sizeof(parent));
  for (int i = 0; i < maxV; i++) {
    edges[i].clear();
  }
  cnt.clear();
  global_counter = 0;
  memset(to_decrease, 0, sizeof(to_decrease));

  pow2[0] = 1;
  for (int i = 1; i < maxV; i++) {
    pow2[i] = pow2[i-1] * 2 % mod;
    if (pow2[i] >= mod) {
      pow2[i] -= mod;
    }
  }

  int n = arr.size();
  for (int i = 0; i < n; i++) {
    cnt[arr[i]]++;
  }
  
  vector<Triple> triples = Pythagora::generate();
  for (auto it: triples) {
    int a = it.a;
    int b = it.b;
    if (cnt.find(a) == cnt.end() || cnt.find(b) == cnt.end()) {
      continue;
    }
    edges[a].push_back(b);
    edges[b].push_back(a);
  }

  for (int i = 0; i < maxV; i++) {
    if (!edges[i].empty()) {
      sort(edges[i].begin(), edges[i].end());
    }
  }
}

typedef unordered_map<int, pair<int, int> > Payload;

struct Response {
  int cnt_0;
  int cnt_1;
  Payload payload;

  Response(int _cnt_0, int _cnt_1, Payload  _payload) :
    cnt_0(_cnt_0), cnt_1(_cnt_1), 
    payload(_payload) {

    }
};

inline void payload_add(Payload& payload, int key, int val_0, int val_1) {
  auto it = payload.find(key);
  if (it == payload.end()) {
    payload.insert(make_pair(key, make_pair(val_0, val_1)));
  } else {
    mod_mul(it->second.first, val_0 + val_1);
    mod_mul(it->second.second, val_0);
  }
}
Response dfs(int node) {
  visited[node] = 1;
  int cnt_0 = 1;

  int node_cnt = pow2[cnt[node]];
  mod_add(node_cnt, -1);
  int cnt_1 = node_cnt;
  Payload load;
  vector<Response> responses;

  for (auto it: edges[node]) {
    if (visited[it]) continue;
    parent[it] = node;
    responses.push_back(dfs(it)); 
  }
  for (int i = 0; i < responses.size(); i++) {
    int resp_cnt_0 = responses[i].cnt_0;
    for (auto payload: responses[i].payload) {
      if (to_decrease[payload.first] == node) {
        resp_cnt_0 -= payload.second.first;
      } else {
        payload_add(load, payload.first, payload.second.first, payload.second.second);
      }
    }
    mod_mul(cnt_1, resp_cnt_0);
    mod_mul(cnt_0, responses[i].cnt_0 + responses[i].cnt_1);
  }

  for (auto it: load) {
    int key = it.first;
    int val_0 = 1;
    int val_1 = node_cnt;
    for (int i = 0; i < responses.size(); i++) {
      bool found = false;
      for (auto payload: responses[i].payload) {
        if (payload.first == key) {
          mod_mul(val_0, payload.second.first + payload.second.second);
          mod_mul(val_1, payload.second.first);
          found = true;
        }
      }
      if (!found) {
        mod_mul(val_0, responses[i].cnt_0 + responses[i].cnt_1);
        mod_mul(val_1, responses[i].cnt_0);
      }
    }
    load[key] = make_pair(val_0, val_1);
  }

  for (auto it: edges[node]) {
    if (visited[it] == 1 && it != parent[node]) {
      to_decrease[global_counter] = it;
      payload_add(load, global_counter, 0, cnt_1);
      ++global_counter;
    }
  }

  cout << node << " " << cnt_0 << " " << cnt_1;
  if (load.size() > 0) {
    cout << " set = ";
    for (auto it: load) {
      cout << it.first << " " << it.second.first << " " << it.second.second << " | ";
    }
  }
  cout << endl;
  visited[node] = 2;
  return Response(cnt_0, cnt_1, load);
}

int solve(vector<int>& arr) {
  int ans = 1;

  for (auto it: cnt) {
    int node = it.first;
    if (visited[node]) continue;
    Response result = dfs(node);
    mod_mul(ans, result.cnt_0 + result.cnt_1);
  }
  mod_add(ans, -1);
  return ans;
}

bool is_in[25];
int solve_brute(vector<int>& arr) {
  int n = arr.size();
  int cnt = 0;
  for (int p = 1; p < 1<<n; p++) {
    memset(is_in, 0, sizeof(is_in));
    for (int i = 0; i < n; i++) {
      if (p & 1<<i) {
        is_in[i] = true;
      }
    }
    int ok = 1;
    for (int i = 0; i < n && ok; i++) {
      for (int j = i+1; j < n && ok; j++) {
        if (is_in[i] && is_in[j]
            && binary_search(edges[arr[i]].begin(), edges[arr[i]].end(), arr[j])) {
          ok = 0;
        }
      }
    }
    if (ok) {
      cout << "{";
      string s;
      for (int i = 0; i < n; i++) {
        if (is_in[i]) s += to_string(arr[i]) + ", ";
      }
      s.pop_back();
      s.pop_back();
      cout << s << "} ";
    }
    cnt += ok;
  }
  return cnt;
}

namespace Test {
  void test_init(vector<int>& arr) {
    init(arr);
    for (int i = 0; i < maxV; i++) {
      edges[i].clear();
    }
  }
  void add(int u, int v) {
    edges[u].push_back(v);
    edges[v].push_back(u);
  }
  void test1() {
    cout << "test 1:\n";
    vector<int> arr {1, 2, 3, 4, 5};
    test_init(arr);
    add(1, 2);
    add(2, 3);
    add(3, 1);
    add(1, 4);
    add(3, 4);
    add(4, 5);
    add(5, 1);
    Response res = dfs(1);
    int sol = (res.cnt_0 + res.cnt_1 - 1);
    int corect_sol = solve_brute(arr);
    cout << "\n";
    cout << "Sol = " << sol << " | Correct Sol = " << corect_sol << "\n\n";
  }

  void run_tests() {
    test1();
  }
};

int main() {
  cin.sync_with_stdio();
  Test::run_tests();
  // int n = read_int();
  // assert(n < maxN);
  // vector<int> arr;
  // arr.reserve(n);
  // for (int i = 0; i < n; i++) {
  //   int x = read_int();
  //   assert(x>0 && x < maxV);
  //   arr.push_back(x);
  // }
  // init(arr);
  // cout << solve(arr) << "\n";
  return 0;
}
