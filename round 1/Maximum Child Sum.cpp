#include <iostream>
#include <algorithm>
#include <vector>
using namespace std;

const int maxNodes = 200224;
const int maxBucketSize = 400;
const int maxBuckets = maxNodes;

struct Child {
  int type, id;
  Child(int _type, int _id) : id(_id), type(_type) {}
  // type = {0 - node, 1 - bucket}
  // id - id noda/bucketa
};

vector<Child> child[maxNodes];
long long value[maxNodes];

int parentBucket[maxBuckets];

int bucketSize[maxBuckets];
long long bucketSum[maxBuckets];

int bucket[maxNodes];

long long dfs(int node) {
  long long sum = value[node];
  for (int i = 0; i < child[node].size(); i++) {
    if (child[node][i].type == 1)  {
      sum += bucketSum[child[node][i].id];
    } else {
      sum += dfs(child[node][i].id);
    }
  }
  return sum;
}

inline void propagateBucketSum(int bucket_id, int sum) {
  while (bucket_id != -1) {
    bucketSum[bucket_id] += sum;
    bucket_id = parentBucket[bucket_id];
  }
}

int q;
int type, x, y;
int id = 1;
int nbuckets = 1;

int main() {
  cin.sync_with_stdio(0);

  bucketSize[0]++;
  bucket[1] = 0;
  parentBucket[0] = -1;

  cin >> q;
  while (q--) {
    cin >> type >> x;
    if (type == 1) {
      cin >> y;
      id++;

      if (bucketSize[bucket[x]] == maxBucketSize) {
        child[x].push_back(Child(1, nbuckets));
        bucketSize[nbuckets]++;
        bucket[id] = nbuckets;
        parentBucket[nbuckets] = bucket[x];
        propagateBucketSum(nbuckets, y);
        nbuckets++;
      } else {
        child[x].push_back(Child(0, id));
        propagateBucketSum(bucket[x], y);
        bucketSize[bucket[x]]++;
        bucket[id] = bucket[x];
      }
      value[id] = y;

    } else {
      long long ans = 0;
      for (int i = 0; i < child[x].size(); i++) {
        if (child[x][i].type == 1) {
          ans = max(ans, bucketSum[child[x][i].id]);
        } else {
          ans = max(ans, dfs(child[x][i].id));
        }
      }
      cout << ans << "\n";
    }
  }
  return 0;
}
