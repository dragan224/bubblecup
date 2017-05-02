#include <iostream>
#include <algorithm>

using namespace std;

const int maxN = 100224;

int t;
int n;
int a[maxN];
int dies[maxN], killer[maxN];

int main() {
  cin >> t;
  while (t--) {
    cin >> n;
    for (int i = 0; i < n; i++) {
      cin >> a[i];
      killer[i] = -1;
      dies[i] = 0;
    }
    int moves = 0;
    for (int i = 1; i < n; i++) {
      if (a[i-1] > a[i]) {
        killer[i] = i-1;
        dies[i] = 1;
      } else {
        int idx = killer[i-1];
        int maxdies = dies[i-1];
        while (idx != -1) {
          if (a[i] < a[idx]) {
            killer[i] = idx;
            dies[i] = maxdies + 1;
            break;
          } else {
            if (idx > 0) {
              maxdies = max(maxdies, dies[idx]);
            }
            idx = killer[idx];
          }
        }
      }
      moves = max(moves, dies[i]);
    }
    // for (int i = 0; i < n; i++) {
    //   cout << dies[i] << " ";
    // }
    // cout << endl;
    cout << moves << endl;
  }
  return 0;
}
