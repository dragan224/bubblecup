#include <iostream>
#include <cstdio>
#include <algorithm>
#include <unordered_map>
#include <cassert>
using namespace std;

const int base = 5449;
const int maxN = 102;
const int maxSpellLength = 102;
const int inf = 1 << 30;

string ans;
unordered_map<int, int> pos;
string spell, heta;
int n;

int main() {
  cin.sync_with_stdio(0);
  cin >> heta;
  ans.reserve(heta.size());
  
  cin >> n; 
  for (int i = 0; i < n; i++) {
    cin >> spell;
    assert(spell.size() < maxSpellLength);
    int hash = 0;
    reverse(spell.begin(), spell.end());
    for (int j = 0; j < spell.size(); j++) {
      hash = hash * base + spell[j] - '0';
    }
    if (pos.find(hash) == pos.end()) {
      pos[hash] = i;
    }
  }

  for (int i = 0; i < heta.size(); i++) {
    ans.push_back(heta[i]);

    while (1) {
      int erase = 0;
      int hash = 0;
      int spell_index = inf;

      for (int j = ans.size() - 1; j >= 0; j--) {
        if (ans.size() - j >= maxSpellLength) break;
        
        hash = hash * base + ans[j] - '0';
        
        if (pos.find(hash) != pos.end()) {
          int curr_index = pos[hash];
          if (curr_index < spell_index) {
            spell_index = curr_index;
            erase = ans.size() - j;
          }
        }
      }

      if (erase == 0) break;
      while (erase--) {
        ans.pop_back();
      }
    }
  }
  cout << ans << "\n";
}
