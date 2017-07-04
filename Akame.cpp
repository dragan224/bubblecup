#include <iostream>
#include <algorithm>
#include <cassert>
#include <vector>
using namespace std;

struct line {
  double angle;

};

bool compare(line A, line B) {
  return A.angle < B.angle;
}
