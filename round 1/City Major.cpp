#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <ctime>
#include <cstdlib>

using namespace std;

const double eps = 1e-9;
const double pi = 3.14159265358979323846;

// Parametri binarne pretrage
const double max_radius = 150224; // upper bound za binarnu pretragu (10^5 * sqrt(2))
const double precision = 1e-5; // preciznost binarne pretrage, while(abs(hi-lo) > precision)


struct Point {
  double x;
  double y;
  Point(double _x, double _y) : x(_x), y(_y) {}
};

struct Segment {
  Point A;
  Point B;
  Segment(Point _A, Point _B) : A(_A), B(_B) {}
};

inline int sgn(double x) {
  if (x < 0) return -1;
  return 1;
}

struct Circle {
  double radius;
  Point center;
  Circle(Point _center, double _radius) : center(_center), radius(_radius) {}
};

namespace Geometry {

vector<Point> removeDuplicatePoints(vector<Point> v) {
  vector<Point> res;
  res.push_back(v[0]);
  for (int i = 1; i < v.size(); i++) {
    if (v[i].x != v[i-1].x || v[i].y != v[i-1].y) {
      res.push_back(v[i]);
    }
  }
  while (res.size() > 1 && res.back().x == res[0].x && res.back().y == res[0].y) {
    res.pop_back();
  }
  return res;
}

bool comparexy(Point A, Point B) {
  if (A.x != B.x) {
    return A.x < B.x;
  }
  return A.y < B.y;
}

Point intersectLine(Point p0, Point p1, Point p2, Point p3) {
  // Ne podzrava pararelleno ili kada se poklapaju prave
  double s1_x, s1_y, s2_x, s2_y;
  s1_x = p1.x - p0.x;     s1_y = p1.y - p0.y;
  s2_x = p3.x - p2.x;     s2_y = p3.y - p2.y;

  double s, t;
  s = (-s1_y * (p0.x - p2.x) + s1_x * (p0.y - p2.y)) / (-s2_x * s1_y + s1_x * s2_y);
  t = ( s2_x * (p0.y - p2.y) - s2_y * (p0.x - p2.x)) / (-s2_x * s1_y + s1_x * s2_y);

  return Point(p0.x + (t * s1_x), p0.y + (t * s1_y));
}

inline double ccw(Point A, Point B, Point C) {
  double vx = A.x - B.x, vy = A.y - B.y;
  double ux = C.x - B.x, uy = C.y - B.y;
  return vx*uy - vy*ux;
}

void sortCCW(vector<Point>& Points) {
  //sorita tacke po ccw redosledu
  sort(Points.begin(), Points.end(), comparexy);
  if (Points.size() <= 2) return;

  vector<Point> upper_half;
  vector<Point> bottom_half;

  bottom_half.push_back(Points[0]);
  for (int i = 1; i < Points.size()-1; i++) {
    if (ccw(Points[0], Points[i], Points.back()) <= 0) {
      bottom_half.push_back(Points[i]);
    } else {
      upper_half.push_back(Points[i]);
    }
  }
  bottom_half.push_back(Points.back());
  reverse(upper_half.begin(), upper_half.end());
  Points.clear();
  Points.insert(Points.end(), bottom_half.begin(), bottom_half.end());
  Points.insert(Points.end(), upper_half.begin(), upper_half.end());
}

double area(const vector<Point>& Points) {
  if (Points.size() < 3) {
    return 0;
  }
  double area = 0;
  for (int i = 0; i < Points.size(); i++) {
    int j = i + 1;
    if (j == Points.size()) j = 0;
    area += Points[i].x*Points[j].y - Points[i].y*Points[j].x;
  }
  return abs(area) * 0.5;
}

Point centroid(const vector<Point>& Points) {
  double sum_x = 0;
  double sum_y = 0;
  for (int i = 0; i < Points.size(); i++) {
    sum_x += Points[i].x;
    sum_y += Points[i].y;
  }
  return Point(sum_x/Points.size(), sum_y/Points.size());
}

} // Geometry


namespace Delunay {

static const double EPSILON=0.0000000001f;

class Vector2d
{
public:
  Vector2d(double x,double y)
  {
    Set(x,y);
  };
  double GetX(void) const { return mX; };
  double GetY(void) const { return mY; };
  void  Set(double x,double y)
  {
    mX = x;
    mY = y;
  };
private:
  double mX;
  double mY;
};

typedef std::vector< Vector2d > Vector2dVector;

class Triangulate
{
public:
  static double Area(const Vector2dVector &contour) {
    int n = contour.size();
    double A=0.0f;
    for(int p=n-1,q=0; q<n; p=q++)
    {
      A+= contour[p].GetX()*contour[q].GetY() - contour[q].GetX()*contour[p].GetY();
    }
    return A*0.5f;
  }

  static bool InsideTriangle(double Ax, double Ay,
                      double Bx, double By,
                      double Cx, double Cy,
                      double Px, double Py) {
    double ax, ay, bx, by, cx, cy, apx, apy, bpx, bpy, cpx, cpy;
    double cCROSSap, bCROSScp, aCROSSbp;

    ax = Cx - Bx;  ay = Cy - By;
    bx = Ax - Cx;  by = Ay - Cy;
    cx = Bx - Ax;  cy = By - Ay;
    apx= Px - Ax;  apy= Py - Ay;
    bpx= Px - Bx;  bpy= Py - By;
    cpx= Px - Cx;  cpy= Py - Cy;

    aCROSSbp = ax*bpy - ay*bpx;
    cCROSSap = cx*apy - cy*apx;
    bCROSScp = bx*cpy - by*cpx;

    return ((aCROSSbp >= 0.0f) && (bCROSScp >= 0.0f) && (cCROSSap >= 0.0f));
  }

  static bool Snip(const Vector2dVector &contour,int u,int v,int w,int n,int *V) {
    int p;
    double Ax, Ay, Bx, By, Cx, Cy, Px, Py;

    Ax = contour[V[u]].GetX();
    Ay = contour[V[u]].GetY();

    Bx = contour[V[v]].GetX();
    By = contour[V[v]].GetY();

    Cx = contour[V[w]].GetX();
    Cy = contour[V[w]].GetY();

    if ( EPSILON > (((Bx-Ax)*(Cy-Ay)) - ((By-Ay)*(Cx-Ax))) ) return false;

    for (p=0;p<n;p++)
    {
      if( (p == u) || (p == v) || (p == w) ) continue;
      Px = contour[V[p]].GetX();
      Py = contour[V[p]].GetY();
      if (InsideTriangle(Ax,Ay,Bx,By,Cx,Cy,Px,Py)) return false;
    }
    return true;
  }

  static bool Process(const Vector2dVector &contour,
                      Vector2dVector &result) {
  /* allocate and initialize list of Vertices in polygon */

    int n = contour.size();
    if ( n < 3 ) return false;

    int *V = new int[n];

    /* we want a counter-clockwise polygon in V */

    if ( 0.0f < Area(contour) )
      for (int v=0; v<n; v++) V[v] = v;
    else
      for(int v=0; v<n; v++) V[v] = (n-1)-v;

    int nv = n;

    /*  remove nv-2 Vertices, creating 1 triangle every time */
    int count = 2*nv;   /* error detection */

    for(int m=0, v=nv-1; nv>2; )
    {
      /* if we loop, it is probably a non-simple polygon */
      if (0 >= (count--))
      {
        //** Triangulate: ERROR - probable bad polygon!
        return false;
      }

      /* three consecutive vertices in current polygon, <u,v,w> */
      int u = v  ; if (nv <= u) u = 0;     /* previous */
      v = u+1; if (nv <= v) v = 0;     /* new v    */
      int w = v+1; if (nv <= w) w = 0;     /* next     */

      if ( Snip(contour,u,v,w,nv,V) )
      {
        int a,b,c,s,t;

        /* true names of the vertices */
        a = V[u]; b = V[v]; c = V[w];

        /* output Triangle */
        result.push_back( contour[a] );
        result.push_back( contour[b] );
        result.push_back( contour[c] );

        m++;

        /* remove v from remaining polygon */
        for(s=v,t=v+1;t<nv;s++,t++) V[s] = V[t]; nv--;

        /* resest error detection counter */
        count = 2*nv;
      }
    }
    delete V;
    return true;
  }
};

} // Delunay

int n, p;
double x, y;

void printv(vector<Point> output) {
  for (int i = 0; i < output.size(); i++) {
    cout << output[i].x << " " << output[i].y << "\n";
  }
}

vector<vector<Point> > triangulate(vector<Point> Points) {
  Points = Geometry::removeDuplicatePoints(Points);
  Delunay::Vector2dVector a;
  for (int i = 0; i < Points.size(); i++) {
    a.push_back(Delunay::Vector2d(Points[i].x, Points[i].y));
  }

  Delunay::Vector2dVector result;
  Delunay::Triangulate::Process(a, result);

  vector<vector<Point> > triangles;
  int tcount = result.size() / 3;
  for (int i = 0; i < tcount; i++) {
    const Delunay::Vector2d &p1 = result[i*3+0];
    const Delunay::Vector2d &p2 = result[i*3+1];
    const Delunay::Vector2d &p3 = result[i*3+2];
    vector<Point> triangle;
    triangle.push_back(Point(p1.GetX(), p1.GetY()));
    triangle.push_back(Point(p2.GetX(), p2.GetY()));
    triangle.push_back(Point(p3.GetX(), p3.GetY()));
    triangles.push_back(triangle);
  }
  return triangles;
}

// MARKO START
inline double sqr(double x) {
    return x * x;
}

// double sgn(double x) {
//     return (x < 0.0) ? -1 : 1;
// }

Point displacement(const Point& initialPoint, const Point& constPoint) {
    return Point(constPoint.x - initialPoint.x, constPoint.y - initialPoint.y);
}

double wedgeProduct(const Point& firstFactor, const Point& secondFactor) {
    return firstFactor.x * secondFactor.y - firstFactor.y * secondFactor.x;
}

double dotProduct(const Point& firstFactor, const Point& secondFactor) {
    return firstFactor.x * secondFactor.x + firstFactor.y * secondFactor.y;
}

double crossProduct(const Point& start, const Point& next, const Point& last) {
    return wedgeProduct(displacement(start, next), displacement(start, last));
}

double onLineSegment(const Segment& segment, const Point& point) {
    return abs(crossProduct(segment.A, segment.B, point)) < eps;
}

bool onTheSameSide(const Segment& segment, const Point& point1, const Point& point2) {
    return sgn(crossProduct(segment.A, segment.B, point1)) ==
        sgn(crossProduct(segment.A, segment.B, point2));
}

double squareDistance(const Point& a, const Point& b) {
    return sqr(a.x - b.x) + sqr(a.y - b.y);
}

double distance(const Point& a, const Point& b) {
    return sqrt(squareDistance(a, b));
}

double distanceToBase(const Point& a) {
    return distance(Point(0, 0), a);
}

bool inCircle(const Circle& circle, const Point& point) {
    return squareDistance(circle.center, point) <= sqr(circle.radius);
}

double triangleArea(const Point& a, const Point& b, const Point& c) {
    double ret = 0.5 * abs(a.x * b.y + b.x * c.y + c.x * a.y - a.x * c.y - b.x * a.y - c.x * b.y);
    assert(ret >= 0.0);
    return ret;
}

// double pizzaArea(const Point& a, const Point& b, double radius) {
//     double angle = atan2(a.y, a.x) - atan2(b.y, b.x);
//     if (angle < 0.0) angle += 2 * pi;
//     angle = min(angle, 2 * pi - angle);
//     //double angle = acos(dotProduct(a, b) / (distanceToBase(a) * distanceToBase(b)));
//     assert(angle >= 0.0);
//     double ret = (angle / (2.0 * pi)) * pi * sqr(radius);
//     assert(ret >= 0.0);
//     return ret;
// }

double pizzaArea(const Point& a, const Point& b, double radius) {
    double angle = atan2(a.y, a.x) - atan2(b.y, b.x);
    if (angle < 0.0) angle += 2 * pi;
    angle = min(angle, 2 * pi - angle);
    //double angle = acos(dotProduct(a, b) / (distanceToBase(a) * distanceToBase(b)));
    assert(angle >= 0.0);
    double ret = 0.5 * angle * sqr(radius);
    assert(ret >= 0.0);
    return ret;
}

bool inSegment(Point a, Point b, Point i) {
    return min(a.x, b.x) <= i.x && max(a.x, b.x) >= i.x && 
           min(a.y, b.y) <= i.y && max(a.y, b.y) >= i.y;
}

// bool inSegment(Point a, Point b, Point i) {
//     return min(a.x, b.x) < i.x + eps && max(a.x, b.x) + eps > i.x && 
//            min(a.y, b.y) < i.y + eps && max(a.y, b.y) + eps > i.y;
// }

double calculateSegmentArea(const Segment& segment, const Circle& circle) {
    // Segment inside circle
    if (inCircle(circle, segment.A) && inCircle(circle, segment.B)) {
        return triangleArea(circle.center, segment.A, segment.B);
    }

    // Do some math
    Point segmentA = displacement(circle.center, segment.A);
    Point segmentB = displacement(circle.center, segment.B);
    Point d = displacement(segmentA, segmentB);
    double dr2 = dotProduct(d, d);
    double r2 = sqr(circle.radius);
    double D = wedgeProduct(segmentA, segmentB);
    double discriminant = r2 * dr2 - D * D;

    // No intersection
    if (discriminant < 0 || abs(discriminant) < eps) {
        return pizzaArea(segmentA, segmentB, circle.radius);
    }

    // Find intersections
    double y1 = (-D * d.x + abs(d.y) * sqrt(discriminant)) / dr2;
    double y2 = (-D * d.x - abs(d.y) * sqrt(discriminant)) / dr2;
    double x1 = (D * d.y + sgn(d.y) * d.x * sqrt(discriminant)) / dr2;
    double x2 = (D * d.y - sgn(d.y) * d.x * sqrt(discriminant)) / dr2;

    Point intersectionA = Point(x1, y1);
    Point intersectionB = Point(x2, y2);

    bool validA = false, validB = false;

    // cout << "Intersections: " << endl;
    // cout << intersectionA.x << " " << intersectionA.y << endl;
    // cout << intersectionB.x << " " << intersectionB.y << endl;

    if (inSegment(segmentA, segmentB, intersectionA)) {
        validA = true;
    }
    if (inSegment(segmentA, segmentB, intersectionB)) {
        validB = true;
    }
    // cout << validA << " " << validB << endl;
    if (!validA && !validB) {
        // No intersection with segment
        return pizzaArea(segmentA, segmentB, circle.radius);
    } else if (!validA) {
        // One point in, one out, intersectionB is valid
        if (inCircle(circle, segment.A)) {
            return triangleArea(Point(0, 0), segmentA, intersectionB)
                + pizzaArea(intersectionB, segmentB, circle.radius);
        } else {
            return triangleArea(Point(0, 0), segmentB, intersectionB)
                + pizzaArea(intersectionB, segmentA, circle.radius);
        }
    } else if (!validB) {
        // One point in, one out
        if (inCircle(circle, segment.A)) {
            return triangleArea(Point(0, 0), segmentA, intersectionA)
                + pizzaArea(intersectionA, segmentB, circle.radius);
        } else {
            return triangleArea(Point(0, 0), segmentB, intersectionA)
                + pizzaArea(intersectionA, segmentA, circle.radius);
        }
    } else {
        // Both intersections valid, swap them to match segments
        if (squareDistance(segmentA, intersectionA) > squareDistance(segmentA, intersectionB)) {
            Point temp = intersectionA;
            intersectionA = intersectionB;
            intersectionB = temp;
        }

        // Compute area
        return pizzaArea(segmentA, intersectionA, circle.radius)
            + pizzaArea(segmentB, intersectionB, circle.radius)
            + triangleArea(Point(0, 0), intersectionA, intersectionB);
    }
}

double allIn(const vector<Point>& triangle, const Circle& circle) {
    for (int i = 0; i < 3; i++) {
        if (!inCircle(circle, triangle[i])) {
            return false;
        }
    }
    return true;
}

double intersectCircleTriangle(const vector<Point>& triangle, const Circle& circle) {
  assert(triangle.size() == 3);
  double area = 0.0;

  if (allIn(triangle, circle)) {
      return triangleArea(triangle[0], triangle[1], triangle[2]);
  }

  if (circle.radius < eps || abs(triangleArea(triangle[0], triangle[1], triangle[2])) < eps) {
      return 0.0;
  }
  int cnt = 0;
  for (int i = 0; i < 3; i++) {
      Segment side = Segment(triangle[i], triangle[(i + 1) % 3]);

      if (onLineSegment(side, circle.center)) {
          continue;
      }

      double segmentArea = calculateSegmentArea(side, circle);
      // cout << segmentArea << endl;
      assert(segmentArea >= 0);

      if (onTheSameSide(side, circle.center, triangle[(i + 2) % 3])) {
          area += segmentArea;
      } else {
          cnt++;
          area -= segmentArea;
      }
  }
  assert(cnt < 3);
  if (abs(area) < eps) {
      return 0.0;
  }
  // assert(area >= 0);
  return area;
}
// MARKO END

double intersectionArea(vector<vector<Point> > triangles, Circle circle) {
  double areaSum = 0.0;
  for (int i = 0; i < triangles.size(); i++) {
    areaSum += intersectCircleTriangle(triangles[i], circle);
  }
  return areaSum;
}

namespace Tests {
bool testTriangulationAreaSum(int n, int max_c) {

  vector<Point> v;
  for (int i = 0; i < n; i++) {
    int x = rand() % max_c;
    if (rand() % 2) x = -x;
    int y = rand() % max_c;
    if (rand() % 2) y = -y;
    v.push_back(Point(x, y));
  }

  Geometry::sortCCW(v);

  double area = Geometry::area(v);
  vector<vector<Point > > tri = triangulate(v);
  double trisum = 0;
  for (int i = 0; i < tri.size(); i++) {
    trisum += Geometry::area(tri[i]);
  }

  // testiranje kruga kada je ceo poligon unutra kruga
  double circleArea = intersectionArea(tri, Circle(Point(0, 0), max_radius));
  if (abs(area - circleArea) > eps) {
    printv(v);
    exit(1);
  }
  // end

  if (abs(area - trisum) < eps) {
    return true;
  } else {
    cout << area << " " << trisum << endl;
    printv(v);
    cout << endl;
    for (int i = 0; i < tri.size(); i++) {
      printv(tri[i]);
      cout << endl;
    }
    return false;
  }
}

void run_tests(int rounds, int n, int max_c) {
  for (int i = 0; i < rounds; i++) {
    if (!testTriangulationAreaSum(n, max_c)) {
      exit(1);
    }
  } 
}

bool areaCalcTest() {
  vector<Point> v;
  v.push_back(Point(-1, 1));
  v.push_back(Point(-0.6, 0.6));
  v.push_back(Point(-0.6, 1));
  v.push_back(Point(-0.4, 0.8));
  v.push_back(Point(-0.4, 1));
  v.push_back(Point(0.2, 0.6));
  v.push_back(Point(2, 2));
  v.push_back(Point(2, 1));
  v.push_back(Point(1, 0.5));
  v.push_back(Point(1.5, 1));
  v.push_back(Point(1.75, 1.5));
  v.push_back(Point(0.5, 0.5));
  v.push_back(Point(0.75, 0.25));
  v.push_back(Point(1.5, 0.5));
  v.push_back(Point(0.5, 0));
  v.push_back(Point(0.4, 0.5));
  v.push_back(Point(0.25, -0.25));
  v.push_back(Point(1.5, -1.5));
  v.push_back(Point(0.3, -0.6));
  v.push_back(Point(0, 0));
  v.push_back(Point(0, 0.5));
  v.push_back(Point(-0.5, -1));
  v.push_back(Point(-2, -2));
  v.push_back(Point(-0.5, -0.5));
  v.push_back(Point(-1.5, 0));
  v.push_back(Point(-0.5, 0.3));

  double correctArea = 1.57;
  vector<vector<Point> > tri = triangulate(v);
  // for (int i = 0; i < tri.size(); i++) {
  //   printf("intersectPath[Polygon[(%.1lf, %.1lf), (%.1lf, %.1lf), (%.1lf, %.1lf)], circle]\n%lf\n", 
  //     tri[i][0].x, tri[i][0].y, tri[i][1].x, tri[i][1].y, tri[i][2].x, tri[i][2].y, 
  //     intersectCircleTriangle(tri[i], Circle(Point(0, 0), 1)));
  // }
  double area = intersectionArea(tri, Circle(Point(0, 0), 1));

  if (abs(correctArea - area) <= 0.01) {
    return true;
  } else {
    cout << area << " " << correctArea << endl;
    return false;
  }
}

void run() {
  srand(time(NULL));
  run_tests(500, 100, 3);
  run_tests(500, 100, 1);
  run_tests(500, 100, 5);
  run_tests(500, 100, 10);
  run_tests(500, 100, 100);
  run_tests(500, 100, 10000);
  run_tests(500, 100, 100000); 
  assert(areaCalcTest());
}
} //Tests end (trenutno ima samo test za trijangulaciju)

int main() {

  // Tests::run();

  cin.sync_with_stdio(0);
  while (cin >> n) {
    double p;
    cin >> p;
    vector<Point> Points;
    Points.reserve(n);
    for (int i = 0; i < n; i++) {
      cin >> x >> y;
      Points.push_back(Point(x, y));
    }
    Points = Geometry::removeDuplicatePoints(Points);

    double polygonArea = Geometry::area(Points);
    if (Points.size() < 3 || abs(polygonArea) < eps) {
      printf("0.00\n");
      continue;
    }
    vector<vector<Point> > triangles = triangulate(Points);

    double lo = precision, hi = max_radius;
    double r;
    double ratio = 1.0 * p / 100.0;
    while (abs(hi - lo) > precision) {

      r = (hi + lo) / 2.0;

      double circleArea = intersectionArea(triangles, Circle(Point(0, 0), r));

      // cout << r << " " << circleArea << " " << polygonArea << endl;
  
      if (circleArea + eps > ratio * polygonArea) {
        hi = r;
      } else {
        lo = r;
      }
    }
    printf("%.2lf\n", r);
  }
  return 0;
}


