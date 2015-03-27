#ifndef ALGEBRA_H
#define ALGEBRA_H

#include <cmath>
#include <iostream>
#include <sstream>

using namespace std;

inline int nint(double d)
{
   return floor( d + 0.5);
}

class Point2D
{
public:
  Point2D()
  {
    v_[0] = 0.0;
    v_[1] = 0.0;
  }
  Point2D(double x, double y)
  { 
    v_[0] = x;
    v_[1] = y;
  }
  Point2D(const Point2D& other)
  {
    v_[0] = other.v_[0];
    v_[1] = other.v_[1];
  }

  Point2D& operator =(const Point2D& other)
  {
    v_[0] = other.v_[0];
    v_[1] = other.v_[1];
    return *this;
  }

  const double &x() const { return v_[0]; }
  const double &y() const { return v_[1]; }

  double length2() {
    return v_[0]*v_[0] + v_[1] * v_[1];
  }
  
  double length() {
    return sqrt(length2());
  }

  int int_length() {
    return nint(length());
  }

  double& operator[](size_t idx) 
  {
    return v_[ idx ];
  }
  double operator[](size_t idx) const 
  {
    return v_[ idx ];
  }

private:
  double v_[2];
};


inline Point2D operator +(const Point2D& a, const Point2D& b)
{
  return Point2D(a[0]+b[0], a[1]+b[1]);
}

inline Point2D operator -(const Point2D& a, const Point2D& b)
{
  return Point2D(a[0]-b[0], a[1]-b[1]);
}

inline std::ostream& operator <<(std::ostream& os, const Point2D& p)
{
  return os << "p<" << p[0] << "," << p[1] << ">";
}

#endif



