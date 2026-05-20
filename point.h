#ifndef _POINT_H_
#define _POINT_H_

#include<math.h>

class Point
{
public:
	double x;
	double y;
	Point(const Point& a) :x(a.x), y(a.y) {}
	Point(double x = 0.0f, double y = 0.0f) : x(x), y(y) {}
	Point& operator=(const Point& a) {
		x = a.x; y = a.y;
		return *this;
	}
	//重载比较运算符
	bool operator ==(const Point& a)const {
		return x == a.x && y == a.y;
	}
	bool operator !=(const Point& a) const {
		return x != a.x || y != a.y;
	}
	//加减法
	Point operator +(const Point& a) const {
		return Point(x + a.x, y + a.y);
	}
	Point operator -(const Point& a) const {
		return Point(x - a.x, y - a.y);
	}
	//点乘
	double operator *(const Point& a) const {
		return x * a.x + y * a.y;
	}
	//标量乘、除法
	Point operator *(double a) const {
		return Point(x * a, y * a);
	}
	Point operator /(double a) const {
		double oneOverA = 1.0f / a; // 没有对除零检查
		return Point(x * oneOverA, y * oneOverA);
	}

	void translate(double ix, double iy) {
		x += ix;
		y += iy;
	}
};

inline double PointMag(const Point& a) {
	return sqrt(a.x * a.x + a.y * a.y);
}

#endif // !_POINT_H_

