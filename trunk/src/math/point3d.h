#ifndef DEFORMATION_POINT3D_H 
#define DEFORMATION_POINT3D_H 
#include "vector3d.h"
#include "matrix3d.h"
#include <cassert>

namespace meshtalent {
namespace math {

template <typename T>
class Point3d {
public:
	// constructors
	Point3d() {
		memset(p, 0, 3 * sizeof(T));
	}
	explicit Point3d(const T arr[3]) {
		p[0] = arr[0];
		p[1] = arr[1];
		p[2] = arr[2];
	}
	Point3d(const T& x, const T& y, const T& z) {
		p[0] = x;
		p[1] = y;
		p[2] = z;
	}

	// access operators
	T& operator[] (int i) {
		assert(i >= 0 && i < 3);
	   	return p[i]; 
	}
	const T& operator[] (int i) const {
		assert(i >= 0 && i < 3);
	   	return p[i]; 
	}

	// operators
	Point3d<T>& operator+= (const Vector3d<T>& V) { // add a vector
		p[0] += V[0];
		p[1] += V[1];
		p[2] += V[2];
		return *this;
	}
	Point3d<T>& operator-= (const Vector3d<T>& V) { // minus a vector
		p[0] -= V[0];
		p[1] -= V[1];
		p[2] -= V[2];
		return *this;
	}
	// this function actually do the action "Point3d = Matrix3d * Point3d"
	Point3d<T>& operator*= (const Matrix3d<T>& M) { // let a matrix time self
		Point3d<T> P(*this);
		p[0] = M[0][0] * P.p[0] + M[0][1] * P.p[1] + M[0][2] * P.p[2];
		p[1] = M[1][0] * P.p[0] + M[1][1] * P.p[1] + M[1][2] * P.p[2];
		p[2] = M[2][0] * P.p[0] + M[2][1] * P.p[1] + M[2][2] * P.p[2];
		return *this;
	}
	
	// friend operators
	friend const Point3d<T> operator+ (const Point3d& P, 
			const Vector3d<T>& V) { // point add vector
		Point3d<T> P2(P);
		P2.p[0] += V[0];
		P2.p[1] += V[1];
		P2.p[2] += V[2];
		return P2;
	}
	friend const Point3d<T> operator- (const Point3d& P, 
			const Vector3d<T>& V) { // point minus vector
		Point3d<T> P2(P);
		P2.p[0] -= V[0];
		P2.p[1] -= V[1];
		P2.p[2] -= V[2];
		return P2;
	}
	friend const Point3d<T> operator* (const Matrix3d<T>& M, 
			const Point3d& P) { // matrix time point
		Point3d<T> P2(P);
		P2.p[0] = M[0][0] * P.p[0] + M[0][1] * P.p[1] + M[0][2] * P.p[2];
		P2.p[1] = M[1][0] * P.p[0] + M[1][1] * P.p[1] + M[1][2] * P.p[2];
		P2.p[2] = M[2][0] * P.p[0] + M[2][1] * P.p[1] + M[2][2] * P.p[2];
		return P2;
	}
	friend const Vector3d<T> operator- (const Point3d<T>& lhs, 
			const Point3d<T>& rhs) { // point - point = vector 
		Vector3d<T> V;
		V[0] = lhs[0] - rhs[0];
		V[1] = lhs[1] - rhs[1];
		V[2] = lhs[2] - rhs[2];
		return V;
	}
	friend bool operator== (const Point3d& lhs, 
			const Point3d& rhs) { // equality
		return (lhs.p[0] == rhs.p[0]) && 
			   (lhs.p[1] == rhs.p[1]) &&
			   (lhs.p[2] == rhs.p[2]);
	}
	// output
	friend std::ostream& operator<< (std::ostream& os, const Point3d& P) {
		os << std::setw(10) << std::setprecision(6) << P.p[0] << ' ' 
			<< std::setw(10) << std::setprecision(6) << P.p[1] << ' ' 
			<< std::setw(10) << std::setprecision(6) << P.p[2]; 
		return os;
	}
private:
	T p[3];
};


} // end of namespace math
} // end of namespace meshtalent

#endif // DEFORMATION_POINT3D_H 
