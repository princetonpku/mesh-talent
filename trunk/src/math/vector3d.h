#ifndef DEFORMATION_VECTOR3D_H 
#define DEFORMATION_VECTOR3D_H 
#include <cmath>
#include <memory>
#include <cstring>
#include <ostream>
#include <iomanip>
#include <cassert>
namespace meshtalent {
namespace math {

template <typename T>
class Vector3d {
public:
	// constructors
	Vector3d() {
		memset(v, 0, 3 * sizeof(T));
	}
	explicit Vector3d(const T arr[3]) {
		v[0] = arr[0];
		v[1] = arr[1];
		v[2] = arr[2];
	}

	// to add a template constructor 
	// and a template copy constructor

	// let the compiler do the bit copy
	/*Vector3d(const Vector3d& V) {
		v[0] = V.v[0];
		v[1] = V.v[1];
		v[2] = V.v[2];
	}*/
	Vector3d(const T& x, const T& y, const T& z) {
		v[0] = x;
		v[1] = y;
		v[2] = z;
	}

	// access operators
	T& operator[] (int i) {
		assert(i >= 0 && i < 3);
	   	return v[i]; 
	}
	const T& operator[] (int i) const {
		assert(i >= 0 && i < 3);
	   	return v[i]; 
	}

	// magnitude and normalize
	T magsqr() const {
		return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
	}
	T mag() const {
		return (T)sqrt(magsqr());
	}
	const Vector3d<T> normalize() const {
		Vector3d<T> V(*this);
		T magrec = 1.0 / V.mag();
		V *= magrec;
		return V;
	}

	// operators
	/*Vector3d<T>& operator= (const Vector3d& V) { // assignment
		v[0] = V.v[0];
		v[1] = V.v[1];
		v[2] = V.v[2];
		return *this;
	}*/
	const Vector3d<T> operator- () const { // unary nagation
		Vector3d<T> V(*this);
		V.v[0] = -V.v[0];
		V.v[1] = -V.v[1];
		V.v[2] = -V.v[2];
		return V;
	}
	Vector3d<T>& operator+= (const Vector3d& V) { // add 
		v[0] += V.v[0];
		v[1] += V.v[1];
		v[2] += V.v[2];
		return *this;
	}
	Vector3d<T>& operator-= (const Vector3d& V) { // minus
		v[0] -= V.v[0];
		v[1] -= V.v[1];
		v[2] -= V.v[2];
		return *this;
	}
	Vector3d<T>& operator*= (const T& s) { // scalar multiply
		v[0] *= s;
		v[1] *= s;
		v[2] *= s;
		return *this;
	}

	// friend operators
	friend const Vector3d<T> operator+ (const Vector3d& lhs,
		   	const Vector3d& rhs) { // vector add
		Vector3d<T> V(lhs);
		V += rhs;
		return V;
	}
	friend const Vector3d<T> operator- (const Vector3d& lhs,
		   	const Vector3d& rhs) { // vector minus
		Vector3d<T> V(lhs);
		V -= rhs;
		return V;
	}
	friend const Vector3d<T> operator* (const T& s,
		   	const Vector3d& V) { // scalar multiply
		Vector3d<T> V2(V);
		V2 *= s;
		return V2;
	}
	friend const Vector3d<T> operator* (const Vector3d& V,
		   	const T& s) { // scalar multiply
		Vector3d<T> V2(V);
		V2 *= s;
		return V2;
	}
	friend T operator* (const Vector3d& lhs,
		   	const Vector3d& rhs) { // vector dot product
		return lhs.v[0] * rhs.v[0] + 
			   lhs.v[1] * rhs.v[1] + 
			   lhs.v[2] * rhs.v[2]; 
	}
	friend const Vector3d<T> operator% (const Vector3d& lhs, 
			const Vector3d& rhs) { // vector cross product
		Vector3d<T> V;
		V.v[0] = lhs.v[1] * rhs.v[2] - lhs.v[2] * rhs.v[1];
		V.v[1] = lhs.v[2] * rhs.v[0] - lhs.v[0] * rhs.v[2];
		V.v[2] = lhs.v[0] * rhs.v[1] - lhs.v[1] * rhs.v[0];
		return V;
	}
	friend bool operator== (const Vector3d& lhs, 
			const Vector3d& rhs) { // equality
		return (lhs.v[0] == rhs.v[0]) && 
			   (lhs.v[1] == rhs.v[1]) &&
			   (lhs.v[2] == rhs.v[2]);
	}
	// output
	friend std::ostream& operator<< (std::ostream& os, const Vector3d& V) {
		os << std::setw(10) << std::setprecision(6) << V.v[0] << ' ' 
		   << std::setw(10) << std::setprecision(6) << V.v[1] << ' ' 
		   << std::setw(10) << std::setprecision(6) << V.v[2]; 
		return os;
	}
	
private:
	T v[3];
}; // end of Vector3d
	
} // end of namespace math
} // end of namespace meshtalent

#endif // DEFORMATION_VECTOR3D_H 
