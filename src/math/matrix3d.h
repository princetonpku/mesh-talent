#ifndef DEFORMATION_MATRIX3D_H
#define DEFORMATION_MATRIX3D_H
#include <memory>
#include <cfloat>
#include "vector3d.h"
namespace meshtalent {
namespace math {

template <typename T>
class Matrix3d {
public:
	// constructors
	Matrix3d() {

	}
	explicit Matrix3d(const T* arr) {
		memcpy(row, arr, 3 * 3 * sizeof(T));
	}
	explicit Matrix3d(const T arr[3][3]) {
		memcpy(row, arr, 3 * 3 * sizeof(T));
	}
	/*Matrix3d(const Matrix3d& M) {
		memcpy(row, M.row, 3 * 3 * sizeof(T));
	}*/
	Matrix3d(const T& a11, const T& a12, const T& a13, 
			 const T& a21, const T& a22, const T& a23, 
			 const T& a31, const T& a32, const T& a33) {
		row[0][0] = a11;
		row[0][1] = a12;
		row[0][2] = a13;
		row[1][0] = a21;
		row[1][1] = a22;
		row[1][2] = a23;
		row[2][0] = a31;
		row[2][1] = a32;
		row[2][2] = a33;
	}

	// access operators
	Vector3d<T>& operator[] (int i) { 
		assert(i >= 0 && i < 3);
		return row[i]; 
	}
	const Vector3d<T>& operator[] (int i) const { 
		assert(i >= 0 && i < 3);
		return row[i]; 
	}

	// operators
	double determinant() const {
		return row[0][0] * row[1][1] * row[2][2] + 
			   row[0][1] * row[1][2] * row[2][0] + 
			   row[0][2] * row[1][0] * row[2][1] -
			   row[0][2] * row[1][1] * row[2][0] -
			   row[0][1] * row[1][0] * row[2][2] -
			   row[0][0] * row[1][2] * row[2][1];
	}
	const Matrix3d<T> transpose() const { // transpose the matrix
		Matrix3d<T> M(*this);
		using std::swap;
		swap(M.row[0][1], M.row[1][0]);
		swap(M.row[0][2], M.row[2][0]);
		swap(M.row[1][2], M.row[2][1]);
		return M;
	}
	const Matrix3d<T> invert() const { // invert the matrix
		Matrix3d<T> invM;
		T d;

		d = row[0][0]*row[1][1]*row[2][2] + row[0][1]*row[1][2]*row[2][0] +
			row[0][2]*row[2][1]*row[1][0] - row[0][2]*row[1][1]*row[2][0] - 
			row[0][1]*row[1][0]*row[2][2] - row[0][0]*row[2][1]*row[1][2];
		assert(fabs(d) > FLT_EPSILON);

		invM[0][0] = (row[1][1]*row[2][2] - row[1][2]*row[2][1]) / d;
		invM[0][1] = (row[0][2]*row[2][1] - row[0][1]*row[2][2]) / d;
		invM[0][2] = (row[0][1]*row[1][2] - row[0][2]*row[1][1]) / d;
		invM[1][0] = (row[1][2]*row[2][0] - row[1][0]*row[2][2]) / d;
		invM[1][1] = (row[0][0]*row[2][2] - row[0][2]*row[2][0]) / d;
		invM[1][2] = (row[0][2]*row[1][0] - row[0][0]*row[1][2]) / d;
		invM[2][0] = (row[1][0]*row[2][1] - row[1][1]*row[2][0]) / d;
		invM[2][1] = (row[0][1]*row[2][0] - row[0][0]*row[2][1]) / d;
		invM[2][2] = (row[0][0]*row[1][1] - row[0][1]*row[1][0]) / d;

		return invM;	
	}
	const Matrix3d<T> operator-() const { // unary nagation
		row[0] = -row[0];
		row[1] = -row[1];
		row[2] = -row[2];
	}
	Matrix3d<T>& operator+= (const Matrix3d& M) { // add
		row[0] += M.row[0];
		row[1] += M.row[1];
		row[2] += M.row[2];
		return *this;
	}
	Matrix3d<T>& operator-= (const Matrix3d& M) { // add
		row[0] -= M.row[0];
		row[1] -= M.row[1];
		row[2] -= M.row[2];
		return *this;
	}
	Matrix3d<T>& operator*= (const T& s) { // add
		row[0] *= s;
		row[1] *= s;
		row[2] *= s;
		return *this;
	}

	friend const Matrix3d<T> operator+ (const Matrix3d& lhs, 
			const Matrix3d& rhs) { // add two matrix
		Matrix3d<T> M(lhs);
		M += rhs;
		return M;
	}
	friend const Matrix3d<T> operator- (const Matrix3d& lhs, 
			const Matrix3d& rhs) { // minus two matrix
		Matrix3d<T> M(lhs);
		M -= rhs;
		return M;
	}
	friend const Matrix3d<T> operator* (const Matrix3d& M, 
			const T& s) { // multiple a value
		Matrix3d<T> M2(M);
		M2 *= s;
		return M2;
	}
	friend const Matrix3d<T> operator* (const T& s, 
			const Matrix3d& M) { // multiple a value
		Matrix3d<T> M2(M);
		M2 *= s;
		return M2;
	}
	friend const Vector3d<T> operator* (const Matrix3d& M, 
			const Vector3d<T>& V) { // multiple a matrix and a vector
		Vector3d<T> V2;
		V2[0] = M.row[0] * V;
		V2[1] = M.row[1] * V;
		V2[2] = M.row[2] * V;
		return V2;
	}
	friend bool operator== (const Matrix3d& lhs, const Matrix3d& rhs) { // check equality
		return (lhs.row[0] == rhs.row[0]) && (lhs.row[1] == rhs.row[1]) &&
			(lhs.row[2] == rhs.row[2]);
	}

	// output
	friend std::ostream& operator<< (std::ostream& os, const Matrix3d& M) {
		os << M.row[0] << '\n'
		   << M.row[1] << '\n'
		   << M.row[2] << '\n';
		return os;
	}
	
private:
	Vector3d<T> row[3];
}; // end of class Matrix3d

} // end of namespace math
} // end of namespace meshtalent

#endif // DEFORMATION_MATRIX3D_H
