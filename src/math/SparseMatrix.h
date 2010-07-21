#ifndef SPARSEMATRIX_H
#define SPARSEMATRIX_H

#include <vector>
#include <utility>
#include <cassert>

namespace meshtalent {

template <typename VTYPE, typename STYPE>
class SparseMatrix {
public:
	typedef VTYPE valuetype;
	typedef STYPE sizetype;
	typedef std::pair<valuetype, sizetype> valuesizepair;
public:
	SparseMatrix();
	bool empty() const {
		assert(colptr.size() > 0 && colptr[0] == 0);
		return values.empty() && rowind.empty() && colptr.size() == 1;
	};
	void clear() {
		values.clear();
		rowind.clear();
		colptr.clear();
		colptr.push_back(0);
	}
	std::size_t colsize() const {
		return colptr.size() - 1;
	};
	void addACol(const std::vector<valuesizepair>& acol);
public:
	const std::vector<VTYPE>& getValues() const {
		return values;
	}
	const std::vector<STYPE>& getRowind() const {
		return rowind;
	}
	const std::vector<STYPE>& getColptr() const {
		return colptr;
	}
	std::vector<VTYPE>& getValues() {
		return values;
	}
	std::vector<STYPE>& getRowind() {
		return rowind;
	}
	std::vector<STYPE>& getColptr() {
		return colptr;
	}
private:
	// used by addACol.
	class AddAnElement {
	public:
		AddAnElement(std::vector<valuetype>& _values, std::vector<sizetype>& _rowind)
			: values(_values), rowind(_rowind) {}
		void operator() (const valuesizepair& e) {
			values.push_back(e.first);
			rowind.push_back(e.second);
		}
	private:
		std::vector<valuetype>& values;
		std::vector<sizetype>& rowind;
	};
private:
	std::vector<VTYPE> values;
	std::vector<STYPE> rowind;
	std::vector<STYPE> colptr;
};

template <typename VTYPE, typename STYPE>
SparseMatrix<VTYPE, STYPE>::SparseMatrix()
{
	colptr.push_back(0);
}

template <typename VTYPE, typename STYPE>
void SparseMatrix<VTYPE, STYPE>::addACol(const std::vector<valuesizepair>& acol)
{
	for_each(acol.begin(), acol.end(), AddAnElement(values, rowind));
	colptr.push_back(values.size());
}

} // end of namespace meshtalent

#endif // SPARSEMATRIX_H
