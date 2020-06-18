#ifndef MATRIX_H
#define MATRIX_H

// Standartbibliothek
#include <vector>    // vector
#include <cstdint>   // uint...
#include <stdexcept> // exceptions
#include <cmath>     // pow

// undefine to disable rang checking
#define RANGE_CHECK

class overdetermined : public std::domain_error
{
public:
    overdetermined()
        : std::domain_error("solution is over-determined")
    {
    }
};

class underdetermined : public std::domain_error
{
public:
    underdetermined()
        : std::domain_error("solution is under-determined")
    {
    }
};

// DECLARATIONS

template <class T>
class kahn_sum
{ // implements Kahn Summation method
public:
    kahn_sum() : sum(0.0), cor(0.0) {}
    kahn_sum<T> &operator+=(const T &val)
    {
        T old_sum = sum;
        T next = val - cor;
        cor = ((sum += next) - old_sum) - next;
        return *this;
    }
    kahn_sum<T> &operator-=(const T &val)
    {
        T old_sum = sum;
        T next = val + cor;
        cor = ((sum -= val) - old_sum) + next;
        return *this;
    }
    operator T &() { return sum; }

private:
    T sum; // running sum
    T cor; // correction term
};

template <class T>
class Matrix
{
    std::vector<T> elements; // array of elements (private)

protected:
    // range check function for matrix access
    void range_check(uint32_t i, uint32_t j) const;

public:
    const uint32_t rows; // number of rows
    const uint32_t cols; // number of columns

    T &operator()(uint32_t i, uint32_t j)
    {
#ifdef RANGE_CHECK
        range_check(i, j);
#endif
        return elements[i * cols + j];
    }

    const T &operator()(uint32_t i, uint32_t j) const
    {
#ifdef RANGE_CHECK
        range_check(i, j);
#endif
        return elements[i * cols + j];
    }

    const T &element(uint32_t i, uint32_t j) const
    {
#ifdef RANGE_CHECK
        range_check(i, j);
#endif
        return elements[i * cols + j];
    }

    T &element(uint32_t i, uint32_t j)
    {
#ifdef RANGE_CHECK
        range_check(i, j);
#endif
        return elements[i * cols + j];
    }

    // constructors
    Matrix(uint32_t rows, uint32_t columns, const T *elements = 0);
    Matrix(const Matrix<T> &M);
    // destructor
    ~Matrix();

    // assignment
    Matrix<T> &operator=(const Matrix<T> &M);

    // comparison
    bool operator==(const Matrix<T> &M) const;
    bool iszero() const;
    bool operator!() const
    {
        return iszero();
    }

    // scalar multiplication/division
    Matrix<T> &operator*=(const T &a);
    Matrix<T> operator*(const T &a) const
    {
        return Matrix<T>(*this).operator*=(a);
    }
    Matrix<T> &operator/=(const T &a);
    Matrix<T> operator/(const T &a)
    {
        return Matrix<T>(*this).operator/=(a);
    }
    Matrix<T> operator-() const;
    Matrix<T> operator+() const;

    // addition/subtraction
    Matrix<T> &operator+=(const Matrix<T> &M);
    Matrix<T> &operator-=(const Matrix<T> &M);
    Matrix<T> operator+(const Matrix<T> &M) const
    {
        return Matrix<T>(*this).operator+=(M);
    }
    Matrix<T> operator-(const Matrix<T> &M) const
    {
        return Matrix<T>(*this).operator-=(M);
    }

    // matrix multiplication
    Matrix<T> operator*(const Matrix<T> &M) const;
    Matrix<T> &operator*=(const Matrix<T> &M)
    {
        return *this = *this * M;
    }

    // matrix division
    Matrix<T> leftdiv(const Matrix<T> &) const;
    Matrix<T> rightdiv(const Matrix<T> &D) const
    {
        return transpose().leftdiv(D.transpose()).transpose();
    }
    Matrix<T> operator/(const Matrix<T> &D) const
    {
        return rightdiv(D);
    }
    Matrix<T> &operator/=(const Matrix<T> &M)
    {
        return *this = *this / M;
    }

    // determinants
    Matrix<T> minor(uint32_t i, uint32_t j) const;
    T det() const;
    T minor_det(uint32_t i, uint32_t j) const;

    // these member functions are only valid for squares
    Matrix<T> inverse() const;
    Matrix<T> pow(uint32_t exp) const;
    Matrix<T> identity() const;
    bool isidentity() const;

    // vector operations
    Matrix<T> getrow(uint32_t j) const;
    Matrix<T> getcol(uint32_t i) const;
    Matrix<T> &setcol(uint32_t j, const Matrix<T> &C);
    Matrix<T> &setrow(uint32_t i, const Matrix<T> &R);
    Matrix<T> delrow(uint32_t i) const;
    Matrix<T> delcol(uint32_t j) const;

    Matrix<T> transpose() const;
    Matrix<T> operator~() const
    {
        return transpose();
    }
};

// DEFINITIONS

template <class T>
Matrix<T>::Matrix(uint32_t rows, unsigned cols, const T *elements) : rows(rows), cols(cols), elements(rows * cols, T(0.0))
{
    if (rows == 0 | cols == 0)
        throw std::range_error("attempt to create a degenerate matrix");

    // initialize from array
    if (elements)
        for (uint32_t i = 0; i < rows * cols; i++)
            this->elements[i] = elements[i];
}

template <class T>
Matrix<T>::Matrix(const Matrix<T> &cp) : rows(cp.rows), cols(cp.cols), elements(cp.elements) {}

template <class T>
Matrix<T>::~Matrix() {} // atm not neccessary, because there is no custom allocation or something

template <class T>
Matrix<T> &Matrix<T>::operator=(const Matrix<T> &cp)
{
    if (cp.rows != rows && cp.cols != cols)
        throw std::domain_error("matrix op= not of same order");
    for (uint32_t i = 0; i < rows * cols; i++)
        elements[i] = cp.elements[i];
    return *this;
}

template <class T>
void Matrix<T>::range_check(uint32_t i, uint32_t j) const
{
    if (rows <= i)
        throw std::range_error("matrix access row out of range");
    if (cols <= j)
        throw std::range_error("matrix access col out of range");
}

// TODO:
template <class T>
bool Matrix<T>::operator==(const Matrix<T> &A) const
{
    const Matrix<T> &N = *this;

    // check dimension:
    if (N.rows != A.rows | N.cols != A.cols)
        return false;

    // every field must be equal to its equivalent in the other matrix
    for (uint32_t i = 0; i < rows * cols; i++)
        if (N.elements[i] != A.elements[i])
            return false;
    return true;
} // O(rows * cols)

template <class T>
bool Matrix<T>::iszero() const
{
    const Matrix<T> &N = *this;
    for (uint32_t i = 0; i < rows * cols; i++)
        if (N.elements[i] != 0)
            return false;
    return true;
} // O(n)

template <class T>
Matrix<T> &Matrix<T>::operator*=(const T &a)
{
    for (uint32_t i = 0; i < rows * cols; i++)
        this->elements[i] *= a;
    return this;
} // O(n) // not yet tested!

template <class T>
Matrix<T> &Matrix<T>::operator/=(const T &a)
{
    for (uint32_t i = 0; i < rows * cols; i++)
        this->elements[i] /= a;
    return this;
} // O(n) // not yet tested!

template <class T>
Matrix<T> Matrix<T>::operator*(const Matrix<T> &M) const {} // matrix multiplication

template <class T>
Matrix<T> &Matrix<T>::operator+=(const Matrix<T> &M)
{
    for (uint32_t i = 0; i < rows * cols; i++)
        this->elements[i] += M.elements[i];
    return *this;
} // O(n) // not yet tested!

template <class T>
Matrix<T> &Matrix<T>::operator-=(const Matrix<T> &M)
{
    for (uint32_t i = 0; i < rows * cols; i++)
        this->elements[i] -= M.elements[i];
    return *this;
} // O(n) // not yet tested!

template <class T>
Matrix<T> Matrix<T>::minor(unsigned i, unsigned j) const
{
    // returns minor matrix without j-th row  and i-th col

#ifdef RANGE_CHECK
    range_check(i, j);
#endif
    /*T minor_arr[(rows - 1) * (cols - 1)];

    uint32_t counter = 0;
    for (uint32_t n = 0; n < rows * cols; n++)
    {
        if (!((i == n / cols) || (j == n % cols)))
        {
            // pick up in minor matrix:
            minor_arr[counter] = this->elements[rows * j + i];
            counter++;
        }
    }

    Matrix<T> M(this->rows - 1, this->cols - 1, minor_arr);

    return M;*/

    const Matrix<T> &N = *this;
    Matrix<T> _modM(N);

    _modM = _modM.delcol(i);
    _modM = _modM.delrow(j);

    return _modM;

    // if it doesn't work, analyse the following code to fix mine:
    /*
    void minor(float b[100][100],float a[100][100],int i,int n){
	int j,l,h=0,k=0;
	for(l=1;l<n;l++)
		for( j=0;j<n;j++){
			if(j == i)
				continue;
			b[h][k] = a[l][j];
			k++;
			if(k == (n-1)){
				h++;
				k=0;
			}
		}
        */

} // O(n²) // not yet tested! // runtime isnt correct!

template <class T>
T Matrix<T>::minor_det(unsigned i, unsigned j) const
{
#ifdef RANGE_CHECK
    range_check(i, j);
#endif

    const Matrix<T> &N = *this;

    /*Matrix<T> _modM(N);

    _modM = _modM.delcol(i);
    _modM = _modM.delrow(j);*/

    return std::pow(-1, i + j) * N(i, j);
} // O(1) // not yet tested!

template <class T>
T Matrix<T>::det() const
{
    if (rows != cols)
        throw std::domain_error("matrix determinant: incompatible orders");

    const Matrix<T> &N = *this; // N is a square matrix

    if (N.rows == 1)
        return N(0, 0);
    else
        return std::pow(-1, i + j) * (N.delcol(0).delrow(0)).det(); // recursive!
} // Laplace Expansion O(n!) !! -> efficient for max. 5x5 matrices! TODO: Implement a more efficient algorithm for bigger matrices

template <class T>
Matrix<T> Matrix<T>::leftdiv(const Matrix<T> &D) const
{
    const Matrix<T> &N = *this;

    if (N.rows != D.rows)
        throw std::domain_error("matrix divide: incompatible orders");

    Matrix<T> Q(D.cols, N.cols); // quotient matrix

    if (N.cols > 1)
    {
        // break this down for each column of the numerator
        for (uint32_t j = 0; j < Q.cols; j++)
            Q.setcol(j, N.getcol(j).leftdiv(D)); // note: recursive
        return Q;
    }

    // from here on, N.col == 1

    if (D.rows < D.cols)
        throw underdetermined();

    if (D.rows > D.cols)
    {
        bool solution = false;
        for (uint32_t i = 0; i < D.rows; i++)
        {
            Matrix<T> D2 = D.delrow(i); // delete a row from the matrix
            Matrix<T> N2 = N.delrow(i);
            Matrix<T> Q2(Q);
            try
            {
                Q2 = N2.leftdiv(D2);
            }
            catch (underdetermined x)
            {
                continue; // try again with next row
            }
            if (!solution)
            {
                // this is our possible solution
                solution = true;
                Q = Q2;
            }
            else
            {
                // do the solutions agree=
                if (Q != Q2)
                    throw overdetermined();
            }
        }
        if (!solution)
            throw underdetermined();
        return Q;
    }

    // D.rows == D.cols && N.cols == 1
    // use Kramer's Rule
    //
    // D is a square matrix of order N x N
    // N is a matrix of order N x 1

    const T T0(0.0); // additive identity

    if (D.cols <= 3)
    {
        T ddet = D.det();
        if (ddet == T0)
            throw underdetermined();
        for (uint32_t j = 0; j < D.cols; j++)
        {
            Matrix<T> A(D); // make a copy of the D matrix
            // replace column with numerator vector
            A.setcol(j, N);
            Q(j, 0) = A.det() / ddet;
        }
    }
    else
    {
        // this method optimizes the determinant calculations by saving a cofactors used in calculating the denominator determinant.
        kahn_sum<T> sum;
        std::vector<T> cofactors(D.cols); // save cofactors
        for (uint32_t j = 0; j < D.cols; j++)
        {
            T c = D.minor_det(0, j);
            cofactors[j] = c;
            T a = D(0, j);
            if (a != T0)
            {
                a *= c;
                if (j % 2)
                    sum -= a;
                else
                    sum += a;
            }
        }
        T ddet = sum;
        if (ddet == T0)
            throw underdetermined();
        for (uint32_t j = 0; j < D.cols; j++)
        {
            Matrix<T> A(D);
            A.setcol(j, N);
            kahn_sum<T> ndet;
            for (uint32_t k = 0; k < D.cols; k++)
            {
                T a = A(0, k);
                if (a != T0)
                {
                    if (k == j)
                        a *= cofactors[k]; // use previously calculated cofactor
                    else
                        a *= A.minor_det(0, k); // calculate minor's determinant
                    if (k % 2)
                        ndet -= a;
                    else
                        ndet += a;
                }
            }
            Q(j, 0) = T(ndet) / ddet;
        }
    }
    return Q;
}

template <class T>
Matrix<T> Matrix<T>::inverse() const // TODO!
{
    //if (this->det == 0)
    //    throw std::logic_error("singularity! det(M) == 0! matrix cannot be inverted");
    const Matrix<T> &N = *this;

    T _det = N.det();

    if (_det == 0)
        throw std::logic_error("cannot calculate inverse matrix with det() == 0");

    T _num = 1/_det;

    // here: det() prooves also if matrix is squared!


    // find adjoint matrix
    // Matrix<T> adjoint() // use method from leftdiv inkl. kahan summation method!
}

template <class T>
Matrix<T> Matrix<T>::getrow(unsigned i) const
{
#ifdef RANGE_CHECK
    range_check(i, 0);
#endif

    const Matrix<T> &N = *this;

    Matrix<T> _row(1, N.cols); // one-row-matrix

    // adds just elements that are in the i-th row
    uint32_t counter = 0;
    for (uint32_t n = i; n < cols; n++)
        _row.elements[counter++] = N(n, i);

    return _row;
} // O(cols) -> O(n) if N x N Matrix // not yet tested!

template <class T>
Matrix<T> Matrix<T>::getcol(unsigned j) const
{
#ifdef RANGE_CHECK
    range_check(0, j);
#endif

    const Matrix<T> &N = *this;

    Matrix<T> _col(N.rows, 1); // one-col-matrix

    // adds just elements that are in the j-th column
    uint32_t counter = 0;
    for (uint32_t n = j; n < rows; n++)
        _col.elements[counter++] = N(j, n);

    return _col;
} // O(rows) -> O(n) if N x N Matrix // not yet tested!

template <class T>
Matrix<T> Matrix<T>::delcol(unsigned j) const
{
#ifdef RANGE_CHECK
    range_check(0, j);
#endif

    const Matrix<T> &N = *this;

    // creates Matrix<> with number of fields like N only one column less
    Matrix<T> _modM(N.rows, N.cols - 1);
    //T _modarr[N.rows*(N.cols-1)] = { 0 };
    uint32_t counter = 0; // counter for arr

    for (uint32_t n = 0; n < rows * cols; n++)
        if (!(n % cols == j)) // add every element thats NOT in j-th col to Matrix
            _modM.elements[counter++] = N.elements[n];
    // _modarr[counter++] = N.elements[n];

    return _modM;
} // O(n²) // not yet tested!

template <class T>
Matrix<T> Matrix<T>::delrow(unsigned i) const
{
#ifdef RANGE_CHECK
    range_check(i, 0);
#endif

    const Matrix<T> &N = *this;

    // creates Matrix<> with number of fields like N only one row less
    Matrix<T> _modM(N.rows - 1, N.cols);

    uint32_t counter = 0; // counter for arr

    for (uint32_t n = 0; n < rows * cols; n++)
        if (!(n / rows == i)) // add every element thats NOT in the i-th row to Matrix
            _modM.elements[counter++] = N.elements[n];

    return _modM;
} // O(n²) // not yet tested!

template <class T>
Matrix<T> &Matrix<T>::setcol(unsigned j, const Matrix<T> &C)
{
#ifdef RANGE_CHECK
    range_check(0, j);
#endif

    const Matrix<T> &N = *this;
    Matrix<T> _modM(N);

    // check dimensions of C:
    if (C.cols != N.cols || C.rows != 1)
        throw std::domain_error("matrix set col: incompatible matrix C");

    // everythings fine, lets go!
    for (uint32_t n = 0; n < N.rows; n++)
        _modM(j, n) = C(0, n);

    return _modM;
} // O(N.rows) -> O(n) if N is squared matrix // not yet tested!

template <class T>
Matrix<T> &Matrix<T>::setrow(unsigned i, const Matrix<T> &R)
{
#ifdef RANGE_CHECK
    range_check(i, 0);
#endif

    const Matrix<T> &N = *this;
    Matrix<T> _modM(N);

    // check dimensions of R
    if (R.row != N.rows || R.cols != 1)
        throw std::domain_error("matrix set row: incompatible matrix R");

    // everythings fine, lets go!
    for (uint32_t n = 0; n < N.cols; n++)
        _modM(n, i) = R(n, 0);

    return _modM;
}

template <class T>
Matrix<T> Matrix<T>::identity() const
{
    if (rows != cols)
        throw std::domain_error("matrix identity: incompatible orders");

    uint8_t _identity_arr[rows * cols] = {0}; // uint8_t is attempt to minimize using storage ; not sure if it compiles!

    for (uint32_t i = 0; i < rows * cols; i += (rows + 1))
        _identity_arr[i] = 1;

    Matrix<T> _identity(rows, cols, _identity_arr);
    return _identity;
} // O(rows * cols) // not yet tested!

template <class T>
bool Matrix<T>::isidentity() const
{
    return (this->identity() == this);
} // O((rows*cols)^2) // not yet tested!

template <class T>
Matrix<T> Matrix<T>::pow(uint32_t exp) const
{
    const Matrix<T> &N = *this;

    if (rows != cols)
        throw std::domain_error("matrix pow: incompatible orders");

    if (exp == 1)            // termination condition & trivial solution
        return N;            // A^1 = A
    else if (exp == 0)       // for completeness to cover all possibilites. shouldn't slow recursion
        return N.identity(); // A^0 = E
    else
    {                            // A^n = A*A*...(n-times)...*A*A
        return N * pow(exp - 1); // RECURSIVE!
    }
} // O(exp) // not yet tested!

template <class T>
Matrix<T> Matrix<T>::transpose() const
{
    // in-situ is no option, because original matrix will be used farther
    const Matrix<T> &N = *this;

    Matrix<T> _transp(cols, rows, &(this->elements[0])); // this was discussed in Scott Meyers' Effective STL, that you can do &vec[0] to get the address of the first element of an std::vector (https://stackoverflow.com/questions/4289612/getting-array-from-stdvector)

#pragma omp parallel for // helpful?
    for (uint32_t n = 0; n < rows * cols; n++)
    {
        uint32_t i = n / rows;
        uint32_t j = n % cols;
        //_clone.elements[n] = this->elements[rows * j + i];
        _transp.elements[n] = N(i, j); // TESTING!
    }

} // O(rows * cols) // not yet tested!

#endif // matrix.h

// Grobe Struktur und leftdiv Funktion aus:
// Quelle: https://www.drdobbs.com/a-c-matrix-template-class/184403323?pgno=1
