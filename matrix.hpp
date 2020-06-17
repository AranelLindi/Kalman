#ifndef MATRIX_H
#define MATRIX_H

// Quelle: https://www.drdobbs.com/a-c-matrix-template-class/184403323?pgno=1

// Standartbibliothek
#include <vector>    // vector
#include <cstdint>   // uint...
#include <stdexcept> // exceptions

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
Matrix<T>::~Matrix() {}

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

    // every field must be equal to its equivalent
    for (uint32_t i = 0; i < rows * cols; i++)
        if (N.elements[i] != A.elements[i])
            return false;
    return true;
}

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

/*template <class T>
Matrix<T> Matrix<T>::operator-(const Matrix<T>& M) const {} // Argument selbstständig hinzugefügt, richtig? */
// already defined in template class

/*template <class T>
Matrix<T> Matrix<T>::operator+(const Matrix<T>& M) const {} */
// already defined in template class

/*template <class T>
Matrix<T> Matrix<T>::operator*(const T& a) const {} // scalar multiplication */
// already defined in template class

template <class T>
Matrix<T> Matrix<T>::operator*(const Matrix<T> &M) const {} // matrix multiplication

template <class T>
Matrix<T> &Matrix<T>::operator+=(const Matrix<T> &M)
{
    for (uint32_t i = 0; i < rows * cols; i++)
        this->elements[i] += M.elements[i];
    return this;
} // O(n) // not yet tested!

template <class T>
Matrix<T> &Matrix<T>::operator-=(const Matrix<T> &M)
{
    for (uint32_t i = 0; i < rows * cols; i++)
        this->elements[i] -= M.elements[i];
    return this;
} // O(n) // not yet tested!

template <class T>
Matrix<T> Matrix<T>::minor(unsigned i, unsigned j) const {}

template <class T>
T Matrix<T>::minor_det(unsigned i, unsigned j) const {}

template <class T>
T Matrix<T>::det() const {}

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
        // this method optimizes the determinant calculations
        // by saving a cofactors used in calculating the
        // denominator determinant.
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
    

}

template <class T>
Matrix<T> Matrix<T>::getrow(unsigned i) const {}

template <class T>
Matrix<T> Matrix<T>::getcol(unsigned j) const {}

template <class T>
Matrix<T> Matrix<T>::delcol(unsigned j) const {}

template <class T>
Matrix<T> Matrix<T>::delrow(unsigned i) const {}

template <class T>
Matrix<T> &Matrix<T>::setcol(unsigned j, const Matrix<T> &C) {}

template <class T>
Matrix<T> &Matrix<T>::setrow(unsigned i, const Matrix<T> &R) {}

template <class T>
Matrix<T> Matrix<T>::identity() const
{
    if (rows != cols)
        throw std::logic_error("calc identity requires N x N - matrix (squared matrix)");

    uint8_t _identity_arr[rows * cols] = {0}; // uint8_t is attempt to minimize using storage ; not sure if it compiles!

    for (uint32_t i = 0; i < rows * cols; i += (rows + 1))
        _identity_arr[i] = 1;

    Matrix<T> _identity(rows, cols, _identity_arr);
    return _identity;
} // O(rows) / O(cols) // not yet tested!

template <class T>
bool Matrix<T>::isidentity() const { return 0; } // return-dummy for compiler

template <class T>
Matrix<T> Matrix<T>::pow(uint32_t exp) const
{
    if (rows != cols)
        throw std::logic_error("calc pow requires N x N - matrix (squared matrix)");

    if (exp == 0)
        return this->identity(); // A^0 = E
    else if (exp == 1)
        return this; // A^1 = A
    else
    {                                  // A^n = A*A*A*...(n -times)...*A*A
        return (*this) * pow(exp - 1); // recursive!
    }
} // O(exp) // not yet tested!

template <class T>
Matrix<T> Matrix<T>::transpose() const
{
    // in-situ is no option, because original matrix will be used farther
    Matrix<T> _clone(cols, rows, &(this->elements[0])); // this was discussed in Scott Meyers' Effective STL, that you can do &vec[0] to get the address of the first element of an std::vector (https://stackoverflow.com/questions/4289612/getting-array-from-stdvector)

#pragma omp parallel for // helpful?
    for (uint32_t n = 0; n < rows * cols; n++)
    {
        uint32_t i = n / cols;
        uint32_t j = n % cols;
        _clone.elements[n] = this->elements[rows * j + i];
    }
    return _clone;
} // O(rows * cols) // not yet tested!

#endif // matrix.h