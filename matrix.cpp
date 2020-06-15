#include "matrix.h"

template <class T>
Matrix<T>::Matrix(uint32_t rows, unsigned cols, const T *elements = 0) : rows(rows), cols(cols), elements(rows * cols, T(0.0))
{
    if (rows == 0 | cols == 0)
        throw std::range_error("attempt to create a degenerate matrix");

    // initialize from array
    if (elements)
        for (uint32_t i = 0; i < rows * cols; i++)
            this->element[i] = element[i];
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
// operator==(const matrix<T>& A) const, iszero() const,
// operator*=(const T& a), operator/=(const T& a), operator-() const,
// operator+() const, operator*(const T& a, const matrix<T>& M),
// operator+=(const matrix<T>& M), operator-=(const matrix<T>& M),
// minor(unsigned i, unsigned j) const,
// minor_det(unsigned i, unsigned j) const, det() const -- mb
// operator*( const matrix<T>& B) const ... mb

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
        for (uint32_t j = 0; j < Q.cols, j++)
            Q.setcol(j, N.getcol(j).leftdiv(D)); // note: recursive
        return Q;
    }

    // from here on, N.col == 1

    if (D.rows < D.cols)
        throw underdetermined();

    if (D.rows > D.cols)
    {
        bool solution = false;
        for (uint32_t i = 0; j < D.rows; i++)
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
        vector<T> cofactors(D.cols); // save cofactors
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

// TODO:
// inverse() const, getrow(unsigned i) const,
// getcol(unsigned j) const, delcol(unsigned j) const,
// delrow(unsigned i) const, setcol(unsigned j, const matrix<T>& C)
// setrow(unsigned i, const matrix<T>& R), identity() const,
// isidentity() const, pow(int exp) const,
// pow(const matrix<T>& M, int exp) --mb
