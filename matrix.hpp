/*
 partial source: https://www.drdobbs.com/a-c-matrix-template-class/184403323?pgno=1
*/


#ifndef MATRIX_H
#define MATRIX_H

// Standartbibliothek
#include <vector>    // vector
#include <cstdint>   // uint types
#include <stdexcept> // exceptions
#include <cmath>     // pow

// undefine to disable rang checking
#define RANGE_CHECK

// undefine to disable error calculation
#define ROUNDTOZERO

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
class kahan_sum
{ // implements Kahan Summation method
public:
    kahan_sum() : sum(0.0), cor(0.0) {}
    kahan_sum<T> &operator+=(const T &val)
    {
        T old_sum = sum;
        T next = val - cor;
        cor = ((sum += next) - old_sum) - next;
        return *this;
    }
    kahan_sum<T> &operator-=(const T &val)
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
    Matrix(uint32_t columns, uint32_t rows, const T *elements = 0);
    Matrix(const Matrix<T> &M);
    // destructor
    ~Matrix();

    // assignment
    Matrix<T> &operator=(const Matrix<T> &M);

    // comparison
    bool operator==(const Matrix<T> &M) const;
    bool operator!=(const Matrix<T> &M) const;
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

template <class T> // important: enter elements line-by-line NOT column-by-column
Matrix<T>::Matrix(uint32_t cols, uint32_t rows, const T *elements) : rows(rows), cols(cols), elements(rows * cols, T(0.0))
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
Matrix<T>::~Matrix() {} // atm not neccessary, bc there is no custom allocation or something

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

    // every field must be equal to its equivalent in the other matrix:
    for (uint32_t i = 0; i < rows * cols; i++)
        if (N.elements[i] != A.elements[i])
            return false;
    return true;
} // O(N²) // not yet tested!

template <class T>
bool Matrix<T>::operator!=(const Matrix<T> &A) const
{
    const Matrix<T> &N = *this;

    // check dimension:
    if (N.rows == A.rows | N.cols == A.cols)
        return false;

    // every field must be equal to its equivalent in the other matrix:
    for (uint32_t i = 0; i < rows * cols; i++)
        if (N.elements[i] == A.elements[i])
            return false;

    return true;
} // O(N²) // not yet tested!

template <class T>
bool Matrix<T>::iszero() const
{
    const Matrix<T> &N = *this;
    for (uint32_t i = 0; i < rows * cols; i++)
        if (N.elements[i] != 0)
            return false;
    return true;
} // O(n) // not yet tested!

template <class T>
Matrix<T> &Matrix<T>::operator*=(const T &a)
{
    Matrix<T> &N = *this;

    for (uint32_t i = 0; i < rows * cols; i++)
        N.elements[i] *= a;
    return N;
} // O(n) // not yet tested!

template <class T>
Matrix<T> &Matrix<T>::operator/=(const T &a)
{
    Matrix<T> &N = *this;

    if (a == 0)
        throw std::logic_error("matrix scalar division: divide by zero");

    for (uint32_t i = 0; i < rows * cols; i++)
        N.elements[i] /= a;
    return N;
} // O(n) // not yet tested!

template <class T>
Matrix<T> Matrix<T>::operator*(const Matrix<T> &M) const
{
    const Matrix<T> &N = *this;

    // condition multiplication:
    if (N.cols != M.rows)
        throw std::logic_error("matrix multiplication: incompatible orders");

    Matrix<T> O(M.cols, N.rows);

    for (uint32_t i = 0; i < N.rows; i++)
        for (uint32_t j = 0; j < M.cols; j++)
        {
            kahan_sum<T> sum;

            for (uint32_t k = 0; k < M.rows; k++)
                sum += N(i, k) * M(k, j);

#ifdef ROUNDTOZERO
            // (not sure if this is help and useful... watch it!)
            // intension: entire instruction is evaluated at compile time
            // thus an optimization is guaranteed
            if (constexpr(std::is_same<T, float>))
            {
                if (abs(sum) < 10e-6)
                    sum = 0;
            }
            else if (constexpr(std::is_same<T, double>))
            {
                if (abs(sum) < 10e-12)
                    sum = 0;
            }
            else if (constexpr(std::is_same<T, long double>))
            {
                if (abs(sum) < 10e-16)
                    sum = 0;
            }
            // Anmerkung:
            // Wegen Ungenauigkeit bei Gleitkommazahlen, entstehen bei Multiplikation teilweise ungenaue Ergebnisse.
            // So müsste mathematisch gesehen eine 0 herauskommen, es wird aber eine sehr kleine Zahl (10^(-16) berechnet.
            // Prüfen ob hier mittels eines Verfahrens gerundet werden kann.
#endif

            O(i, j) = sum;
        }

    return O;
} // O(n*m*p) (O(n³)) // better solution: strassen algorithm: O(n^2.808) // matrix multiplication

template <class T>
Matrix<T> &Matrix<T>::operator+=(const Matrix<T> &M)
{
    const Matrix<T> &N = *this;

    if (N.cols != M.cols || N.rows != M.rows)
        throw std::domain_error("matrix addition: incompatible orders");

    // from here on: both matrices have same dimension

    // so simply add component by component
    for (uint32_t i = 0; i < rows * cols; i++)
        this->elements[i] += M.elements[i];
    return *this;
} // O(n) // not yet tested!

template <class T>
Matrix<T> &Matrix<T>::operator-=(const Matrix<T> &M)
{
    const Matrix<T> &N = *this;

    if (N.cols != M.cols || N.rows != M.rows)
        throw std::domain_error("matrix subtraction: incompatible orders");

    // from here on: both matrices have same dimension

    // so simply substract component by component
    for (uint32_t i = 0; i < rows * cols; i++)
        this->elements[i] -= M.elements[i];

    return *this;
} // O(n) // not yet tested!

template <class T>
Matrix<T> Matrix<T>::minor(uint32_t i, uint32_t j) const // returns minor matrix (= original matrix without j-th row  and i-th col)
{
    const Matrix<T> &N = *this;

    // ********************************************
    //Matrix<T> M(N.rows - 1, N.cols - 1);

    // add elements only that are not effected in i-j-range
    /*uint32_t counter = 0;
    for (uint32_t _j = 0; _j < N.cols; _j++)
        for (uint32_t _i = 0; _i < N.rows; _i++)
            if (!(_i == i | _j == j))
                M.elements[counter++] = N(_j, _i);*/

    // code above is just if the one-line-instruction doesn't work (keep for testing purposes!)
    // ********************************************

    // N withouth i-th column and j-th row
    return N.delcol(i).delrow(j);
} // O(N²)

template <class T>
T Matrix<T>::minor_det(uint32_t i, uint32_t j) const // returns det of minor_matrix (max 3x3 matrices!)
{
    const Matrix<T> &N = *this;

    if (N.rows != N.cols)
        throw std::domain_error("minor_det: incompatible orders");

    // from here on: N is square matrix

    if (N.rows <= 3)
    {
        return (((i + j) % 2 ? -1 : 1) * (N.minor(i, j)).det());
    }
    else
    {
        // matrix is bigger than 3x3 and must be reduced:           DO LATER!
        kahan_sum<T> sum;

        for (uint32_t i = 0; i < N.rows; i++)
            sum += N.minor_det(0, i);

        return sum;
    }
} // WC: O(n!), BC: O(1) -> function is more efficient the smaller the matrix is // not yet tested!

template <class T>
T Matrix<T>::det() const //
{
    const Matrix<T> &N = *this;

    if (N.rows != N.cols)
        throw std::domain_error("det: incompatible orders");

    // from here on: N is square matrix

    if (N.rows == 3)
    { // 3x3
        kahan_sum<T> sum;
        // gleiches auch mal ohne kahan probieren und vergleichen!
        sum += N(0, 0) * N(1, 1) * N(2, 2);
        sum += N(0, 1) * N(1, 2) * N(2, 0);
        sum += N(0, 2) * N(1, 0) * N(2, 1);
        sum -= N(0, 2) * N(1, 1) * N(2, 0);
        sum -= N(0, 0) * N(1, 2) * N(2, 1);
        sum -= N(0, 1) * N(1, 0) * N(2, 2);

        return sum;
    }
    else if (N.rows == 2)
    {                                                                         // 2x2
        return N.elements[0] * N.elements[3] - N.elements[2] * N.elements[1]; // das auch mal mit kahan Summe probieren und vergleichen!
    }
    else // 1x1 because all other cases are excluded by the previous code of this function an calling code
    {
        return N.elements[0];
    }
} // O(1) // not yet tested!

template <class T>
Matrix<T> Matrix<T>::leftdiv(const Matrix<T> &D) const
{
    const Matrix<T> &N = *this;

    if (N.rows != D.rows)
        throw std::domain_error("matrix divide: incompatible orders");

    //Matrix<T> Q(D.cols, N.cols); // quotient matrix
    Matrix<T> Q(N.cols, D.cols); // debug!

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
            auto test = A.det() / ddet; // debug
            Q(0, j) = A.det() / ddet;   // debug: Q(0, j)
        }
    }
    else
    {
        // this method optimizes the determinant calculations by saving a cofactors used in calculating the denominator determinant.
        kahan_sum<T> sum;
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
            kahan_sum<T> ndet;
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
    const Matrix<T> &N = *this;

    T _det = N.det();

    if (_det == 0)
        throw std::logic_error("inverse: cannot calculate inverse matrix with det() == 0");

    // here: det() checks also if N is squared!

    // find adjoint matrix
    Matrix<T> adjoint(N.cols, N.rows);

    // calculate and save determinant for every minor matrix
    for (uint32_t j = 0; j < N.cols; j++)
        for (uint32_t i = 0; i < N.rows; i++)
            adjoint(i, j) = N.minor_det(i, j) / _det;

    return adjoint;
} // O(n!) // not yet tested!

template <class T>
Matrix<T> Matrix<T>::getrow(uint32_t i) const
{
    const Matrix<T> &N = *this;

    Matrix<T> M(N.cols, 1); // create one col-matrix (= vector)

    uint32_t counter = 0; // counts cols
    for (uint32_t f = i * N.cols; counter < N.cols; f++)
        M.elements[counter++] = N.elements[f];

    return M;
} // O(cols) -> O(n) if N x N Matrix // tested!

template <class T>
Matrix<T> Matrix<T>::getcol(uint32_t j) const
{
    const Matrix<T> &N = *this;

    Matrix<T> M(1, N.rows); // create one row-matrix

    uint32_t counter = 0; // counts rows
    for (uint32_t f = j; counter < N.rows; f += N.cols)
        M.elements[counter++] = N.elements[f];

    return M;
} // O(rows) -> O(n) if N x N Matrix // tested! erneut testen! hab es aus versehen geändert aber wiederhergestellt. sicher gehen!!!!

template <class T>
Matrix<T> Matrix<T>::delcol(uint32_t j) const
{
    const Matrix<T> &N = *this;

    Matrix<T> M(N.cols - 1, N.rows); // one column less than N

    uint32_t counter = 0; // counts cols
    for (uint32_t f = 0; counter < N.cols - 1; f++)
        if (f != j)
            M.setcol(counter++, N.getcol(f));

    return M;
} // O(n²) // not yet tested! seems to work! more testing!!

template <class T>
Matrix<T> Matrix<T>::delrow(uint32_t i) const
{
    const Matrix<T> &N = *this;

    Matrix<T> M(N.cols, N.rows - 1); // one row less than N

    uint32_t counter = 0; // counts rows
    for (uint32_t f = 0; counter < N.rows - 1; f++)
        if (f != i)
            M.setrow(counter++, N.getrow(f));

    return M;
} // O(n²) // not yet tested! seems to work! more testing!!!

template <class T>
Matrix<T> &Matrix<T>::setcol(uint32_t j, const Matrix<T> &C)
{
    Matrix<T> &N = *this;

    uint32_t counter = 0; // counts rows
    for (uint32_t f = j; counter < N.rows; f += N.cols)
        N.elements[f] = C.elements[counter++];

    return N;
} // O(N.rows) -> O(n) if N is squared matrix // tested!

template <class T>
Matrix<T> &Matrix<T>::setrow(uint32_t i, const Matrix<T> &R)
{
    Matrix<T> &N = *this;

    uint32_t counter = 0; // counts cols
    for (uint32_t f = i * N.cols; counter < N.cols; f++)
        N.elements[f] = R.elements[counter++];

    return N;
} // O(N.cols) // tested!

template <class T>
Matrix<T> Matrix<T>::identity() const
{
    const Matrix<T> &N = *this;

    if (N.rows != N.cols)
        throw std::domain_error("matrix identity: incompatible orders");

    Matrix<T> M(N.cols, N.rows); // create zero-matrix with N dimension

    // from here: M is a square matrix

    for (uint32_t i = 0; i < M.rows; i++)
        M(i, i) = 1;

    return M;
} // O(N) // not yet tested!

template <class T>
bool Matrix<T>::isidentity() const
{
    return (identity() == this);
} // O(N²) // not yet tested!

template <class T>
Matrix<T> Matrix<T>::pow(uint32_t exp) const
{
    const Matrix<T> &N = *this;

    if (rows != cols)
        throw std::domain_error("matrix pow: incompatible orders");

    // use recursion to calc pow:
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
    const Matrix<T> &N = *this;

    Matrix<T> Tr(N.rows, N.cols); // swap number of cols & rows

    // case distinction:
    if (N.rows == 1 || N.cols == 1) // just copy
        for (uint32_t i = 0; i < N.rows * N.cols; i++)
            Tr.elements[i] = N.elements[i];
    else // swap elements (i, j) -> (j, i), except for elements with i==j but costs are negligible
        for (uint32_t i = 0; i < N.cols; i++)
            for (uint32_t j = 0; j < N.rows; j++)
                Tr(j, i) = N(i, j);

    return Tr;
} // O(N²) // not yet tested!

#endif // matrix.hpp