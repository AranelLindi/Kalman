#ifndef matrix_h
#define matrix_h

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
private:
    std::vector<T> elements; // array of elements

public:
    const uint32_t rows; // number of rows
    const uint32_t cols; // number of columns

protected:
    // range check function for matrix access
    void range_check(uint32_t i, uint32_t j) const;

public:
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

public:
    // constructors
    Matrix(uint32_t rows, uint32_t columns, const T *elements = 0);
    Matrix(const Matrix<T> &);
    // destructor
    ~Matrix();

    // assignment
    Matrix<T> &operator=(const Matrix<T> &);

    // comparison
    bool operator==(const Matrix<T> &) const;
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
    Matrix<T> &operator+=(const Matrix<T> &);
    Matrix<T> &operator-=(const Matrix<T> &);
    Matrix<T> operator+(const Matrix<T> &M) const
    {
        return Matrix<T>(*this).operator+=(M);
    }
    Matrix<T> operator-(const Matrix<T> &M) const
    {
        return Matrix<T>(*this).operator-=(M);
    }

    // matrix multiplication
    Matrix<T> operator*(const Matrix<T> &)const;
    Matrix<T> &operator*=(const Matrix<T> &M)
    {
        return *this = *this * M;
    }

    // matrix division
    Matrix<T> leftdiv(const Matrix<T> &) const;
    Matrix<T> rightdiv(const Matri<T> &D) const
    {
        return transpose().leftdiv(D.transpose()).transpose();
    }
    Matrix<T> operator/(const Matrix<T> &D) const
    {
        return rightdiv(D);
    }
    Matrix<T> &operator/=(const Matrix<T>)
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
    Matrix<T> identidy() const;
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

#endif // matrix.h