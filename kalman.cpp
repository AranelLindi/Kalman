// HEADER
#include "kalman.h"

// defines
#define dt 0      // Zeitschritt
#define rounds 10 // iterations to make

// lambdas (constexpr in this case needs C++17! only possible until variables will be used!)
constexpr auto sq = [](double i) noexcept -> double { return (i * i); };         // square
constexpr auto cu = [](double i) noexcept -> double { return (sq(i) * i); };     // cubic
constexpr auto bq = [](double i) noexcept -> double { return (sq(i) * sq(i)); }; // biquadrate

// define const 1D Arrays to "convert" them into "2D" matrices
//                   x  y   x'  y'
const double x_[] = {0, 0, 10, 0};
const double P_[] = {10, 0, 0, 0,
                     0, 10, 0, 0,
                     0, 0, 10, 0,
                     0, 0, 0, 10};
const double A_[] = {1, 0, dt, 0,
                     0, 1, 0, dt,
                     0, 0, 1, 0,
                     0, 0, 0, 1};
const double H_[] = {0, 0, 1, 0,
                     0, 0, 0, 1};
const double R_[] = {10, 0,
                     0, 10};
const double Q_[] = {0.25 * bq(dt), 0.25 * bq(dt), 0.5 * cu(dt), 0.5 * cu(dt),
                     0.25 * bq(dt), 0.25 * bq(dt), 0.5 * cu(dt), 0.5 * cu(dt),
                     0.5 * cu(dt), 0.5 * cu(dt), sq(dt), sq(dt),
                     0.5 * cu(dt), 0.5 * cu(dt), sq(dt), sq(dt)};
const double I_[] = {1, 0, 0, 0,
                     0, 1, 0, 0,
                     0, 0, 1, 0,
                     0, 0, 0, 1};

// initialize matrices
Matrix<double> x(4, 1, x_);
Matrix<double> P(4, 4, P_);
Matrix<double> A(4, 4, A_);
Matrix<double> H(2, 4, H_);
Matrix<double> R(2, 2, R_);
Matrix<double> Q(4, 4, Q_);
Matrix<double> I(4, 4, I_); // Identity matrix

Matrix<double> Z(4, 4); // exact size still unknown! TODO (gets filling later...)
Matrix<double> S(4, 4); // exact size still unknown! TODO (gets filling later...)
Matrix<double> K(4, 4); // exact size still unknown! TODO (gets filling later...)
Matrix<double> y(4, 4); // exact size still unknown! TODO (gets filling later...)

int main(void)
{
    // doesnt work at the moment:
    for (size_t i = 0; i < rounds; i++)
    {
        // Prediction
        x = A * x;                // Prädizierter Zustand aus Bisherigem und System
        P = ((A * P) * (~A)) + Q; // Prädizieren der Kovarianz // set brackets to keep correkt order (multiplication)

        // Correction
        // Z =  // Z is getting there new record from measurements (sensors)
        y = Z - (H * x);              // Innovation aus Messwertdifferenz
        S = (((H * P) * (~H)) + R);   // Innovationskovarianz // set brackets to keep correct order (multiplication)
        K = (P * (~H) * S.inverse()); // Filter-Matrix (Kalman-Gain)

        x = x + (K * y);       // aktualisieren des Systemzustands
        P = (I - (K * H)) * P; // aktualisieren der Kovarianz
    }

    // sizeof(Matrix<double>) == 32 Bytes

    //Matrix<int> i(2,2, (const int[]){0,1,2,3}); // passing anonymous array as argument (C++11)
}