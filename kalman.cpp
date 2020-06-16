// HEADER
#include "kalman.h"

// defines
#define dt 0      // Zeitschritt
#define rounds 10 // iterations to make

// lambdas
const auto bq = [](double i) -> double { return (i * i * i * i); }; // biquadrate
const auto cu = [](double i) -> double { return (i * i * i); };     // cubic
const auto sq = [](double i) -> double { return (i * i); };         // square

// define const 1D Arrays to "convert" them into 2D matrixes
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

// initialize matrixes
auto *x = new Matrix<double>(4, 1, x_);
auto *P = new Matrix<double>(4, 4, P_);
auto *A = new Matrix<double>(4, 4, A_);
auto *H = new Matrix<double>(2, 4, H_);
auto *R = new Matrix<double>(2, 2, R_);
auto *Q = new Matrix<double>(4, 4, Q_);
auto *I = new Matrix<double>(4, 4, I_); // Identity matrix

auto *Z = new Matrix<double>(4, 4); // exact size still unknown! TODO (gets filling later...)
auto *S = new Matrix<double>(4, 4); // exact size still unknown! TODO (gets filling later...)
auto *K = new Matrix<double>(4, 4); // exact size still unknown! TODO (gets filling later...)
auto *y = new Matrix<double>(4, 4); // exact size still unknown! TODO (gets filling later...)

int main(void)
{
    // doesnt work at the moment:
    for (size_t i = 0; i < rounds; i++)
    {
        // Prediction
        *x = (*A) * (*x);                      // Prädizierter Zustand aus Bisherigem und System
        *P = (((*A) * (*P)) * (~(*A))) + (*Q); // Prädizieren der Kovarianz // set brackets to keep correkt order (multiplication)

        // Correction
        // Z =  // Z is getting there new record from measurements (sensors)
        *y = (*Z) - ((*H) * (*x));               // Innovation aus Messwertdifferenz
        *S = ((((*H) * (*P)) * (~(*H))) + (*R)); // Innovationskovarianz // set brackets to keep correct order (multiplication)
        *K = ((*P) * (~(*H)) * (*S).inverse());  // Filter-Matrix (Kalman-Gain)

        *x = *x + ((*K) * (*y));            // aktualisieren des Systemzustands
        *P = ((*I) - ((*K) * (*H))) * (*P); // aktualisieren der Kovarianz
    }
}