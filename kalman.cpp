// HEADER
#include "kalman.h"

// defines
#define dt 1      // Zeitschritt
#define rounds 10 // Anzahl Iterationen

// debugging
#define plot(x) std::cout << x << std::endl;

// lambdas
const auto sq = [](double i) noexcept -> double { return (i * i); };         // square
const auto cu = [](double i) noexcept -> double { return (sq(i) * i); };     // cubic
const auto bq = [](double i) noexcept -> double { return (sq(i) * sq(i)); }; // biquadrate

// define const 1D Arrays to "convert" them into "2D" matrices
//                   x  y   x'  y'
const double x_[] = {0, 0, 10, 0}; // Systemzustand

const double P_[] = {10, 0, 0, 0,
                     0, 10, 0, 0,
                     0, 0, 10, 0,
                     0, 0, 0, 10}; // Kovarianzmatrix
const double A_[] = {1, 0, dt, 0,
                     0, 1, 0, dt,
                     0, 0, 1, 0,
                     0, 0, 0, 1}; // Dynamikmatrix
const double H_[] = {0, 0, 1, 0,
                     0, 0, 0, 1}; // Messmatrix
const double R_[] = {10, 0,
                     0, 10}; // Messrauschkovarianzmatrix
const double Q_[] = {0.25 * bq(dt), 0.25 * bq(dt), 0.5 * cu(dt), 0.5 * cu(dt),
                     0.25 * bq(dt), 0.25 * bq(dt), 0.5 * cu(dt), 0.5 * cu(dt),
                     0.5 * cu(dt), 0.5 * cu(dt), sq(dt), sq(dt),
                     0.5 * cu(dt), 0.5 * cu(dt), sq(dt), sq(dt)}; // Prozessrauschkovarianzmatrix
const double I_[] = {1, 0, 0, 0,
                     0, 1, 0, 0,
                     0, 0, 1, 0,
                     0, 0, 0, 1}; // Einheitsmatrix

// Beteiligte Matrizen:
Matrix<double> x(1, 4, x_); // Zustandsvektor/Systemzustand (Ist nichts bekannt, kann hier 0 eingetragen werden, bekannte Werte x y x' y' können dem Filter mitgeteilt werden, sofern bekannt)
Matrix<double> P(4, 4, P_); // Unsicherheit/Kovarianzmatrix (Ist man sich zu Beginn über die Zustände ganz sicher, können hier niedrige Werte eingetragen werden. Weiß man zu Beginn nicht, wie die Werte des Zustandsvektor sind, hier sehr große Werte (~ 10^6) initialisieren)
Matrix<double> A(4, 4, A_); // Dynamik (Sagt aus, wohin sich der Zustandsvektor im nächsten Schritt entwickeln wird) (dt!)
Matrix<double> H(4, 2, H_); // Messmatrix (Was wird gemessen und in welchem Verhältnis steht es zum Zustandsvektor)
Matrix<double> R(2, 2, R_); // Messrauschkovarianzmatrix (Messunsicherheit, gibt an, wie verlässlich die Messwerte sind) (Ist Sensor sehr genau, hier kleine Werte einsetzen, ansonsten eher große Werte)
Matrix<double> Q(4, 4, Q_); // Prozessrauschkovarianzmatrix (Rauschen/Störungen zwischen einer Berechnungsiteration) (Lässt sich am einfachsten berechnen, in dem man den Vektor G aufstellt, und diesen nachfolgend mit der angenommenden Standardabweichung für die angreifende Beschleunigung multipliziert, z.B. sigma = 8 m/s^2 :=> Q = G * ~G * sigma² mit G = [0.5*dt², 0.5*dt², dt, dt])
Matrix<double> I(4, 4, I_); // Einheitsmatrix

// Einwirken externer Steuergrößen (z.B. Lenken, Bremsen, etc.) ist über die Steuermatrix B möglich, wird hier aber weggelassen!

Matrix<double> Z(1, 2); // aktuelle Messwerte (sind hier nur x und y, deswegen nur 2 Felder groß!)
Matrix<double> S(2, 2); // Residualkovarianz
Matrix<double> K(2, 4); // Kalman-Gain
Matrix<double> y(1, 2); // Innovation des Zustandsvektors x mit Messmatrix H

int main(void)
{
    /* Exemplarische Vorgehensweise:
        folgende Matrix soll abgebildet werden:
        [ 1 2 3 ]
        [ 4 5 6 ]

        das entspricht:

        Matrix<int> M(3, 2, (const int[]){1, 2, 3, 4, 5, 6});

        Konstruktor:
        erstes Argument: Anzahl Spalten
        zweites Argument: Anzahl Zeilen
        drittes Argument: Zeilenweise eingeben (optional)

        
        Matrix<T> ist nullbasiert:

        M(0,0) liefert erstes Element
    */

    for (size_t i = 0; i < rounds; i++)
    {
        // Prediction
        x = A * x;              // Prädizierter Zustand aus Bisherigem und System
        P = (A * P * (~A)) + Q; // Prädizieren der Kovarianz // set brackets to keep correct order (multiplication)

        // Correction
        // Z =  // Z is getting there new record from measurements (sensors)
        y = Z - (H * x);                // Innovation aus Messwertdifferenz
        S = ((H * P * (~H)) + R);       // Innovationskovarianz // set brackets to keep correct order (multiplication)
        K = (P * (~H) * (S.inverse())); // Filter-Matrix (Kalman-Gain)

        x = x + (K * y);       // aktualisieren des Systemzustands
        P = (I - (K * H)) * P; // aktualisieren der Kovarianz

        plot(x); // Ausgabe des neuen Zustandsvektors
    }
}

// Quelle und gute Erklärung: https://www.cbcity.de/das-kalman-filter-einfach-erklaert-teil-2

// Gute Seiten für maschinelle Matrix Operationen: https://rechneronline.de/lineare-algebra/matrizen.php (unterschiedlich große Matrizen)