/// HEADER
#include "kalman.h"

/// Makros
#define dt 1     // Zeitschritt
#define rounds 6 // Anzahl Iterationen

// debugging
#define plot(x) std::cout << x << std::endl;

//// Lambdas
const auto sq = [](double i) constexpr noexcept -> double { return (i * i); };         // square
const auto cu = [](double i) constexpr noexcept -> double { return (sq(i) * i); };     // cubic
const auto bq = [](double i) constexpr noexcept -> double { return (sq(i) * sq(i)); }; // biquadrate

/// Konstante Arrays, die bequem dargestellt später die Matrizen befüllen:
//                   x  y   x'  y'
const double x_[] = {0, 0, 10, 0}; // Systemzustand (zu Beginn; wird durch Filter laufend aktualisiert)


const double P_[] = {10, 0, 0, 0,
                     0, 10, 0, 0,
                     0, 0, 10, 0,
                     0, 0, 0, 10}; // Kovarianzmatrix
const double A_[] = {1, 0, dt, 0,  /* x = x' * dt + x */
                     0, 1, 0, dt,  /* y = y' * dt + y */
                     0, 0, 1, 0,   /* x' = x' */
                     0, 0, 0, 1};  /* y' = y' */    // Dynamikmatrix
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
Matrix<double> x(1, 4, x_);       // Zustandsvektor/Systemzustand (Ist nichts bekannt, kann hier 0 eingetragen werden, bekannte Werte x y x' y' können dem Filter mitgeteilt werden, sofern bekannt)
Matrix<double> P(4, 4, P_);       // Unsicherheit/Kovarianzmatrix (Ist man sich zu Beginn über die Zustände ganz sicher, können hier niedrige Werte eingetragen werden. Weiß man zu Beginn nicht, wie die Werte des Zustandsvektor sind, hier mit sehr großen Werte (~ 10^6) initialisieren)
const Matrix<double> A(4, 4, A_); // Dynamik (Sagt aus, wohin sich der Zustandsvektor im nächsten Schritt entwickeln wird) (dt!)
const Matrix<double> H(4, 2, H_); // Messmatrix (Was wird gemessen und in welchem Verhältnis steht es zum Zustandsvektor)
const Matrix<double> R(2, 2, R_); // Messrauschkovarianzmatrix (Messunsicherheit, gibt an, wie verlässlich die Messwerte sind) (Ist Sensor sehr genau, hier kleine Werte einsetzen, ansonsten eher große Werte)
const Matrix<double> Q(4, 4, Q_); // Prozessrauschkovarianzmatrix (Rauschen/Störungen zwischen einer Berechnungsiteration) (Lässt sich am einfachsten berechnen, in dem man den Vektor G aufstellt, und diesen nachfolgend mit der angenommenden Standardabweichung für die angreifende Beschleunigung multipliziert, z.B. sigma = 8 m/s^2 :=> Q = G * ~G * sigma² mit G = [0.5*dt², 0.5*dt², dt, dt])
const Matrix<double> I(4, 4, I_); // Einheitsmatrix
// Konstante Matrizen können durchaus auch im Betrieb geändert werden und sind hier nur aus Effizienzgründen const!

/* 
 Zur Wahl von Q und R:
   - Dynamisch: reagiert schnell auf Änderungen -> große Werte in Q und kleine in R
   - Glättung: filtert Rauschen sehr gut weg -> kleine Werte in Q und große in R
*/

// Einwirken externer Steuergrößen (z.B. Lenken, Bremsen, etc.) ist über die Steuermatrix B möglich, wird hier aber weggelassen.

// Beteiligte Matrizen (werden im Filterprozess geändert):
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

        Entspricht:

        Matrix<int> M(3, 2, (const int[]){1, 2, 3, 4, 5, 6});

        Konstruktor:
        1. Arg: Anzahl Spalten
        2. Arg: Anzahl Zeilen
        3. Arg: Zeilenweise eingeben (optional, wird sonst mit Null initialisiert)

        Matrix<T> ist nullbasiert:

        M(0,0) liefert erstes Element
    */

    // Plot Ausgangslage:
    plot(x);

    for (uint32_t i{0}; i < rounds; i++)
    {
        // 1.) Vorhersage (Prediction):
        x = A * x;              // Prädizierter Zustand aus Bisherigem und System
        P = (A * P * (~A)) + Q; // Prädizieren der Kovarianz // set brackets to keep correct order (multiplication)

        // 2.) Korrektur (Correction):
        //  ********************
        //  hier: neue Messwerte abfragen und an Z übergeben
        //  ********************

        y = Z - (H * x);                // Innovation aus Messwertdifferenz
        S = ((H * P * (~H)) + R);       // Innovationskovarianz
        K = (P * (~H) * (S.inverse())); // Filter-Matrix (Kalman-Gain)

        // 3.) Zustand updaten
        x = x + (K * y);       // aktualisieren des Systemzustands
        P = (I - (K * H)) * P; // aktualisieren der Kovarianz

        // 4.) Neuer Zustand
        plot(x); // Ausgabe des neuen Zustandsvektors auf der Konsole
    }
}

// Quelle und gute Erklärung: https://www.cbcity.de/das-kalman-filter-einfach-erklaert-teil-2

// Gute Seiten für maschinelle Matrix Operationen: https://rechneronline.de/lineare-algebra/matrizen.php (unterschiedlich große Matrizen)
