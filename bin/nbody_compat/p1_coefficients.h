// Stable P1 softened-potential Taylor coefficients for falcON's scalar kernel.
// For Q=r^2+eps^2 and H=eps^2/2, psi=M*(Q^-1/2 + H*Q^-3/2).
// D_n=(-1/r d/dr)^n psi = T_n*[1+(2*n+1)*H/Q],
// T_n=M*(2*n-1)!!*Q^(-n-1/2). Do not form T_(n+1) just to multiply it by H:
// at eps=0.1 pc in kpc units it can overflow float while D_n is finite.
#ifndef OCEN_P1_COEFFICIENTS_H
#define OCEN_P1_COEFFICIENTS_H
#include <cmath>

template<int Order, typename Real>
inline bool ocen_p1_coefficients(Real inverse_q, Real half_eps2, Real* derivative)
{
    const double x = static_cast<double>(inverse_q);
    const double correction = static_cast<double>(half_eps2) * x;
    double term = static_cast<double>(derivative[0]) * std::sqrt(x);
    for (int n = 0; n <= Order; ++n) {
        derivative[n] = static_cast<Real>(term * (1.0 + (2*n+1)*correction));
        if (!std::isfinite(static_cast<double>(derivative[n]))) return false;
        if (n < Order) term *= (2*n+1)*x;
    }
    return true;
}
#endif
