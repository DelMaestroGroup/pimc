#include "sms_physics.h"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace sms {
namespace {
constexpr double pi = 3.141592653589793238462643383279502884;
constexpr double boltzmann = 1.380649e-23;       // J/K
constexpr double planck = 6.62607015e-34;        // J s
constexpr double epsilon0 = 8.8541878128e-12;    // F/m
constexpr double debye = 3.33564e-30;           // C m

void requireFinite(double value, const char *name) {
    if (!std::isfinite(value))
        throw std::invalid_argument(std::string(name) + " must be finite");
}
} // namespace

Couplings couplingsFromLab(const LabParameters &p) {
    requireFinite(p.dipole_debye, "SMS dipole");
    requireFinite(p.rabi_mhz, "SMS Rabi frequency");
    requireFinite(p.detuning_mhz, "SMS detuning");
    requireFinite(p.xi_fraction, "SMS xi fraction");
    if (p.dipole_debye <= 0 || p.rabi_mhz <= 0 || p.detuning_mhz < 0 ||
            p.xi_fraction < 0 || p.xi_fraction >= 1)
        throw std::invalid_argument("SMS requires dipole>0, Omega>0, delta>=0 and 0<=xi_fraction<1");

    const double dipole = p.dipole_debye*debye;
    const double dipolar = dipole*dipole/(4*pi*epsilon0);
    const double ratio = p.detuning_mhz/p.rabi_mhz;
    const double dressing = 1 + ratio*ratio;
    // hbar*Omega = h*(Omega/2pi). SI m^3 and m^6 become A^3 and A^6.
    const double c3 = dipolar/(12*dressing)/(boltzmann*1e-30);
    const double c6 = dipolar*dipolar/
        (8*std::pow(dressing,1.5)*planck*p.rabi_mhz*1e6)/(boltzmann*1e-60);
    if (!std::isfinite(c3) || !std::isfinite(c6) || c3 <= 0 || c6 <= 0)
        throw std::invalid_argument("SMS parameters produce non-finite or underflowed couplings");
    const double xi = p.xi_fraction*pi/4;
    return {ratio,xi,std::sin(2*xi),c3,c6};
}

ValueGradient interaction(const std::array<double,3> &r, const Couplings &c) {
    for (double x : r)
        requireFinite(x, "SMS separation");
    const double x2 = r[0]*r[0], y2 = r[1]*r[1], z2 = r[2]*r[2];
    const double r2 = x2+y2+z2;
    if (r2 == 0) {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        return {std::numeric_limits<double>::infinity(),{nan,nan,nan}};
    }
    const double s = c.sin_two_xi;
    const std::array<double,3> a3{-1+3*s,-1-3*s,2};
    const double q3 = a3[0]*x2+a3[1]*y2+2*z2;
    // Factor r^4-[z^2+s(x^2-y^2)]^2 to avoid polar-axis division and
    // cancellation. This is exactly the spherical SMS r^-6 angular factor.
    const std::array<double,3> am{1-s,1+s,0}, ap{1+s,1-s,2};
    const double qm = am[0]*x2+am[1]*y2;
    const double qp = ap[0]*x2+ap[1]*y2+2*z2;
    const double inv2 = 1/r2;
    const double inv5 = inv2*inv2*std::sqrt(inv2);
    const double inv10 = inv5*inv5;
    ValueGradient result{c.c3_kelvin_angstrom3*q3*inv5 +
                         c.c6_kelvin_angstrom6*qm*qp*inv10,{}};
    for (int axis = 0; axis < 3; ++axis) {
        result.gradient[axis] = c.c3_kelvin_angstrom3*r[axis]*inv5*
            (2*a3[axis]-5*q3*inv2) + c.c6_kelvin_angstrom6*r[axis]*inv10*
            (2*(am[axis]*qp+ap[axis]*qm)-10*qm*qp*inv2);
    }
    return result;
}
} // namespace sms
