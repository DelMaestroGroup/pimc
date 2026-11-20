#ifndef SMS_PHYSICS_H
#define SMS_PHYSICS_H

#include <array>

// Single-microwave-shielded pair potential, in Del's kelvin/angstrom units.
namespace sms {

struct LabParameters {
    double dipole_debye;
    double rabi_mhz;       // Omega/(2*pi), not angular MHz
    double detuning_mhz;   // delta/(2*pi), blue detuning >= 0
    double xi_fraction;    // xi/(pi/4), in [0,1)
};

struct Couplings {
    double detuning_ratio;
    double xi_radians;
    double sin_two_xi;
    double c3_kelvin_angstrom3;
    double c6_kelvin_angstrom6;
};

struct ValueGradient {
    double value;
    std::array<double,3> gradient; // K/angstrom
};

Couplings couplingsFromLab(const LabParameters &parameters);

// Bare pair kernel: no wrapping, cutoff, images or tail correction here.
// At coincident points V=+infinity and the gradient is undefined (NaN).
ValueGradient interaction(const std::array<double,3> &separation,
                          const Couplings &couplings);

} // namespace sms
#endif
