#include "potential.h"
#include "factory_potential.h"

#include <cmath>
#include <stdexcept>

#if NDIM == 3
REGISTER_INTERACTION_POTENTIAL("sms", SMSPotential, GET_SETUP(),
    setup.params["sms_dipole_debye"].as<double>(),
    setup.params["sms_omega_mhz"].as<double>(),
    setup.params["sms_delta_mhz"].as<double>(),
    setup.params["sms_xi_fraction"].as<double>())

SMSPotential::SMSPotential(double dipole, double omega, double delta, double xi) :
    couplings(sms::couplingsFromLab({dipole,omega,delta,xi})),
    cutoffSquared(constants()->rc2()) {
    if (!std::isfinite(cutoffSquared) || cutoffSquared <= 0)
        throw std::invalid_argument("SMS cutoff must be finite and positive");
}

double SMSPotential::V(const dVec &r) {
    // Apply the same cutoff in all-pair energy and nearest-neighbor updates.
    if (dot(r,r) >= cutoffSquared)
        return 0;
    return sms::interaction(r,couplings).value;
}

dVec SMSPotential::gradV(const dVec &r) {
    if (dot(r,r) >= cutoffSquared)
        return dVec{};
    return sms::interaction(r,couplings).gradient;
}

double SMSPotential::grad2V(const dVec &) {
    throw std::logic_error("SMS Laplacian is not implemented; use primitive action and energy estimator");
}
#endif
