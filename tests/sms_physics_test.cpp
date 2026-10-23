#include "sms_physics.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>

namespace {
void require(bool condition, const char *message) {
    if (!condition) {
        std::cerr << message << '\n';
        std::exit(EXIT_FAILURE);
    }
}
void near(double a, double b, double tolerance, double scale = 0) {
    require(std::isfinite(a) && std::isfinite(b), "non-finite comparison");
    if (std::abs(a-b) > tolerance*std::max({std::abs(a),std::abs(b),scale})) {
        std::cerr << "actual=" << a << " expected=" << b << '\n';
        require(false, "comparison failed");
    }
}
template<class F> void invalid(F operation) {
    try { operation(); }
    catch (const std::invalid_argument &) { return; }
    require(false, "invalid parameter was accepted");
}

// Independent spherical expression; used away from the polar coordinate axis.
double spherical(const std::array<double,3> &r, const sms::Couplings &c) {
    const double length = std::hypot(r[0],r[1],r[2]);
    const double phi = std::atan2(r[1],r[0]);
    const double ct = r[2]/length, st2 = 1-ct*ct;
    const double f = c.sin_two_xi*std::cos(2*phi);
    return c.c3_kelvin_angstrom3*(3*ct*ct-1+3*f*st2)/std::pow(length,3) +
        c.c6_kelvin_angstrom6*st2*(1-f*f+(1-f)*(1-f)*ct*ct)/std::pow(length,6);
}
void gradient(const std::array<double,3> &r, const sms::Couplings &c) {
    const double length = std::hypot(r[0],r[1],r[2]);
    const double scale = c.c3_kelvin_angstrom3/std::pow(length,4) +
                         c.c6_kelvin_angstrom6/std::pow(length,7);
    const auto v = sms::interaction(r,c);
    const double step = length*1e-5;
    for (int a = 0; a < 3; ++a) {
        auto plus = r, minus = r;
        plus[a] += step;
        minus[a] -= step;
        near(v.gradient[a],(sms::interaction(plus,c).value-
             sms::interaction(minus,c).value)/(2*step),2e-7,scale);
    }
}
} // namespace

int main() {
    const sms::LabParameters lab{4.6,4,22,0.3};
    const auto c = sms::couplingsFromLab(lab);
    near(c.detuning_ratio,5.5,1e-15);
    near(c.xi_radians,0.23561944901923448,1e-15);
    near(c.c3_kelvin_angstrom3,408.69644992609824,2e-14);
    // Exact SI h differs slightly from the previously rounded hbar constant.
    near(c.c6_kelvin_angstrom6,87552088116.64389,2e-14);
    auto doubled = lab;
    doubled.dipole_debye *= 2;
    const auto d = sms::couplingsFromLab(doubled);
    near(d.c3_kelvin_angstrom3,4*c.c3_kelvin_angstrom3,1e-14);
    near(d.c6_kelvin_angstrom6,16*c.c6_kelvin_angstrom6,1e-14);
    doubled = lab;
    doubled.rabi_mhz *= 2;
    doubled.detuning_mhz *= 2;
    const auto twiceFrequency = sms::couplingsFromLab(doubled);
    near(twiceFrequency.c3_kelvin_angstrom3,c.c3_kelvin_angstrom3,1e-14);
    near(twiceFrequency.c6_kelvin_angstrom6,c.c6_kelvin_angstrom6/2,1e-14);

    std::mt19937 generator(12031);
    std::uniform_real_distribution<double> position(-3000,3000);
    for (double fraction : {0.0,0.3,0.8,0.999}) {
        auto parameters = lab;
        parameters.xi_fraction = fraction;
        const auto coupling = sms::couplingsFromLab(parameters);
        for (int n = 0; n < 200; ++n) {
            const std::array<double,3> r{position(generator),position(generator),position(generator)};
            const std::array<double,3> reverse{-r[0],-r[1],-r[2]};
            const auto v = sms::interaction(r,coupling), w = sms::interaction(reverse,coupling);
            const double length = std::hypot(r[0],r[1],r[2]);
            const double scale = coupling.c3_kelvin_angstrom3/std::pow(length,3) +
                                 coupling.c6_kelvin_angstrom6/std::pow(length,6);
            near(v.value,spherical(r,coupling),3e-13,scale);
            near(v.value,w.value,1e-15,scale);
            for (int a = 0; a < 3; ++a)
                near(v.gradient[a],-w.gradient[a],1e-15,scale/length);
            gradient(r,coupling);
        }
        // Exact axes and arbitrarily close to the polar axis: no rho^-2 branch.
        for (const std::array<double,3> r : {std::array<double,3>{2000,0,0},
                {0,2000,0},{0,0,2000},{0,0,-2000},{1e-10,-1e-10,2000}})
            gradient(r,coupling);
        const auto axis = sms::interaction({0,0,2000},coupling);
        near(axis.value,2*coupling.c3_kelvin_angstrom3/std::pow(2000.,3),1e-14);
        near(axis.gradient[2],-6*coupling.c3_kelvin_angstrom3/std::pow(2000.,4),1e-14);
        require(axis.gradient[0] == 0 && axis.gradient[1] == 0,"polar gradient");
    }
    const auto origin = sms::interaction({0,0,0},c);
    require(std::isinf(origin.value) && origin.value > 0,"coincident pair not rejected");
    require(std::isnan(origin.gradient[0]),"undefined coincident gradient not signaled");
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    for (double value : {nan,inf,-inf}) {
        invalid([&] { sms::couplingsFromLab({value,4,22,0.3}); });
        invalid([&] { sms::couplingsFromLab({4.6,value,22,0.3}); });
        invalid([&] { sms::couplingsFromLab({4.6,4,value,0.3}); });
        invalid([&] { sms::couplingsFromLab({4.6,4,22,value}); });
        invalid([&] { sms::interaction({value,0,1},c); });
    }
    invalid([] { sms::couplingsFromLab({0,4,22,0.3}); });
    invalid([] { sms::couplingsFromLab({4.6,0,22,0.3}); });
    invalid([] { sms::couplingsFromLab({4.6,4,-22,0.3}); });
    invalid([] { sms::couplingsFromLab({4.6,4,22,-0.1}); });
    invalid([] { sms::couplingsFromLab({4.6,4,22,1}); });
    invalid([] { sms::couplingsFromLab({1e300,4,22,0.3}); });
    invalid([] { sms::couplingsFromLab({1e-300,4,22,0.3}); });
    require(sms::couplingsFromLab({4.6,4,0,0}).c3_kelvin_angstrom3 > 0,"zero detuning");
    std::cout << "SMS unit, spherical-reference, gradient, symmetry and input tests passed\n";
}
