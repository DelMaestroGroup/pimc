# Single-microwave-shielded molecules

`--interaction=sms` adds an opt-in anisotropic effective pair potential. It
does not change moves, worm weights, estimators, boundary conditions or state
serialization. There is no new third-party dependency.

## Model and conventions

The elliptically polarized effective potential is Eq. (6.3) in
[X.-Y. Chen, *Microwave-shielded ultracold polar molecules*, Chapter 6](https://edoc.ub.uni-muenchen.de/32631/1/Chen_Xing-Yan.pdf).
For the circular limit and derivation of the effective shielding model, see
[F. Deng et al., Phys. Rev. Lett. 130, 183001 (2023)](https://doi.org/10.1103/PhysRevLett.130.183001).
These are effective-model references, not validation of every experimental
parameter accepted by this implementation. In particular, the approximation
must be assessed against coupled-channel calculations for the chosen dressing
parameters; a mathematical domain check is not a physical validity guarantee.

The fixed polarization axes are x, y and z, with the microwave propagating along
z. Let \(s=\sin(2\xi)\), \(f=s\cos(2\phi)\), and \(\eta=\delta/\Omega\).
The implemented pair energy is

\[
V(\mathbf r)=\frac{C_3}{r^3}
 [3\cos^2\theta-1+3f\sin^2\theta]
+\frac{C_6}{r^6}\sin^2\theta
 [1-f^2+(1-f)^2\cos^2\theta].
\]

With SI dipole moment \(d\), \(D=d^2/(4\pi\epsilon_0)\), the coefficients are

\[
C_{3,\mathrm{SI}}=\frac{D}{12(1+\eta^2)},\qquad
C_{6,\mathrm{SI}}=\frac{D^2}{8\hbar\Omega(1+\eta^2)^{3/2}}.
\]

Divide by \(k_B\,10^{-30}\) and \(k_B\,10^{-60}\), respectively, to obtain
K angstrom^3 and K angstrom^6. The implementation uses exact SI h and k_B,
and \(\hbar\Omega=h(\Omega/2\pi)\), avoiding a hidden extra factor of 2*pi.

The code evaluates an equivalent Cartesian expression. With \(R=x^2+y^2+z^2\)
and \(Q=z^2+s(x^2-y^2)\), the r^-6 term is
\(C_6(R-Q)(R+Q)/R^5\). Both factors are evaluated directly as quadratic
polynomials, without dividing by x^2+y^2 or subtracting nearly equal R and Q.
The analytic gradient is in K/angstrom. Coincident particles have infinite
pair energy; their gradient is undefined and returned as NaN.

## Inputs

All four SMS parameters are required; none is inferred from a molecule name.

| Option | Unit / domain |
| --- | --- |
| `sms_dipole_debye` | Permanent molecular dipole, Debye; positive |
| `sms_omega_mhz` | Omega/(2*pi), MHz; positive |
| `sms_delta_mhz` | Blue-detuning magnitude delta/(2*pi), MHz; non-negative |
| `sms_xi_fraction` | xi/(pi/4), dimensionless; 0 <= value < 1 |

For example, 13.5 degrees means `--sms_xi_fraction=0.3`, **not** 13.5.
The linear-polarization endpoint xi=pi/4 is excluded because the shielding
core disappears along an attractive direction. Red-detuned states are not
supported by this interface.

All other Del units are unchanged: `temperature` is in K, `mass` in atomic mass
units and cell dimensions/`potential_cutoff` in angstroms. Thus 7 nK is `7e-9`,
and a 20 micrometer box has `--size=200000`. No new nondimensional simulation
units or global temperature/length aliases are introduced.

## Boundaries and truncation

Del supplies the minimum-image displacement in periodic directions. This PR
does **not** add extra images, self-image energies, Ewald sums, anisotropic tail
corrections or new geometry. It uses the existing prism or cylinder boundary
conditions.

The adapter explicitly returns zero for r >= `potential_cutoff`, so full-pair
energy measurements and nearest-neighbor updates use the same truncated
Hamiltonian. Inside the cutoff it returns the unshifted pair potential. The
cutoff is sharp: the gradient is only the ordinary derivative away from this
surface, not a distributional cutoff correction. The existing default cutoff
is retained (the cell's last side length). Set it deliberately when using
anisotropic cells. A sufficiently large cutoff includes every minimum-image
pair, but **still does not implement an infinite dipolar periodic sum**.

Long-range/thermodynamic-limit applications require a separately specified
summation prescription and convergence study. Results of this PR must not be
treated as interchangeable with finite-image runs in the downstream fork.

## Supported execution

This first contribution supports 3D CPU finite-temperature PIMC using
`--action=primitive`. Pair-product, higher-order actions, GPU and PIGS uses of
SMS are not enabled. The existing radial `virial` estimator is rejected for
SMS; use `energy`. The implementation deliberately does not advertise an
unimplemented Laplacian as zero. Non-SMS configurations keep their behavior.

Minimal small-system *smoke test*, not a converged physical calculation:

Build note for upstream revision `184c800`: GCC 13 exposes a pre-existing
missing `<complex>` include in `include/estimator.h`. Until that header is
fixed upstream, configure with `-DCMAKE_CXX_FLAGS="-include complex"` for GCC
or Clang. This workaround is a build setting, not an SMS source modification;
apply it equally to both builds when running the regression comparison.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DNDIM=3 -DGPU_BACKEND=none
cmake --build build -j2
ctest --test-dir build --output-on-failure

# Run in a fresh directory to avoid mixing outputs.
mkdir sms-smoke
cd sms-smoke
../build/pimc.e --canonical --geometry=prism --size=10000 \
  --number_particles=4 --mass=155.8 --temperature=7e-9 \
  --number_time_slices=16 --update_length=4 --action=primitive \
  --interaction=sms --sms_dipole_debye=4.6 \
  --sms_omega_mhz=4 --sms_delta_mhz=22 --sms_xi_fraction=0.3 \
  --external=free --number_eq_steps=0 --bin_size=2 --number_bins_stored=2 \
  --update="center of mass" --update=bisection --estimator=energy
```

The SMS unit test has no dependency beyond C++17's standard library. Optional
CLI integration tests use Python's standard library. They exercise prism and
cylinder runs, required/invalid options, metadata, cutoff-consistent energy
and text-state continuation. A build of unmodified upstream can additionally
be compared using:

```bash
python3 tests/sms_integration.py build/pimc.e --baseline /path/to/upstream/pimc.e
```

This compares same-seed free/harmonic final states (including stored RNG) and
energy rows, and reads upstream text checkpoints as initial states for SMS.
Short smoke tests and analytic tests do not establish many-body equilibration,
imaginary-time convergence, or physical accuracy of the effective model.
