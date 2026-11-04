// InjectedPolKernel.h
//
// Polarized two-body hyperon decay kernel + polarization models for the InjectedPolO2MC afterburner.
// TEST CODE -- see ../README.md. Not an official ALICE tool.
//
// Header-only, no ROOT / Pythia / O2 dependency on purpose, and C++11-compatible (some `pythia8-config --cxxflags`
// add -std=c++11, which overrides an earlier -std=c++17 on the command line):
//   - it is #included by the O2 external generator macro (JIT-compiled by Cling inside o2-sim),
//   - and compiled stand-alone by ../tests/testKernel.cxx (plain g++) for unit tests.
//
// Physics conventions (must stay identical to the analysis, see
// O2Physics PWGLF/Tasks/Strangeness/lambdaJetPolarizationIonsDerived.cxx):
//   * dN/dcos(xi*) = (1 + alpha |P*| cos(xi*)) / 2, xi* = angle between the baryon daughter momentum in the
//     hyperon rest frame and the polarization vector P* (also in the rest frame).
//   * Polarization axes are built from LAB-frame vectors and used directly as quantization axes in the
//     hyperon rest frame. This is exactly what the measured ring observable projects on:
//       R = (3/alpha) < p_hat*_baryon . (j_hat x p_Lambda) / |j_hat x p_Lambda| >
//   * Ring axis t_hat = (j_hat x p_Lambda) / |j_hat x p_Lambda|  (same order of the cross product as the analysis).
//
// Kernel originally from PhD_codes/vorticity_hydro/lambda_pol_toy_model.cxx
// (sample_P_angle_proton_from_cos_xi + Lambda_decay), re-written to:
//   * use a numerically stable inverse CDF (no cancellation for alpha|P| -> 0, no special case needed),
//   * avoid acos/cos/sin round trips (cos xi* is used directly, sin xi* = sqrt(1 - cos^2)),
//   * return BOTH daughters (baryon and meson), exactly back-to-back in the rest frame,
//   * take the parent mass from its own 4-momentum (exact 4-momentum conservation).

#ifndef INJECTEDPOL_KERNEL_H
#define INJECTEDPOL_KERNEL_H

#include <algorithm>
#include <cmath>

namespace injpol
{

//==================================================================================================
// Constants (PDG 2025). Decay parameters MUST match lambdaJetPolarizationIonsDerived.cxx
//==================================================================================================
namespace constants
{
constexpr double kPi = 3.14159265358979323846;
constexpr double kAlphaLambdaToPPi = 0.749;      // Lambda    -> p pi-      (alpha_-),  as in the consumer task
constexpr double kAlphaAntiLambdaToPPi = -0.758; // AntiLambda -> pbar pi+ (alpha_+),  as in the consumer task
constexpr double kAlphaLambdaToNPi0 = 0.670;     // Lambda    -> n pi0      (alpha_0)
constexpr double kAlphaAntiLambdaToNPi0 = -0.675; // AntiLambda -> nbar pi0 (alpha_0bar)
constexpr double kBRLambdaToPPi = 0.641;         // PDG; the two channels are renormalised to 1 in the generator
constexpr double kBRLambdaToNPi0 = 0.359;
} // namespace constants

//==================================================================================================
// Minimal 3-/4-vector helpers (kept trivial so Cling and g++ treat them identically)
//==================================================================================================
struct V3 {
  double x = 0., y = 0., z = 0.;
  V3() = default;
  V3(double a, double b, double c) : x(a), y(b), z(c) {}
};
inline V3 operator+(const V3& a, const V3& b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
inline V3 operator-(const V3& a, const V3& b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
inline V3 operator*(double s, const V3& a) { return {s * a.x, s * a.y, s * a.z}; }
inline double dot(const V3& a, const V3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline V3 cross(const V3& a, const V3& b) { return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x}; }
inline double norm(const V3& a) { return std::sqrt(dot(a, a)); }
inline V3 unit(const V3& a)
{
  const double n = norm(a);
  return n > 0. ? (1. / n) * a : V3{0., 0., 0.};
}
inline double phiOf(const V3& a) { return std::atan2(a.y, a.x); }

struct P4 {
  V3 p;
  double e = 0.;
  double mass() const { return std::sqrt(std::max(0., e * e - dot(p, p))); }
};

/// Boost a 4-vector given in the parent rest frame into the frame where the parent has 4-momentum `parent`.
inline P4 boostFromRest(const P4& vRest, const P4& parent, double parentMass)
{
  const V3 beta = (1. / parent.e) * parent.p;
  const double b2 = dot(beta, beta);
  if (b2 < 1e-30) {
    return vRest;
  }
  const double gamma = parent.e / parentMass;
  const double bp = dot(beta, vRest.p);
  const double g2 = (gamma - 1.) / b2;
  P4 out;
  out.p = vRest.p + (g2 * bp + gamma * vRest.e) * beta;
  out.e = gamma * (vRest.e + bp);
  return out;
}

/// Boost a lab 4-vector into the rest frame of `parent` (inverse of boostFromRest).
inline P4 boostToRest(const P4& vLab, const P4& parent, double parentMass)
{
  P4 inv = parent;
  inv.p = -1. * parent.p;
  return boostFromRest(vLab, inv, parentMass);
}

/// Two-body break-up momentum in the parent rest frame.
inline double twoBodyMomentum(double M, double m1, double m2)
{
  const double a = M * M - (m1 + m2) * (m1 + m2);
  const double b = M * M - (m1 - m2) * (m1 - m2);
  return (a > 0. && b > 0.) ? std::sqrt(a * b) / (2. * M) : 0.;
}

//==================================================================================================
// Angular sampling
//==================================================================================================
/// Inverse CDF of dN/dc = (1 + A c)/2 on c in [-1, 1], |A| <= 1, u uniform in [0, 1).
/// Same solution as lambda_pol_toy_model.cxx, c = (-1 + sqrt((1-A)^2 + 4 A u)) / A, but rationalised:
///   c = (A - 2 + 4u) / (1 + sqrt((1-A)^2 + 4 A u)),
/// which is exact for A = 0 (c = 2u - 1) and has no catastrophic cancellation for small |A|.
inline double clampUnit(double x) { return std::max(-1., std::min(1., x)); } // C++11-safe (no std::clamp)

inline double sampleCosXi(double A, double u)
{
  A = clampUnit(A);
  const double s = std::sqrt(std::max(0., (1. - A) * (1. - A) + 4. * A * u));
  return clampUnit((A - 2. + 4. * u) / (1. + s));
}

struct TwoBodyDecay {
  P4 baryon;     // lab frame (p, pbar, n, nbar)
  P4 meson;      // lab frame (pi-, pi+, pi0)
  V3 baryonDirRest; // unit vector of the baryon momentum in the parent rest frame (for checks)
};

/// Polarized two-body decay H -> B + M.
///  parent : lab 4-momentum of the hyperon
///  mB, mM : daughter masses
///  polRest: polarization vector P* (rest frame, lab-oriented axes), |P*| <= 1. Zero -> isotropic.
///  alpha  : decay parameter of this channel (sign included: negative for antihyperons)
///  u1, u2 : two independent uniform random numbers in [0, 1)
inline TwoBodyDecay decayPolarized(const P4& parent, double mB, double mM, const V3& polRest, double alpha, double u1, double u2)
{
  const double M = parent.mass();
  const double q = twoBodyMomentum(M, mB, mM);
  const double pMag = std::min(norm(polRest), 1.);

  V3 n{0., 0., 1.};
  double c;
  if (pMag > 0.) {
    n = (1. / norm(polRest)) * polRest;
    c = sampleCosXi(alpha * pMag, u1);
  } else {
    c = 2. * u1 - 1.; // isotropic
  }
  const double s = std::sqrt(std::max(0., 1. - c * c));
  const double phi = 2. * constants::kPi * u2;

  // Orthonormal basis (n, a, b). The arbitrary choice of `a` is harmless because phi is uniform.
  const V3 t = (std::abs(n.z) < 0.99) ? V3{0., 0., 1.} : V3{0., 1., 0.};
  const V3 a = unit(cross(t, n));
  const V3 b = cross(n, a);
  const V3 dir = c * n + (s * std::cos(phi)) * a + (s * std::sin(phi)) * b;

  P4 bRest, mRest;
  bRest.p = q * dir;
  bRest.e = std::sqrt(mB * mB + q * q);
  mRest.p = -1. * bRest.p;
  mRest.e = std::sqrt(mM * mM + q * q);

  TwoBodyDecay out;
  out.baryon = boostFromRest(bRest, parent, M);
  out.meson = boostFromRest(mRest, parent, M);
  out.baryonDirRest = dir;
  return out;
}

/// Lab-frame decay displacement and time for a proper decay length l (same length units as l, c = 1).
inline void decayDisplacement(const P4& parent, double properLength, V3& dx, double& dt)
{
  const double M = parent.mass();
  dx = (properLength / M) * parent.p; // beta*gamma*l * p_hat
  dt = (parent.e / M) * properLength; // gamma*l  (time in length units, e.g. mm/c)
}

//==================================================================================================
// Polarization models
//==================================================================================================
struct EventAxes {
  bool hasLeadP = false; // leading charged particle above threshold
  V3 leadP;              // its momentum (lab)
  bool hasJet = false;   // leading accepted anti-kT jet above threshold
  V3 jet;                // its momentum (lab)
  double psiRP = 0.;     // (random) reaction-plane angle used by global/longitudinal models
};

enum RingProxy { kRingLeadP = 0,
                 kRingJet = 1 };

struct ModelConfig {
  double ring = 0.;         // P0 for P* = P0 t_hat
  int ringProxy = kRingJet; // which axis defines t_hat
  double global = 0.;       // P_G along -y' (normal to the reaction plane, see below)
  double longitudinal = 0.; // A in P*_z = A sin(n (phi_Lambda - psi))
  int longHarmonic = 2;     // n (2 = usual elliptic-flow-driven longitudinal polarization)
  double helicity = 0.;     // P_H along p_hat_Lambda
  double antiSign = 1.;     // multiply the total P* of antihyperons by this (+1: same sign as hyperons)
};

/// Total injected polarization vector P* for a hyperon of lab momentum pL. Returns |P*| <= 1
/// (rescaled if the sum of the models exceeds 1; `clamped` reports it).
inline V3 polarization(const ModelConfig& cfg, const EventAxes& ax, const V3& pL, bool isAnti, bool* clamped = nullptr)
{
  V3 P{0., 0., 0.};
  if (cfg.ring != 0.) {
    const bool useJet = (cfg.ringProxy == kRingJet);
    if ((useJet && ax.hasJet) || (!useJet && ax.hasLeadP)) {
      const V3 jHat = unit(useJet ? ax.jet : ax.leadP);
      const V3 t = cross(jHat, pL); // same order as the analysis: leadingJetUnitVec.Cross(lambdaLike3Vec)
      const double tn = norm(t);
      if (tn > 1e-12) {
        P = P + (cfg.ring / tn) * t;
      }
    }
  }
  if (cfg.global != 0.) {
    // Reaction-plane frame: x' = (cos psi, sin psi, 0), y' = (-sin psi, cos psi, 0).
    // Global polarization along -y' (system angular momentum direction for this convention).
    P = P + cfg.global * V3{std::sin(ax.psiRP), -std::cos(ax.psiRP), 0.};
  }
  if (cfg.longitudinal != 0.) {
    P = P + V3{0., 0., cfg.longitudinal * std::sin(cfg.longHarmonic * (phiOf(pL) - ax.psiRP))};
  }
  if (cfg.helicity != 0.) {
    P = P + cfg.helicity * unit(pL);
  }
  if (isAnti) {
    P = cfg.antiSign * P;
  }
  const double n = norm(P);
  if (clamped) {
    *clamped = (n > 1.);
  }
  return n > 1. ? (1. / n) * P : P;
}

} // namespace injpol

#endif // INJECTEDPOL_KERNEL_H
