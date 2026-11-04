// testKernel.cxx -- stand-alone unit/closure tests of generator/InjectedPolKernel.h (no ROOT, no Pythia).
//
// Build & run:
//   g++ -O2 -std=c++17 -I../generator testKernel.cxx -o testKernel.exe && ./testKernel.exe
// (normally done by run/runTests.sh, which builds into ${OUT_BASE}/bin and keeps the repository clean)
//
// Exit code 0 if all checks pass. Every check compares an estimator with its expectation in units of its
// statistical uncertainty (pull), or with a tight absolute tolerance for exact identities.

#include "InjectedPolKernel.h"

#include <cstdio>
#include <random>
#include <vector>

using namespace injpol;

static int gFailures = 0;

static void checkPull(const char* name, double value, double expected, double sigma, double maxPull = 4.)
{
  const double pull = sigma > 0. ? (value - expected) / sigma : 0.;
  const bool ok = std::abs(pull) < maxPull;
  std::printf("[%s] %-60s value = %+.5f  expected = %+.5f  sigma = %.5f  pull = %+.2f\n", ok ? " OK " : "FAIL", name, value, expected, sigma, pull);
  if (!ok) {
    ++gFailures;
  }
}

static void checkAbs(const char* name, double value, double expected, double tol)
{
  const bool ok = std::abs(value - expected) < tol;
  std::printf("[%s] %-60s value = %+.3e  expected = %+.3e  tol = %.1e\n", ok ? " OK " : "FAIL", name, value, expected, tol);
  if (!ok) {
    ++gFailures;
  }
}

// Running mean and its standard error
struct Mean {
  double s = 0., s2 = 0.;
  long n = 0;
  void add(double x)
  {
    s += x;
    s2 += x * x;
    ++n;
  }
  double mean() const { return s / n; }
  double err() const { return std::sqrt(std::max(0., s2 / n - mean() * mean()) / n); }
};

// Masses (GeV/c^2), PDG
constexpr double mL = 1.115683, mP = 0.93827208816, mPi = 0.13957039, mN = 0.93956542052, mPi0 = 0.1349768;

int main()
{
  std::mt19937_64 rng(12345);
  std::uniform_real_distribution<double> U(0., 1.);

  //------------------------------------------------------------------------------------------------
  // 1) Inverse CDF: <cos xi> = A/3, and agreement with the original toy-model formula
  //------------------------------------------------------------------------------------------------
  std::printf("\n=== 1) sampleCosXi ===\n");
  for (double A : {-0.9, -0.758, -0.3, 0., 1e-9, 0.3, 0.749, 1.0}) {
    Mean m;
    for (int i = 0; i < 2000000; ++i) {
      m.add(sampleCosXi(A, U(rng)));
    }
    char name[128];
    std::snprintf(name, sizeof(name), "<cos xi> = A/3 for A = %+.3g", A);
    checkPull(name, m.mean(), A / 3., m.err());
  }
  {
    // identity with the original (non-rationalised) formula where the latter is numerically safe
    double maxDiff = 0.;
    for (double A : {-0.9, -0.5, 0.2, 0.749}) {
      for (int i = 0; i <= 1000; ++i) {
        const double u = i / 1000.;
        const double orig = (-1.0 + std::sqrt((1.0 - A) * (1.0 - A) + 4.0 * A * u)) / A;
        maxDiff = std::max(maxDiff, std::abs(orig - sampleCosXi(A, u)));
      }
    }
    checkAbs("max |stable - original toy-model inverse CDF|", maxDiff, 0., 1e-12);
    checkAbs("A = 1e-12 reduces to 2u - 1 (u = 0.3)", sampleCosXi(1e-12, 0.3), -0.4, 1e-11);
  }

  //------------------------------------------------------------------------------------------------
  // 2) Exact kinematics: 4-momentum conservation, masses, rest-frame direction
  //------------------------------------------------------------------------------------------------
  std::printf("\n=== 2) decay kinematics ===\n");
  {
    double maxDP = 0., maxDE = 0., maxDm = 0., maxDdir = 0.;
    for (int i = 0; i < 100000; ++i) {
      const double pt = 0.2 + 10. * U(rng), eta = -1.5 + 3. * U(rng), phi = 2. * constants::kPi * U(rng);
      P4 L;
      L.p = {pt * std::cos(phi), pt * std::sin(phi), pt * std::sinh(eta)};
      L.e = std::sqrt(mL * mL + dot(L.p, L.p));
      const V3 pol = 0.8 * unit(V3{U(rng) - 0.5, U(rng) - 0.5, U(rng) - 0.5});
      const auto d = decayPolarized(L, mP, mPi, pol, constants::kAlphaLambdaToPPi, U(rng), U(rng));
      const V3 dp = (d.baryon.p + d.meson.p) - L.p;
      maxDP = std::max(maxDP, norm(dp));
      maxDE = std::max(maxDE, std::abs(d.baryon.e + d.meson.e - L.e));
      maxDm = std::max({maxDm, std::abs(d.baryon.mass() - mP), std::abs(d.meson.mass() - mPi)});
      const P4 bRest = boostToRest(d.baryon, L, L.mass());
      maxDdir = std::max(maxDdir, norm(unit(bRest.p) - d.baryonDirRest));
    }
    checkAbs("max |p_B + p_M - p_Lambda| (GeV)", maxDP, 0., 1e-9);
    checkAbs("max |E_B + E_M - E_Lambda| (GeV)", maxDE, 0., 1e-9);
    checkAbs("max |m_daughter - m_PDG| (GeV)", maxDm, 0., 1e-6);
    checkAbs("max |boostToRest(baryon) dir - generated dir|", maxDdir, 0., 1e-7);
  }

  //------------------------------------------------------------------------------------------------
  // 3) Closure of each polarization model with the analysis-style estimators (truth level)
  //------------------------------------------------------------------------------------------------
  std::printf("\n=== 3) model closure with analysis estimators ===\n");
  const int nLambda = 4000000;
  struct Case {
    const char* name;
    ModelConfig cfg;
    bool anti;
    bool nPi0;
  };
  std::vector<Case> cases;
  {
    ModelConfig c;
    cases.push_back({"none", c, false, false});
    c.ring = 0.30;
    cases.push_back({"ring P0=0.30, Lambda -> p pi-", c, false, false});
    cases.push_back({"ring P0=0.30, AntiLambda -> pbar pi+", c, true, false});
    cases.push_back({"ring P0=0.30, Lambda -> n pi0", c, false, true});
    ModelConfig g;
    g.global = 0.20;
    cases.push_back({"global P_G=0.20", g, false, false});
    ModelConfig l;
    l.longitudinal = 0.40;
    cases.push_back({"longitudinal A=0.40, n=2", l, false, false});
    ModelConfig h;
    h.helicity = -0.25;
    cases.push_back({"helicity P_H=-0.25", h, false, false});
  }

  for (const auto& cs : cases) {
    Mean ring, glob, lng, hel;
    for (int i = 0; i < nLambda; ++i) {
      EventAxes ax;
      ax.hasJet = ax.hasLeadP = true;
      const double jphi = 2. * constants::kPi * U(rng), jeta = -0.5 + U(rng);
      ax.jet = {std::cos(jphi), std::sin(jphi), std::sinh(jeta)};
      ax.leadP = ax.jet;
      ax.psiRP = 2. * constants::kPi * U(rng);

      const double pt = 0.5 + 4.5 * U(rng), eta = -0.8 + 1.6 * U(rng), phi = 2. * constants::kPi * U(rng);
      P4 L;
      L.p = {pt * std::cos(phi), pt * std::sin(phi), pt * std::sinh(eta)};
      L.e = std::sqrt(mL * mL + dot(L.p, L.p));

      const V3 pol = polarization(cs.cfg, ax, L.p, cs.anti);
      double alpha;
      TwoBodyDecay d;
      if (cs.nPi0) {
        alpha = cs.anti ? constants::kAlphaAntiLambdaToNPi0 : constants::kAlphaLambdaToNPi0;
        d = decayPolarized(L, mN, mPi0, pol, alpha, U(rng), U(rng));
      } else {
        alpha = cs.anti ? constants::kAlphaAntiLambdaToPPi : constants::kAlphaLambdaToPPi;
        d = decayPolarized(L, mP, mPi, pol, alpha, U(rng), U(rng));
      }

      // --- Analysis-style estimators: only lab quantities + boost of the reconstructed baryon ---
      const V3 pStarHat = unit(boostToRest(d.baryon, L, L.mass()).p);
      const V3 cr = cross(unit(ax.jet), L.p); // leadingJetUnitVec.Cross(lambdaLike3Vec)
      ring.add(3. / alpha * dot(pStarHat, cr) / norm(cr));
      glob.add(3. / alpha * dot(pStarHat, V3{std::sin(ax.psiRP), -std::cos(ax.psiRP), 0.}));
      lng.add(3. / alpha * pStarHat.z * std::sin(2. * (phiOf(L.p) - ax.psiRP)));
      hel.add(3. / alpha * dot(pStarHat, unit(L.p)));
    }
    std::printf("-- %s\n", cs.name);
    char n1[128];
    std::snprintf(n1, sizeof(n1), "  ring R  [%s]", cs.name);
    checkPull(n1, ring.mean(), cs.cfg.ring, ring.err());
    std::snprintf(n1, sizeof(n1), "  global  [%s]", cs.name);
    checkPull(n1, glob.mean(), cs.cfg.global, glob.err());
    std::snprintf(n1, sizeof(n1), "  long <P_z sin2dphi> = A/2  [%s]", cs.name);
    checkPull(n1, lng.mean(), cs.cfg.longitudinal / 2., lng.err());
    std::snprintf(n1, sizeof(n1), "  helicity  [%s]", cs.name);
    checkPull(n1, hel.mean(), cs.cfg.helicity, hel.err());
  }

  //------------------------------------------------------------------------------------------------
  // 4) Decay length: <proper length> recovered; lab displacement along p_hat
  //------------------------------------------------------------------------------------------------
  std::printf("\n=== 4) decay displacement ===\n");
  {
    const double ctau = 78.9; // mm
    std::exponential_distribution<double> expo(1. / ctau);
    Mean lproper;
    double maxPerp = 0.;
    for (int i = 0; i < 1000000; ++i) {
      P4 L;
      L.p = {1.2 * U(rng), -0.7 * U(rng), 2. * (U(rng) - 0.5)};
      L.e = std::sqrt(mL * mL + dot(L.p, L.p));
      const double l = expo(rng);
      V3 dx;
      double dt;
      decayDisplacement(L, l, dx, dt);
      // proper length = |dx| * m / |p|, and dt^2 - |dx|^2 = l^2 (c = 1)
      lproper.add(norm(dx) * L.mass() / norm(L.p));
      maxPerp = std::max(maxPerp, std::abs(std::sqrt(std::max(0., dt * dt - dot(dx, dx))) - l) / l);
    }
    checkPull("<proper decay length> (mm)", lproper.mean(), ctau, lproper.err());
    checkAbs("max relative |sqrt(dt^2 - dx^2) - l| / l", maxPerp, 0., 1e-9);
  }

  std::printf("\n%s: %d failure(s)\n", gFailures ? "FAILED" : "ALL PASSED", gFailures);
  return gFailures ? 1 : 0;
}
