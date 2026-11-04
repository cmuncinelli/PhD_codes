// standalonePythiaClosure.cxx -- validation step 0 of InjectedPolO2MC (see ../DESIGN.md, section 7).
//
// Runs the afterburner logic of generator_pythia8_injectedPol.C on plain Pythia8 events, WITHOUT O2/Geant4:
// same Config file, same generator-level axes and pre-selection (InjectedPolPythia.h), same per-Lambda
// decay (injpol::decayLambda). Reports
//   - pre-selection efficiency and Pythia time per accepted event (to budget productions),
//   - Lambdas per accepted event and fraction of decays inside a typical V0 radius window,
//   - truth-level estimators (ring R with the configured proxy, global, longitudinal, helicity) computed with
//     analysis-style formulas, to be compared with the injected amplitudes.
//
// Build (needs Pythia8 only; e.g. inside the O2/O2Physics environment, or your own Pythia install):
//   g++ -O2 -I../generator standalonePythiaClosure.cxx -o standalonePythiaClosure $(pythia8-config --cxxflags --libs) -std=c++17
// (-std=c++17 LAST: pythia8-config may add an older -std, and the last one wins.)
// Run:
//   ./standalonePythiaClosure ../generator/pythia8_pp_13600.cfg ../generator/configs/ring_leadP_0p30.cfg 2000 [seed]

#include "InjectedPolKernel.h"
#include "InjectedPolPythia.h"

#include "Pythia8/Pythia.h"

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <string>

namespace
{
struct Mean {
  double s = 0., s2 = 0.;
  long n = 0;
  void add(double x)
  {
    s += x;
    s2 += x * x;
    ++n;
  }
  double mean() const { return n ? s / n : 0.; }
  double err() const { return n > 1 ? std::sqrt(std::max(0., s2 / n - mean() * mean()) / n) : 0.; }
};
void report(const char* name, const Mean& m, double expected)
{
  const double pull = m.err() > 0. ? (m.mean() - expected) / m.err() : 0.;
  std::printf("  %-44s %+.4f +- %.4f   (expected %+.4f, pull %+.2f, N = %ld)\n", name, m.mean(), m.err(), expected, pull, m.n);
}
} // namespace

int main(int argc, char** argv)
{
  if (argc < 4) {
    std::fprintf(stderr, "usage: %s pythia.cfg injpol.cfg nAcceptedEvents [seed]\n", argv[0]);
    return 1;
  }
  const std::string pyCfg = argv[1], injCfg = argv[2];
  const long nWanted = std::atol(argv[3]);
  const int seed = argc > 4 ? std::atoi(argv[4]) : 12345;

  injpol::Config cfg;
  if (!cfg.readFile(injCfg)) {
    return 1;
  }
  cfg.print(std::cout);

  Pythia8::Pythia pythia;
  pythia.readFile(pyCfg);
  pythia.readString("3122:mayDecay = off");
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = " + std::to_string(seed));
  pythia.readString("Next:numberCount = 0");
  if (!pythia.init()) {
    return 1;
  }
  Pythia8::Event work;
  work.init("InjectedPol work event", &pythia.particleData);
  std::unique_ptr<Pythia8::SlowJet> slowJet(injpol::makeSlowJet(cfg));
  const injpol::DecayConsts dc = injpol::makeDecayConsts(cfg, pythia.particleData);
  std::printf("c*tau = %.2f mm, fraction p pi = %.4f\n", dc.ctauMM, dc.fracPPi);

  long nTries = 0, nAccepted = 0, nLambda = 0, nPPi = 0, nInRadius = 0, nClamped = 0, nWithJet = 0;
  Mean ring, glob, lng, hel;
  const auto t0 = std::chrono::steady_clock::now();

  while (nAccepted < nWanted) {
    injpol::AxesInfo info;
    injpol::EventAxes ax;
    // --- same loop as GeneratorPythia8InjectedPol::generateEvent() ---
    while (true) {
      if (++nTries > cfg.maxTries * nWanted) {
        std::fprintf(stderr, "too many tries\n");
        return 1;
      }
      if (!pythia.next()) {
        continue;
      }
      ax = injpol::computeAxes(pythia.event, cfg, *slowJet, work, info);
      if (injpol::preselect(pythia.event, cfg, ax)) {
        break;
      }
    }
    ax.psiRP = cfg.randomPsiRP ? 2. * injpol::constants::kPi * pythia.rndm.flat() : 0.;
    ++nAccepted;
    nWithJet += ax.hasJet;

    const bool useJet = (cfg.model.ringProxy == injpol::kRingJet);
    const bool hasAxis = useJet ? ax.hasJet : ax.hasLeadP;
    const injpol::V3 axis = useJet ? ax.jet : ax.leadP;

    // --- same per-Lambda treatment as GeneratorPythia8InjectedPol::injectPolarizedDecays() ---
    for (int i = 0; i < pythia.event.size(); ++i) {
      const auto& p = pythia.event[i];
      if (!p.isFinal() || std::abs(p.id()) != 3122) {
        continue;
      }
      ++nLambda;
      const bool anti = p.id() < 0;
      injpol::P4 pL;
      pL.p = {p.px(), p.py(), p.pz()};
      pL.e = p.e();
      const injpol::LambdaDecay d = injpol::decayLambda(cfg, ax, pL, anti, dc, pythia.rndm);
      nClamped += d.clamped;
      const double rDecayCM = 0.1 * std::hypot(p.xProd() + d.dx.x, p.yProd() + d.dx.y);
      nInRadius += (rDecayCM > 1.0 && rDecayCM < 30.); // illustrative V0 window (producer: v0radius > 1 cm; lifetime cut 30 cm)
      if (!d.chargedChannel) {
        continue;
      }
      ++nPPi;
      // analysis-style estimators (truth level, full acceptance)
      const injpol::V3 pStarHat = injpol::unit(injpol::boostToRest(d.dec.baryon, pL, pL.mass()).p);
      if (hasAxis) {
        const injpol::V3 cr = injpol::cross(injpol::unit(axis), pL.p);
        if (injpol::norm(cr) > 1e-12) {
          ring.add(3. / d.alpha * injpol::dot(pStarHat, cr) / injpol::norm(cr));
        }
      }
      glob.add(3. / d.alpha * injpol::dot(pStarHat, injpol::V3{std::sin(ax.psiRP), -std::cos(ax.psiRP), 0.}));
      lng.add(3. / d.alpha * pStarHat.z * std::sin(cfg.model.longHarmonic * (injpol::phiOf(pL.p) - ax.psiRP)));
      hel.add(3. / d.alpha * injpol::dot(pStarHat, injpol::unit(pL.p)));
    }
  }
  const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

  std::printf("\n=== summary ===\n");
  std::printf("accepted %ld / %ld Pythia events (efficiency %.3e), %.3f s Pythia+afterburner per accepted event\n",
              nAccepted, nTries, double(nAccepted) / nTries, secs / nAccepted);
  std::printf("accepted events with a fiducial jet: %ld; (anti)Lambdas: %ld (%.2f / event), p pi: %ld, 1 < R_decay < 30 cm: %.3f, |P|>1 clamped: %ld\n",
              nWithJet, nLambda, double(nLambda) / nAccepted, nPPi, nLambda ? double(nInRadius) / nLambda : 0., nClamped);
  std::printf("truth-level estimators (p pi channel, full acceptance):\n");
  report("ring R (configured proxy)", ring, cfg.model.ring);
  report("global <3/a p*.(-y')>", glob, cfg.model.global);
  report("longitudinal <3/a p*_z sin(n dphi)> = A/2", lng, cfg.model.longitudinal / 2.);
  report("helicity <3/a p*.p_hat>", hel, cfg.model.helicity);
  std::printf("(cross-talk between models is expected to be zero only on average over the acceptance;\n"
              " a non-zero ring R in a global/longitudinal run is a genuine geometric leakage worth studying.)\n");
  return 0;
}
