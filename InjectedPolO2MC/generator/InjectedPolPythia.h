// InjectedPolPythia.h
//
// Pythia8-dependent pieces of the InjectedPolO2MC afterburner, shared by
//   - generator_pythia8_injectedPol.C (O2 external generator, JIT-compiled by Cling inside o2-sim), and
//   - ../tests/standalonePythiaClosure.cxx (plain Pythia8 executable, no O2/ROOT).
// TEST CODE -- see ../README.md. Not an official ALICE tool.
//
// Contents:
//   * Config: all user settings, read from a simple "key = value" text file (unknown keys are fatal).
//   * computeAxes(): generator-level jet proxies with the same definitions as the reconstruction-level ones of
//     O2Physics PWGLF/TableProducer/Strangeness/lambdaJetPolarizationIons.cxx (LambdaPolRing-v4):
//       - candidates: final-state charged particles, pT > 0.15 GeV/c, |eta| < 0.9, pion-mass hypothesis;
//       - leading particle: highest-pT candidate, valid if pT > 2 GeV/c;
//       - jets: anti-kT R = 0.4, E-scheme, |eta_jet| < 0.9 - R, no background subtraction (pp),
//         pT > 8 GeV/c by default: a LOOSE generator-level threshold for charged-only jets in pp, kept below the
//         analysis threshold so that the generator-level selection does not sculpt the reconstructed jet sample.
//     Jets are ONLY computed when something needs them (ring model with the jet proxy, or presel.requireJet).
//     Generator-vs-reconstructed jet comparisons are NOT done here: they belong in O2Physics, with particle-level
//     jets found on McParticles (PWGJE/JetFinders/jetFinderMCPCharged.cxx practice).
//     Jets are found with Pythia8::SlowJet, which by default delegates to fjcore (the FastJet core shipped
//     inside Pythia), so anti-kT results are those of FastJet without needing an external FastJet library.
//   * preselect(): cheap generator-level event filter applied BEFORE Geant4.

#ifndef INJECTEDPOL_PYTHIA_H
#define INJECTEDPOL_PYTHIA_H

#include "InjectedPolKernel.h"

#include "Pythia8/Pythia.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>

namespace injpol
{

struct Config {
  // --- polarization models (see InjectedPolKernel.h) ---
  ModelConfig model;

  // --- decay ---
  double ctauMM = -1.;                           // Lambda proper decay length in mm; <= 0: take Pythia's tau0 (78.9 mm, = Geant4 value)
  double brPPi = constants::kBRLambdaToPPi;      // renormalised with brNPi0
  double brNPi0 = constants::kBRLambdaToNPi0;
  bool polarizeNPi0 = true;                      // false: n pi0 channel decayed isotropically

  // --- generator-level jet proxies (defaults = LambdaPolRing-v4 producer defaults) ---
  double trackPtMin = 0.15;
  double trackEtaMax = 0.9;
  bool pionMassHypothesis = true;
  double leadPtMin = 2.0;
  double jetR = 0.4;
  double jetPtMin = 8.0;       // generator-level threshold (charged jets, pp); loose on purpose

  // --- event pre-selection (before Geant4) ---
  int minLambdas = 1;           // number of Lambda + AntiLambda required in the acceptance below
  double lambdaRapMax = 1.0;    // |y| (producer's rapidityCut default)
  double lambdaPtMin = 0.0;
  bool requireLeadP = true;
  bool requireJet = false;
  long maxTries = 100000000;    // safety net against infinite loops

  /// Jet finding is only needed to define the ring axis with the jet proxy, or to pre-select on jets.
  bool needJets() const { return requireJet || (model.ring != 0. && model.ringProxy == kRingJet); }

  // --- reaction plane for global/longitudinal models ---
  bool randomPsiRP = true;      // false: psi_RP = 0 in every event

  int verbose = 0;

  /// Parse "key = value" lines; '#' starts a comment. Returns false (and prints) on unknown keys or bad values.
  bool readFile(const std::string& path)
  {
    std::ifstream in(path);
    if (!in) {
      std::cerr << "[InjectedPol] ERROR: cannot open config file '" << path << "'" << std::endl;
      return false;
    }
    std::string line;
    int lineNo = 0;
    bool ok = true;
    while (std::getline(in, line)) {
      ++lineNo;
      const auto hash = line.find('#');
      if (hash != std::string::npos) {
        line.erase(hash);
      }
      const auto eq = line.find('=');
      if (eq == std::string::npos) {
        if (line.find_first_not_of(" \t\r") != std::string::npos) {
          std::cerr << "[InjectedPol] ERROR: " << path << ":" << lineNo << ": expected 'key = value'" << std::endl;
          ok = false;
        }
        continue;
      }
      std::string key = trim(line.substr(0, eq)), val = trim(line.substr(eq + 1));
      if (!set(key, val)) {
        std::cerr << "[InjectedPol] ERROR: " << path << ":" << lineNo << ": unknown key or bad value '" << key << " = " << val << "'" << std::endl;
        ok = false;
      }
    }
    return ok;
  }

  bool set(const std::string& key, const std::string& val)
  {
    std::istringstream is(val);
    auto D = [&](double& x) { return static_cast<bool>(is >> x); };
    auto I = [&](int& x) { return static_cast<bool>(is >> x); };
    auto L = [&](long& x) { return static_cast<bool>(is >> x); };
    auto B = [&](bool& x) { int v; if (!(is >> v)) return false; x = (v != 0); return true; };
    if (key == "model.ring") return D(model.ring);
    if (key == "model.ringProxy") {
      if (val == "jet") { model.ringProxy = kRingJet; return true; }
      if (val == "leadP") { model.ringProxy = kRingLeadP; return true; }
      return false;
    }
    if (key == "model.global") return D(model.global);
    if (key == "model.longitudinal") return D(model.longitudinal);
    if (key == "model.longHarmonic") return I(model.longHarmonic);
    if (key == "model.helicity") return D(model.helicity);
    if (key == "model.antiSign") return D(model.antiSign);
    if (key == "decay.ctauMM") return D(ctauMM);
    if (key == "decay.brPPi") return D(brPPi);
    if (key == "decay.brNPi0") return D(brNPi0);
    if (key == "decay.polarizeNPi0") return B(polarizeNPi0);
    if (key == "axes.trackPtMin") return D(trackPtMin);
    if (key == "axes.trackEtaMax") return D(trackEtaMax);
    if (key == "axes.pionMassHypothesis") return B(pionMassHypothesis);
    if (key == "axes.leadPtMin") return D(leadPtMin);
    if (key == "axes.jetR") return D(jetR);
    if (key == "axes.jetPtMin") return D(jetPtMin);
    if (key == "presel.minLambdas") return I(minLambdas);
    if (key == "presel.lambdaRapMax") return D(lambdaRapMax);
    if (key == "presel.lambdaPtMin") return D(lambdaPtMin);
    if (key == "presel.requireLeadP") return B(requireLeadP);
    if (key == "presel.requireJet") return B(requireJet);
    if (key == "presel.maxTries") return L(maxTries);
    if (key == "rp.random") return B(randomPsiRP);
    if (key == "verbose") return I(verbose);
    return false;
  }

  void print(std::ostream& os) const
  {
    os << "[InjectedPol] configuration:\n"
       << "  model: ring=" << model.ring << " (proxy " << (model.ringProxy == kRingJet ? "jet" : "leadP") << ")"
       << " global=" << model.global << " longitudinal=" << model.longitudinal << " (n=" << model.longHarmonic << ")"
       << " helicity=" << model.helicity << " antiSign=" << model.antiSign << "\n"
       << "  decay: ctauMM=" << ctauMM << " (<=0: Pythia tau0)  BR(p pi)=" << brPPi << " BR(n pi0)=" << brNPi0
       << " polarizeNPi0=" << polarizeNPi0 << "\n"
       << "  alpha: Lambda->p pi- " << constants::kAlphaLambdaToPPi << ", AntiLambda->pbar pi+ " << constants::kAlphaAntiLambdaToPPi
       << ", Lambda->n pi0 " << constants::kAlphaLambdaToNPi0 << ", AntiLambda->nbar pi0 " << constants::kAlphaAntiLambdaToNPi0 << "\n"
       << "  axes: track pT>" << trackPtMin << " |eta|<" << trackEtaMax << " pionMass=" << pionMassHypothesis
       << " leadP pT>" << leadPtMin << " anti-kT R=" << jetR << " jet pT>" << jetPtMin << " |eta_jet|<" << trackEtaMax - jetR
       << " (jet finding " << (needJets() ? "ON" : "OFF: not needed by this configuration") << ")\n"
       << "  presel: >=" << minLambdas << " (anti)Lambda with |y|<" << lambdaRapMax << " pT>" << lambdaPtMin
       << " requireLeadP=" << requireLeadP << " requireJet=" << requireJet << " maxTries=" << maxTries << "\n"
       << "  psiRP random=" << randomPsiRP << std::endl;
  }

  static std::string trim(const std::string& s)
  {
    const auto b = s.find_first_not_of(" \t\r");
    if (b == std::string::npos) return "";
    const auto e = s.find_last_not_of(" \t\r");
    return s.substr(b, e - b + 1);
  }
};

/// Extra per-event diagnostics (go into the MC event header in the O2 generator).
struct AxesInfo {
  int nCandidates = 0;
  double leadPt = 0., leadEta = 0., leadPhi = 0.;
  double jetPt = 0., jetEta = 0., jetPhi = 0.;
  int nJets = 0;
};

/// Build generator-level jet proxies from a Pythia event (final-state particles only).
/// `slowJet` (see makeSlowJet) and `work` are owned by the caller and reused event by event;
/// `work` must have been initialised once with work.init("injpol", &pythia.particleData).
inline EventAxes computeAxes(const Pythia8::Event& event, const Config& cfg,
                             Pythia8::SlowJet& slowJet, Pythia8::Event& work, AxesInfo& info)
{
  constexpr double kMassPi = 0.13957039;
  EventAxes ax;
  info = AxesInfo{};

  work.reset();
  int iLead = -1;
  double leadPt = 0.;
  for (int i = 0; i < event.size(); ++i) {
    const auto& p = event[i];
    if (!p.isFinal() || !p.isCharged()) continue;
    if (p.pT() < cfg.trackPtMin || std::abs(p.eta()) > cfg.trackEtaMax) continue;
    const double m = cfg.pionMassHypothesis ? kMassPi : p.m();
    const double e = std::sqrt(m * m + p.pAbs2());
    work.append(211, 91, 0, 0, p.px(), p.py(), p.pz(), e, m);
    ++info.nCandidates;
    if (p.pT() > leadPt) {
      leadPt = p.pT();
      iLead = i;
    }
  }
  if (iLead >= 0 && leadPt > cfg.leadPtMin) {
    const auto& p = event[iLead];
    ax.hasLeadP = true;
    ax.leadP = {p.px(), p.py(), p.pz()};
    info.leadPt = p.pT();
    info.leadEta = p.eta();
    info.leadPhi = p.phi();
  }

  if (cfg.needJets() && info.nCandidates > 0 && slowJet.analyze(work)) {
    // jets are ordered in decreasing pT; take the leading one inside the fiducial acceptance
    for (int j = 0; j < slowJet.sizeJet(); ++j) {
      const Pythia8::Vec4 pj = slowJet.p(j);
      if (std::abs(pj.eta()) > cfg.trackEtaMax - cfg.jetR) continue;
      if (pj.pT() < cfg.jetPtMin) break;
      ++info.nJets;
      if (!ax.hasJet) {
        ax.hasJet = true;
        ax.jet = {pj.px(), pj.py(), pj.pz()};
        info.jetPt = pj.pT();
        info.jetEta = pj.eta();
        info.jetPhi = pj.phi();
      }
    }
  }
  return ax;
}

/// SlowJet configured like the analysis: anti-kT (power -1), radius R, all particles of `work` used
/// (selection already done in computeAxes), true masses of `work` entries (= pion hypothesis if chosen).
/// The jet pT threshold is applied after the acceptance cut in computeAxes, so pTjetMin here is 0.
inline Pythia8::SlowJet* makeSlowJet(const Config& cfg)
{
  return new Pythia8::SlowJet(-1, cfg.jetR, 0., 1e3, /*select=*/1, /*massSet=*/2);
}

/// Count (anti)Lambdas left undecayed by Pythia inside the pre-selection acceptance.
inline int countLambdas(const Pythia8::Event& event, const Config& cfg)
{
  int n = 0;
  for (int i = 0; i < event.size(); ++i) {
    const auto& p = event[i];
    if (!p.isFinal() || std::abs(p.id()) != 3122) continue;
    if (std::abs(p.y()) > cfg.lambdaRapMax || p.pT() < cfg.lambdaPtMin) continue;
    ++n;
  }
  return n;
}

inline bool preselect(const Pythia8::Event& event, const Config& cfg, const EventAxes& ax)
{
  if (cfg.requireLeadP && !ax.hasLeadP) return false;
  if (cfg.requireJet && !ax.hasJet) return false;
  return countLambdas(event, cfg) >= cfg.minLambdas;
}

//==================================================================================================
// One (anti)Lambda decay: channel choice, polarization, kinematics and decay vertex.
// Shared by the O2 generator and the stand-alone harness so both consume random numbers identically:
//   flat (channel), flat (cos xi*), flat (phi*), exp (proper decay length).
//==================================================================================================
struct DecayConsts {
  double mP = 0., mN = 0., mPiCh = 0., mPi0 = 0.; // GeV/c^2
  double ctauMM = 78.9;                         // mm
  double fracPPi = 0.641;                       // P(p pi channel)
};

inline DecayConsts makeDecayConsts(const Config& cfg, Pythia8::ParticleData& pd)
{
  DecayConsts dc;
  dc.mP = pd.m0(2212);
  dc.mN = pd.m0(2112);
  dc.mPiCh = pd.m0(211);
  dc.mPi0 = pd.m0(111);
  dc.ctauMM = cfg.ctauMM > 0. ? cfg.ctauMM : pd.tau0(3122);
  const double brSum = cfg.brPPi + cfg.brNPi0;
  dc.fracPPi = brSum > 0. ? cfg.brPPi / brSum : -1.; // < 0 flags an invalid configuration
  return dc;
}

struct LambdaDecay {
  bool chargedChannel = true;
  int pdgB = 0, pdgM = 0; // baryon / meson daughters
  double alpha = 0.;      // decay parameter used
  V3 pol;                 // injected P* actually used (zero for an unpolarized n pi0 decay if so configured)
  bool clamped = false;   // |sum of models| > 1 was rescaled to 1
  TwoBodyDecay dec;       // lab 4-momenta
  V3 dx;                  // lab displacement of the decay vertex (mm)
  double dt = 0.;         // lab decay time (mm/c)
};

inline LambdaDecay decayLambda(const Config& cfg, const EventAxes& ax, const P4& pL, bool anti, const DecayConsts& dc, Pythia8::Rndm& rndm)
{
  LambdaDecay out;
  out.pol = polarization(cfg.model, ax, pL.p, anti, &out.clamped);
  out.chargedChannel = rndm.flat() < dc.fracPPi;
  double mB, mM;
  if (out.chargedChannel) {
    out.pdgB = anti ? -2212 : 2212;
    out.pdgM = anti ? 211 : -211;
    mB = dc.mP;
    mM = dc.mPiCh;
    out.alpha = anti ? constants::kAlphaAntiLambdaToPPi : constants::kAlphaLambdaToPPi;
  } else {
    out.pdgB = anti ? -2112 : 2112;
    out.pdgM = 111;
    mB = dc.mN;
    mM = dc.mPi0;
    out.alpha = anti ? constants::kAlphaAntiLambdaToNPi0 : constants::kAlphaLambdaToNPi0;
    if (!cfg.polarizeNPi0) {
      out.pol = V3{0., 0., 0.};
    }
  }
  const double u1 = rndm.flat();
  const double u2 = rndm.flat();
  out.dec = decayPolarized(pL, mB, mM, out.pol, out.alpha, u1, u2);
  decayDisplacement(pL, dc.ctauMM * rndm.exp(), out.dx, out.dt);
  return out;
}

} // namespace injpol

#endif // INJECTEDPOL_PYTHIA_H
