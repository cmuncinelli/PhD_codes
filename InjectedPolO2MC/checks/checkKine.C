// checkKine.C -- validation step 1 of InjectedPolO2MC (see ../DESIGN.md, section 7).
//
// Inspects the o2-sim kinematics (<prefix>_Kine.root + <prefix>_MCHeader.root) produced with the
// injected-polarization generator and checks, for every generator-level (anti)Lambda:
//   [B1] it was NOT transported (toBeDone == false) and keeps HepMC status 1;
//   [B2] it is a physical primary (o2::mcutils::MCTrackNavigator::isPhysicalPrimary);
//   [B3] it has exactly 2 daughters (tracks whose mother is the Lambda), none created by Geant4 processes
//        other than kPDecay (anything else would mean Geant4 also transported the Lambda);
//   [B4] daughters: process kPDecay, isPrimary() == false, NOT physical primaries, transported;
//   [B5] daughters share one production vertex, displaced from the Lambda vertex ALONG the Lambda momentum,
//        and p_B + p_M = p_Lambda;
//   [P1] truth-level closure: R = 3/alpha <p_hat*_baryon . t_hat> (p pi channel) reproduces the injected
//        ring amplitude, with t_hat built from the generator jet/leadP axis stored in the MC header;
//   [P2] <proper decay length> ~ 78.9 mm.
//
// Run: through run/smokeTestSim.sh or run/runInjectedPol.sh, which call it COMPILED (ACLiC), as is the O2 convention
// for macros using O2 headers. Interpreting it makes ROOT auto-parse the VMC dictionary (TMCProcess) recursively and
// fail. By hand, from the directory with <prefix>_Kine.root, inside the O2 environment:
//   root -l -b -q -e 'gSystem->SetBuildDir("/tmp/aclic", kTRUE); gSystem->Load("libO2CommonUtils"); gSystem->Load("libO2SimulationDataFormat"); gSystem->Load("libO2Steer");' \
//        '/home/users/cicerodm/PhD_codes/InjectedPolO2MC/checks/checkKine.C+("o2sim")'
// where "o2sim" is the o2-sim output prefix (O2DPG productions use "sgn" inside each tfN/ directory).

#if !defined(__CLING__) || defined(__ROOTCLING__)
#include "../generator/InjectedPolKernel.h" // relative to this file: no environment variable needed
#include "SimulationDataFormat/MCEventHeader.h"
#include "SimulationDataFormat/MCGenProperties.h"
#include "SimulationDataFormat/MCTrack.h"
#include "SimulationDataFormat/MCUtils.h"
#include "Steer/MCKinematicsReader.h"
#include "TMCProcess.h"

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>
#endif

namespace
{
struct Counter {
  long pass = 0, fail = 0;
  void add(bool ok) { ok ? ++pass : ++fail; }
};

double headerDouble(const o2::dataformats::MCEventHeader& h, const char* key, bool& ok)
{
  bool valid = false;
  const double v = h.getInfo<double>(key, valid);
  ok = ok && valid;
  return valid ? v : 0.;
}
int headerInt(const o2::dataformats::MCEventHeader& h, const char* key, bool& ok)
{
  bool valid = false;
  const int v = h.getInfo<int>(key, valid);
  ok = ok && valid;
  return valid ? v : 0;
}
} // namespace

void checkKine(const char* prefix = "o2sim", int maxPrint = 5)
{
  using o2::MCTrack;
  using o2::mcutils::MCTrackNavigator;
  using namespace injpol;

  o2::steer::MCKinematicsReader reader(prefix, o2::steer::MCKinematicsReader::Mode::kMCKine);
  const int nEv = reader.getNEvents(0);
  std::printf("checkKine: %d events in %s_Kine.root\n", nEv, prefix);

  Counter b1, b2, b3, b4, b5;
  long nLambda = 0, nPPi = 0, nPrinted = 0;
  double sumR = 0., sumR2 = 0., sumL = 0., sumL2 = 0.;
  long nR = 0, nL = 0;
  double ringAmp = 0.;
  int ringProxy = kRingJet;

  for (int iev = 0; iev < nEv; ++iev) {
    const auto& tracks = reader.getTracks(iev);
    const auto& header = reader.getMCEventHeader(0, iev);

    bool hdrOk = true;
    ringAmp = headerDouble(header, "injpol.model.ring", hdrOk);
    ringProxy = headerInt(header, "injpol.model.ringProxy", hdrOk);
    const bool useJet = (ringProxy == kRingJet);
    const bool hasAxis = headerInt(header, useJet ? "injpol.hasJet" : "injpol.hasLeadP", hdrOk) == 1;
    const double aEta = headerDouble(header, useJet ? "injpol.jet.eta" : "injpol.leadP.eta", hdrOk);
    const double aPhi = headerDouble(header, useJet ? "injpol.jet.phi" : "injpol.leadP.phi", hdrOk);
    if (!hdrOk && iev == 0) {
      std::printf("WARNING: injpol.* header keys missing -- was the injected-polarization generator used?\n");
    }
    const V3 axis{std::cos(aPhi), std::sin(aPhi), std::sinh(aEta)}; // direction only

    // mother -> daughters map built from the mother indices (robust against daughter-index conventions)
    std::vector<std::vector<int>> daughters(tracks.size());
    for (size_t i = 0; i < tracks.size(); ++i) {
      const int m = tracks[i].getMotherTrackId();
      if (m >= 0 && m < (int)tracks.size()) {
        daughters[m].push_back(i);
      }
    }

    for (size_t i = 0; i < tracks.size(); ++i) {
      const MCTrack& L = tracks[i];
      if (std::abs(L.GetPdgCode()) != 3122 || !L.isPrimary()) {
        continue; // only generator-level Lambdas (Geant4 Lambdas from Xi/Omega decays are left alone)
      }
      ++nLambda;
      const bool anti = L.GetPdgCode() < 0;

      const bool ok1 = !L.getToBeDone() && o2::mcgenstatus::getHepMCStatusCode(L.getStatusCode()) == 1;
      const bool ok2 = MCTrackNavigator::isPhysicalPrimary(L, tracks);
      const auto& dau = daughters[i];
      bool ok3 = (dau.size() == 2);
      bool ok4 = ok3, ok5 = ok3;
      for (int d : dau) {
        const MCTrack& D = tracks[d];
        ok3 = ok3 && (D.getProcess() == kPDecay);
        ok4 = ok4 && (D.getProcess() == kPDecay) && !D.isPrimary() && !MCTrackNavigator::isPhysicalPrimary(D, tracks) && D.getToBeDone();
      }
      if (ok3) {
        const MCTrack& D0 = tracks[dau[0]];
        const MCTrack& D1 = tracks[dau[1]];
        const V3 v0{D0.Vx(), D0.Vy(), D0.Vz()}, v1{D1.Vx(), D1.Vy(), D1.Vz()}, vL{L.Vx(), L.Vy(), L.Vz()};
        const V3 pLv{L.Px(), L.Py(), L.Pz()};
        const V3 disp = v0 - vL;
        const double perp = norm(cross(disp, unit(pLv)));
        const V3 psum{D0.Px() + D1.Px(), D0.Py() + D1.Py(), D0.Pz() + D1.Pz()};
        ok5 = norm(v0 - v1) < 1e-6 && perp < 1e-4 && dot(disp, pLv) >= 0. && norm(psum - pLv) < 1e-5;

        // proper decay length (cm -> mm)
        const double lProper = 10. * norm(disp) * L.GetMass() / norm(pLv);
        sumL += lProper;
        sumL2 += lProper * lProper;
        ++nL;

        // truth-level ring observable, p pi channel only
        const MCTrack* B = nullptr;
        for (int d : dau) {
          if (std::abs(tracks[d].GetPdgCode()) == 2212) {
            B = &tracks[d];
          }
        }
        if (B && hasAxis) {
          ++nPPi;
          P4 pL;
          pL.p = pLv;
          pL.e = L.GetEnergy();
          P4 pB;
          pB.p = {B->Px(), B->Py(), B->Pz()};
          pB.e = B->GetEnergy();
          const V3 pStarHat = unit(boostToRest(pB, pL, pL.mass()).p);
          const V3 cr = cross(unit(axis), pLv);
          if (norm(cr) > 1e-12) {
            const double alpha = anti ? constants::kAlphaAntiLambdaToPPi : constants::kAlphaLambdaToPPi;
            const double r = 3. / alpha * dot(pStarHat, cr) / norm(cr);
            sumR += r;
            sumR2 += r * r;
            ++nR;
          }
        }
      }
      b1.add(ok1);
      b2.add(ok2);
      b3.add(ok3);
      b4.add(ok4);
      b5.add(ok5);
      if (!(ok1 && ok2 && ok3 && ok4 && ok5) && nPrinted < maxPrint) {
        ++nPrinted;
        std::printf("  event %d track %zu pdg %d: B1 %d B2 %d B3 %d (ndau %zu) B4 %d B5 %d\n", iev, i, L.GetPdgCode(), ok1, ok2, ok3, dau.size(), ok4, ok5);
      }
    }
  }

  auto line = [](const char* name, const Counter& c) {
    std::printf("[%s] %-72s pass %ld fail %ld\n", c.fail == 0 && c.pass > 0 ? " OK " : "FAIL", name, c.pass, c.fail);
  };
  std::printf("\ngenerator-level (anti)Lambdas: %ld (p pi channel with valid axis: %ld)\n", nLambda, nPPi);
  line("B1 Lambda not transported, HepMC status 1", b1);
  line("B2 Lambda is physical primary", b2);
  line("B3 exactly 2 daughters, all kPDecay (no Geant4 transport of the Lambda)", b3);
  line("B4 daughters: kPDecay, not primary, not physical primary, transported", b4);
  line("B5 common displaced vertex along p_Lambda, momentum conservation", b5);
  if (nR > 1) {
    const double m = sumR / nR, e = std::sqrt(std::max(0., sumR2 / nR - m * m) / nR);
    std::printf("[%s] P1 truth-level ring R = %+.4f +- %.4f   (injected %+.4f, pull %+.2f)\n",
                std::abs(m - ringAmp) < 4. * e ? " OK " : "FAIL", m, e, ringAmp, (m - ringAmp) / e);
  }
  if (nL > 1) {
    const double m = sumL / nL, e = std::sqrt(std::max(0., sumL2 / nL - m * m) / nL);
    std::printf("[info] P2 <proper decay length> = %.2f +- %.2f mm (expected 78.9 mm with the default c*tau)\n", m, e);
  }
}
