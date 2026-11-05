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
//        and p_B + p_M = p_Lambda. Tolerances scale with float precision (o2::MCTrack stores floats), and every
//        failure is printed with the offending quantity;
//   [P1] truth-level closure: R = 3/alpha <p_hat*_baryon . n_hat> (p pi channel), n_hat = (t_hat x p_Lambda)/|...|,
//        reproduces the injected ring amplitude; t_hat = generator TRIGGER direction (jet or leading particle)
//        stored in the MC header;
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

  Counter b1, b2, b3, b4, b5, b3e;
  long nEscapedTotal = 0;
  long nLambda = 0, nPPi = 0, nPrinted = 0, nPrintedB5 = 0;
  double sumR = 0., sumR2 = 0., sumL = 0., sumL2 = 0.;
  long nR = 0, nL = 0;
  long nDecayOutsideTPC = 0, nDecayBeyond10m = 0; // decay radius (3D distance from the PV) > 250 cm / > 1000 cm
  double maxDecayDist = 0.;
  double ringAmp = 0.;
  int trigger = kTriggerJet;
  // o2::MCTrack is MCTrackT<float>: positions/momenta carry ~7 significant digits. Tolerances must scale with the
  // magnitudes involved, otherwise Lambdas decaying metres away fail on rounding alone.
  constexpr double kFloatEps = 1.2e-7;

  for (int iev = 0; iev < nEv; ++iev) {
    const auto& tracks = reader.getTracks(iev);
    const auto& header = reader.getMCEventHeader(0, iev);

    bool hdrOk = true;
    ringAmp = headerDouble(header, "injpol.model.ring", hdrOk);
    {
      // "injpol.model.trigger" since 2026-11-05; older files have "injpol.model.ringProxy" (same meaning/values)
      bool okNew = true;
      trigger = header.hasInfo("injpol.model.trigger") ? headerInt(header, "injpol.model.trigger", okNew)
                                                       : headerInt(header, "injpol.model.ringProxy", hdrOk);
    }
    const bool useJet = (trigger == kTriggerJet);
    const bool hasTrigger = headerInt(header, useJet ? "injpol.hasJet" : "injpol.hasLeadP", hdrOk) == 1;
    const double tEta = headerDouble(header, useJet ? "injpol.jet.eta" : "injpol.leadP.eta", hdrOk);
    const double tPhi = headerDouble(header, useJet ? "injpol.jet.phi" : "injpol.leadP.phi", hdrOk);
    if (!hdrOk && iev == 0) {
      std::printf("WARNING: injpol.* header keys missing -- was the injected-polarization generator used?\n");
    }
    const V3 tHat = unit(V3{std::cos(tPhi), std::sin(tPhi), std::sinh(tEta)}); // trigger direction t_hat
    // Lambdas left undecayed because their decay vertex was outside the world ("escaped"; header key since 2026-11-05)
    const bool hasEscInfo = header.hasInfo("injpol.nEscaped");
    bool okEsc = true;
    const int nEscHeader = hasEscInfo ? headerInt(header, "injpol.nEscaped", okEsc) : 0;
    int nNoDaughters = 0;

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
      if (dau.empty() && ok1 && hasEscInfo) {
        // escaped Lambda: no decay by construction; B1/B2 still apply, B3-B5 do not
        ++nNoDaughters;
        b1.add(ok1);
        b2.add(ok2);
        continue;
      }
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
        const double dist = norm(disp);
        const V3 psum{D0.Px() + D1.Px(), D0.Py() + D1.Py(), D0.Pz() + D1.Pz()};

        // B5 sub-conditions with float-aware tolerances
        const double dVtx = norm(v0 - v1), tolVtx = 1e-6 + 4. * kFloatEps * std::max(norm(v0), norm(v1));
        const double perp = norm(cross(disp, unit(pLv))), tolPerp = 1e-6 + 8. * kFloatEps * (norm(vL) + norm(v0));
        const double along = dot(disp, unit(pLv));
        const double dP = norm(psum - pLv), tolP = 1e-6 + 8. * kFloatEps * (norm(pLv) + norm(V3{D0.Px(), D0.Py(), D0.Pz()}) + norm(V3{D1.Px(), D1.Py(), D1.Pz()}));
        const bool c1 = dVtx < tolVtx, c2 = perp < tolPerp, c3 = along >= -tolPerp, c4 = dP < tolP;
        ok5 = c1 && c2 && c3 && c4;
        if (!ok5 && nPrintedB5 < 20) {
          ++nPrintedB5;
          std::printf("  B5 detail: event %d track %zu pdg %d p=%.2f GeV/c decay distance=%.1f cm | "
                      "|v_d0-v_d1|=%.2e (tol %.1e)%s  perp=%.2e (tol %.1e)%s  along=%.2e%s  |sum p_d - p_L|=%.2e (tol %.1e)%s\n",
                      iev, i, L.GetPdgCode(), norm(pLv), dist, dVtx, tolVtx, c1 ? "" : " <-- FAIL", perp, tolPerp, c2 ? "" : " <-- FAIL",
                      along, c3 ? "" : " <-- FAIL", dP, tolP, c4 ? "" : " <-- FAIL");
        }

        // decay-distance bookkeeping
        maxDecayDist = std::max(maxDecayDist, dist);
        nDecayOutsideTPC += (dist > 250.);
        nDecayBeyond10m += (dist > 1000.);

        // proper decay length (cm -> mm)
        const double lProper = 10. * dist * L.GetMass() / norm(pLv);
        sumL += lProper;
        sumL2 += lProper * lProper;
        ++nL;

        // truth-level ring observable, p pi channel only: R = 3/alpha <p_hat*_B . n_hat>, n_hat = (t_hat x p_L)/|...|
        const MCTrack* B = nullptr;
        for (int d : dau) {
          if (std::abs(tracks[d].GetPdgCode()) == 2212) {
            B = &tracks[d];
          }
        }
        if (B && hasTrigger) {
          ++nPPi;
          P4 pL;
          pL.p = pLv;
          pL.e = L.GetEnergy();
          P4 pB;
          pB.p = {B->Px(), B->Py(), B->Pz()};
          pB.e = B->GetEnergy();
          const V3 pStarHat = unit(boostToRest(pB, pL, pL.mass()).p);
          const V3 n = cross(tHat, pLv);
          if (norm(n) > 1e-12) {
            const double alpha = anti ? constants::kAlphaAntiLambdaToPPi : constants::kAlphaLambdaToPPi;
            const double r = 3. / alpha * dot(pStarHat, n) / norm(n);
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
    if (hasEscInfo) {
      b3e.add(nNoDaughters == nEscHeader);
      nEscapedTotal += nNoDaughters;
    }
  }

  auto line = [](const char* name, const Counter& c) {
    if (c.pass == 0 && c.fail == 0) {
      return;
    }
    std::printf("[%s] %-72s pass %ld fail %ld\n", c.fail == 0 && c.pass > 0 ? " OK " : "FAIL", name, c.pass, c.fail);
  };
  std::printf("\ngenerator-level (anti)Lambdas: %ld (p pi channel with valid axis: %ld)\n", nLambda, nPPi);
  line("B1 Lambda not transported, HepMC status 1", b1);
  line("B2 Lambda is physical primary", b2);
  line("B3 exactly 2 daughters, all kPDecay (no Geant4 transport of the Lambda)", b3);
  line("B4 daughters: kPDecay, not primary, not physical primary, transported", b4);
  line("B5 common displaced vertex along p_Lambda, momentum conservation (float tol.)", b5);
  line("B3e per event: Lambdas without daughters == escaped count in the header (decay vertex outside world)", b3e);
  std::printf("[info] escaped Lambdas (not decayed, decay vertex outside the world): %ld\n", nEscapedTotal);
  if (nLambda > 0 && b3.pass == 0) {
    std::printf("[HINT] NO Lambda has its daughters in this file. If it comes from an O2DPG production, the generator most likely ran\n"
                "       in a SEPARATE step (o2dpg_sim_workflow.py default --event-gen-mode separated): the transport re-reads the events\n"
                "       with GeneratorFromO2Kine, which keeps only primary tracks and drops the daughters. Use --event-gen-mode integrated.\n");
  }
  if (nR > 1) {
    const double m = sumR / nR, e = std::sqrt(std::max(0., sumR2 / nR - m * m) / nR);
    std::printf("[%s] P1 truth-level ring R = %+.4f +- %.4f   (injected %+.4f, pull %+.2f)\n",
                std::abs(m - ringAmp) < 4. * e ? " OK " : "FAIL", m, e, ringAmp, (m - ringAmp) / e);
  }
  if (nL > 1) {
    const double m = sumL / nL, e = std::sqrt(std::max(0., sumL2 / nL - m * m) / nL);
    std::printf("[info] P2 <proper decay length> = %.2f +- %.2f mm (expected 78.9 mm with the default c*tau)\n", m, e);
    std::printf("[info] decay distance from the PV: > 250 cm (outside the TPC) %ld / %ld, > 10 m %ld, max %.0f cm\n",
                nDecayOutsideTPC, nL, nDecayBeyond10m, maxDecayDist);
  }
}
