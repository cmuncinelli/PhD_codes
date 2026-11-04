// checkAODFlags.C -- validation step 2 of InjectedPolO2MC (see ../DESIGN.md, section 7).
//
// Reads McParticles / McCollisions directly from an AO2D.root (no O2Physics needed) and checks that the
// injected-polarization Lambdas look like those of a standard production in the AO2D:
//   [A1] generator-level (anti)Lambda: producedByGenerator, PhysicalPrimary flag set;
//   [A2] it has 2 daughters in the McParticles slice, and both daughters point back to it as their mother;
//   [A3] daughters: ProducedByTransport flag set, getProcess() == 4 (kPDecay), PhysicalPrimary NOT set
//        (this is what PWGLF strangenessbuilder.cxx requires for generated-V0 MC info);
//   [A4] McCollisions.eventPlaneAngle is filled with the psi_RP used by the generator (spread over [0, 2pi)).
//
// Run with plain ROOT:
//   root -l -b -q 'checkAODFlags.C("tf1/AO2D.root")'
//
// Table/column names follow AliceO2 Framework/Core AnalysisDataModel.h (McParticles_001, McCollisions_001).
// Index columns are local to each DF_ directory.

#include "TDirectory.h"
#include "TFile.h"
#include "TKey.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"

#include <cmath>
#include <cstdio>
#include <vector>

void checkAODFlags(const char* fileName = "AO2D.root", int maxPrint = 5)
{
  constexpr unsigned char kProducedByTransport = 0x1, kPhysicalPrimary = 0x4;
  constexpr int kPDecay = 4;

  TFile f(fileName);
  if (f.IsZombie()) {
    std::printf("cannot open %s\n", fileName);
    return;
  }
  long nL = 0, a1 = 0, a2 = 0, a3 = 0, nPrinted = 0, nColl = 0, nPsiNonZero = 0;
  long nDauAll = 0, nDauPhysPrim = 0;
  double psiMin = 1e9, psiMax = -1e9;

  TIter next(f.GetListOfKeys());
  while (auto* key = static_cast<TKey*>(next())) {
    if (std::string(key->GetName()).rfind("DF_", 0) != 0) {
      continue;
    }
    auto* dir = static_cast<TDirectory*>(key->ReadObj());

    if (auto* tc = static_cast<TTree*>(dir->Get("O2mccollision_001"))) {
      TTreeReader rc(tc);
      TTreeReaderValue<float> psi(rc, "fEventPlaneAngle");
      while (rc.Next()) {
        ++nColl;
        nPsiNonZero += (*psi != 0.f);
        psiMin = std::min<double>(psiMin, *psi);
        psiMax = std::max<double>(psiMax, *psi);
      }
    }

    auto* t = static_cast<TTree*>(dir->Get("O2mcparticle_001"));
    if (!t) {
      std::printf("%s: no O2mcparticle_001 tree\n", key->GetName());
      continue;
    }
    // cache the needed columns of this DF
    std::vector<int> pdg, status;
    std::vector<unsigned char> flags;
    std::vector<int> d0, d1, m0; // daughters slice; first mother (-1 if none)
    {
      TTreeReader r(t);
      TTreeReaderValue<int> rPdg(r, "fPdgCode"), rStatus(r, "fStatusCode");
      TTreeReaderValue<unsigned char> rFlags(r, "fFlags");
      TTreeReaderArray<int> rDau(r, "fIndexSlice_Daughters");
      TTreeReaderArray<int> rMoth(r, "fIndexArray_Mothers");
      while (r.Next()) {
        pdg.push_back(*rPdg);
        status.push_back(*rStatus);
        flags.push_back(*rFlags);
        d0.push_back(rDau[0]);
        d1.push_back(rDau[1]);
        m0.push_back(rMoth.GetSize() > 0 ? rMoth[0] : -1);
      }
    }
    for (size_t i = 0; i < pdg.size(); ++i) {
      if (std::abs(pdg[i]) != 3122 || (flags[i] & kProducedByTransport)) {
        continue;
      }
      ++nL;
      const bool ok1 = (flags[i] & kPhysicalPrimary);
      const int nd = (d0[i] >= 0 && d1[i] >= d0[i]) ? d1[i] - d0[i] + 1 : 0;
      bool ok2 = (nd == 2);
      for (int d = d0[i]; ok2 && d <= d1[i]; ++d) {
        ok2 = ok2 && d < (int)pdg.size() && m0[d] == (int)i; // daughter -> mother back-reference
      }
      bool ok3 = ok2;
      for (int d = d0[i]; ok2 && d <= d1[i]; ++d) {
        ++nDauAll;
        nDauPhysPrim += (flags[d] & kPhysicalPrimary) ? 1 : 0;
        ok3 = ok3 && (flags[d] & kProducedByTransport) && status[d] == kPDecay && !(flags[d] & kPhysicalPrimary);
      }
      a1 += ok1;
      a2 += ok2;
      a3 += ok3;
      if (!(ok1 && ok2 && ok3) && nPrinted < maxPrint) {
        ++nPrinted;
        std::printf("  %s particle %zu pdg %d flags 0x%x: A1 %d A2 %d (ndau %d) A3 %d\n", key->GetName(), i, pdg[i], flags[i], ok1, ok2, nd, ok3);
      }
    }
  }

  auto line = [](const char* name, long pass, long tot) {
    std::printf("[%s] %-70s %ld / %ld\n", (tot > 0 && pass == tot) ? " OK " : "FAIL", name, pass, tot);
  };
  std::printf("\n%s: %ld generator-level (anti)Lambdas, %ld MC collisions\n", fileName, nL, nColl);
  line("A1 Lambda producedByGenerator && PhysicalPrimary", a1, nL);
  line("A2 Lambda has 2 daughters in McParticles, pointing back to it", a2, nL);
  line("A3 daughters ProducedByTransport, process 4, not PhysicalPrimary", a3, nL);
  std::printf("[info] daughters flagged PhysicalPrimary: %ld / %ld (must be 0)\n", nDauPhysPrim, nDauAll);
  if (nL > 0 && a3 < nL / 2) {
    std::printf("[HINT] most Lambdas have no proper decay daughters: was the production made with o2dpg_sim_workflow.py\n"
                "       --event-gen-mode integrated? (the default 'separated' mode drops the daughters; DESIGN.md section 4b)\n");
  }
  std::printf("[info] all particles carry FromBackgroundEvent (0x2) in non-embedded productions: the AOD producer flags source 0\n"
              "       as background, and 'sgn' is source 0 when there is no embedding. This is standard and harmless.\n");
  std::printf("[%s] A4 eventPlaneAngle non-zero in %ld / %ld collisions, range [%.3f, %.3f] (expect ~[0, 6.283] if rp.random = 1)\n",
              nPsiNonZero > 0 ? " OK " : "WARN", nPsiNonZero, nColl, psiMin, psiMax);
}
