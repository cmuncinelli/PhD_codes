// generator_pythia8_injectedPol.C
//
// O2 external generator: Pythia8 + polarized (anti)Lambda decay afterburner, for full-detector closure tests
// of Lambda polarization measurements. TEST CODE -- see ../README.md and ../DESIGN.md. NOT an official
// ALICE generator; never use its output as, or mix it with, an official MC production.
//
// Usage (via O2DPG, see ../run/runInjectedPol.sh):
//   (the run/ scripts set INJPOL_ROOT = <REPO_DIR>/InjectedPolO2MC themselves, from their hard-coded REPO_DIR;
//    if you ever set it by hand: absolute path, NO trailing "/")
//   o2dpg_sim_workflow.py ... -gen external -ini ${INJPOL_ROOT}/generator/injectedPol.ini
// with injectedPol.ini containing
//   [GeneratorExternal]
//   fileName = ${INJPOL_ROOT}/generator/generator_pythia8_injectedPol.C
//   funcName = generator_injectedPol("${INJPOL_ROOT}/generator/configs/ring_0p30.cfg")
//   [GeneratorPythia8]
//   config = ${INJPOL_ROOT}/generator/pythia8_pp_13600.cfg
//
// What it does, per event:
//   1. Pythia8 generates pp events (Lambdas are NOT decayed by Pythia) until the generator-level
//      pre-selection passes (>= N Lambdas in acceptance, leading particle / jet if requested).
//   2. Generator-level jet proxies (leading particle, anti-kT jet) mirroring the reconstruction-level ones.
//   3. Every to-be-transported (anti)Lambda gets an injected polarization P* (models in InjectedPolKernel.h),
//      is decayed into p pi- / n pi0 (and charge conjugates) with dN/dcos = (1 + alpha |P*| cos)/2, at a
//      decay vertex sampled from its lifetime. The daughters are appended to the event and handed to Geant4.
//
// MC bookkeeping (see ../DESIGN.md, section 4) is chosen to look like a standard production where
// Geant4 decays the Lambda:
//   - Lambda : HepMC status kept at 1 (=> stays a physical primary), NOT transported, daughters linked;
//   - p, pi  : TMCProcess kPDecay as creation process (=> not "primary", ProducedByTransport in the AO2D,
//              getProcess() == 4 as required by the strangeness builder), displaced production vertex.
//
// Inspired by / adapted from:
//   - AliceO2   Generators/src/EventDecayer.cxx            (EventDecayer::decayOne: bookkeeping of generator-level decays)
//   - O2DPG     MC/config/PWGEM/external/generator/Generator_pythia8_forcedDecays.C (importParticles() override pattern)
//   - PhD_codes vorticity_hydro/lambda_pol_toy_model.cxx   (polarized decay kernel)

R__ADD_INCLUDE_PATH($INJPOL_ROOT/generator)

#include "Generators/GeneratorPythia8.h"
#include "SimulationDataFormat/MCEventHeader.h"
#include "SimulationDataFormat/MCGenProperties.h"
#include "SimulationDataFormat/MCUtils.h"
#include "SimulationDataFormat/ParticleStatus.h"
#include "Pythia8/Pythia.h"
#include "TMCProcess.h"
#include "TParticle.h"
#include "TString.h"
#include "TSystem.h"
#include <fairlogger/Logger.h>

#include "InjectedPolKernel.h"
#include "InjectedPolPythia.h"

#include <cmath>
#include <cstdlib>
#include <memory>
#include <sstream>
#include <string>

namespace o2
{
namespace eventgen
{

class GeneratorPythia8InjectedPol : public GeneratorPythia8
{
 public:
  explicit GeneratorPythia8InjectedPol(const char* cfgPath)
  {
    if (cfgPath && *cfgPath) {
      TString path(cfgPath);
      gSystem->ExpandPathName(path);
      if (!mCfg.readFile(path.Data())) {
        LOG(fatal) << "[InjectedPol] invalid configuration file " << path.Data();
      }
      mCfgPath = path.Data();
    } else {
      LOG(warning) << "[InjectedPol] no configuration file given: running with defaults (NO polarization)";
    }
  }

  ~GeneratorPythia8InjectedPol() override { printSummary(); }

  Bool_t Init() override
  {
    // Lambdas must reach importParticles() undecayed. They already are with the standard
    // 'ParticleDecays:limitTau0 = on / tau0Max = 10' (c tau = 78.9 mm), but we make it explicit.
    mPythia.readString("3122:mayDecay = off"); // also applies to the antiparticle
    if (!GeneratorPythia8::Init()) {
      return kFALSE;
    }
    if (mPythia.particleData.mayDecay(3122)) {
      LOG(fatal) << "[InjectedPol] Pythia is configured to decay Lambdas (check the Pythia card): refusing to run";
    }

    auto& pd = mPythia.particleData;
    mDecayConsts = injpol::makeDecayConsts(mCfg, pd);
    if (!(mDecayConsts.fracPPi > 0.)) {
      LOG(fatal) << "[InjectedPol] invalid branching ratios";
    }

    mWork.init("InjectedPol work event", &mPythia.particleData);
    mSlowJet.reset(injpol::makeSlowJet(mCfg));

    std::ostringstream os;
    mCfg.print(os);
    LOG(info) << os.str();
    LOG(info) << "[InjectedPol] config file: " << (mCfgPath.empty() ? "<none>" : mCfgPath);
    LOG(info) << "[InjectedPol] masses (GeV): p " << mDecayConsts.mP << " n " << mDecayConsts.mN << " pi+- " << mDecayConsts.mPiCh
              << " pi0 " << mDecayConsts.mPi0 << " | Lambda c*tau used = " << mDecayConsts.ctauMM << " mm (Pythia tau0 = " << pd.tau0(3122)
              << " mm) | fraction p pi = " << mDecayConsts.fracPPi;
    return kTRUE;
  }

  /// Pythia event + generator-level pre-selection (before any Geant4 cost).
  Bool_t generateEvent() override
  {
    long tries = 0;
    while (true) {
      if (++tries > mCfg.maxTries) {
        LOG(fatal) << "[InjectedPol] pre-selection not satisfied after " << mCfg.maxTries << " Pythia events";
        return kFALSE;
      }
      if (!GeneratorPythia8::generateEvent()) {
        continue; // Pythia failure: just try again (counted in tries)
      }
      mAxes = injpol::computeAxes(mPythia.event, mCfg, *mSlowJet, mWork, mAxesInfo);
      if (injpol::preselect(mPythia.event, mCfg, mAxes)) {
        break;
      }
    }
    mAxes.psiRP = mCfg.randomPsiRP ? 2. * injpol::constants::kPi * mPythia.rndm.flat() : 0.;
    mEventTries = tries;
    mTotalTries += tries;
    ++mNAccepted;
    return kTRUE;
  }

  Bool_t importParticles() override
  {
    if (!GeneratorPythia8::importParticles()) {
      return kFALSE;
    }
    injectPolarizedDecays();
    return kTRUE;
  }

  void updateHeader(o2::dataformats::MCEventHeader* header) override
  {
    GeneratorPythia8::updateHeader(header);
    using Key = o2::dataformats::MCInfoKeys;
    // Reaction-plane angle used by the global/longitudinal models -> AO2D McCollisions.eventPlaneAngle
    header->SetRotZ(mAxes.psiRP);
    header->putInfo<double>(Key::planeAngle, mAxes.psiRP);
    // Everything else stays in o2sim_MCHeader.root (not propagated to the AO2D)
    header->putInfo<std::string>("injpol.config", mCfgPath);
    header->putInfo<double>("injpol.model.ring", mCfg.model.ring);
    header->putInfo<int>("injpol.model.trigger", mCfg.model.trigger); // 0 = leadP, 1 = jet (t_hat of the ring model)
    header->putInfo<double>("injpol.model.global", mCfg.model.global);
    header->putInfo<double>("injpol.model.longitudinal", mCfg.model.longitudinal);
    header->putInfo<int>("injpol.model.longHarmonic", mCfg.model.longHarmonic);
    header->putInfo<double>("injpol.model.helicity", mCfg.model.helicity);
    header->putInfo<double>("injpol.model.antiSign", mCfg.model.antiSign);
    header->putInfo<int>("injpol.nTries", static_cast<int>(mEventTries));
    header->putInfo<int>("injpol.nDecayed", mEventDecayed);
    header->putInfo<int>("injpol.hasLeadP", mAxes.hasLeadP ? 1 : 0);
    header->putInfo<double>("injpol.leadP.pt", mAxesInfo.leadPt);
    header->putInfo<double>("injpol.leadP.eta", mAxesInfo.leadEta);
    header->putInfo<double>("injpol.leadP.phi", mAxesInfo.leadPhi);
    header->putInfo<int>("injpol.hasJet", mAxes.hasJet ? 1 : 0);
    header->putInfo<double>("injpol.jet.pt", mAxesInfo.jetPt);
    header->putInfo<double>("injpol.jet.eta", mAxesInfo.jetEta);
    header->putInfo<double>("injpol.jet.phi", mAxesInfo.jetPhi);
    header->putInfo<int>("injpol.nJets", mAxesInfo.nJets);
  }

 private:
  void injectPolarizedDecays()
  {
    mEventDecayed = 0;
    const int nOrig = mParticles.size(); // daughters are appended after this index
    for (int i = 0; i < nOrig; ++i) {
      TParticle& lambda = mParticles[i];
      const int pdg = lambda.GetPdgCode();
      if (std::abs(pdg) != 3122 || !lambda.TestBit(ParticleStatus::kToBeDone)) {
        continue;
      }
      const bool anti = (pdg < 0);

      injpol::P4 pL;
      pL.p = {lambda.Px(), lambda.Py(), lambda.Pz()};
      pL.e = lambda.Energy();
      const injpol::LambdaDecay d = injpol::decayLambda(mCfg, mAxes, pL, anti, mDecayConsts, mPythia.rndm);
      if (d.clamped) {
        ++mNClamped;
      }
      // decay vertex, in the generator units of mParticles (Pythia: mm and mm/c)
      const double vx = lambda.Vx() + d.dx.x, vy = lambda.Vy() + d.dx.y, vz = lambda.Vz() + d.dx.z, vt = lambda.T() + d.dt;
      const double weight = lambda.GetWeight();

      // Mother: keep HepMC status 1 (=> physical primary, as in standard productions), do not transport.
      const int iDau = mParticles.size();
      o2::mcutils::MCGenHelper::encodeParticleStatusAndTracking(lambda, 1, o2::mcgenstatus::getGenStatusCode(lambda.GetStatusCode()), false);
      lambda.SetFirstDaughter(iDau);
      lambda.SetLastDaughter(iDau + 1);
      // NOTE: `lambda` must not be used after the push_backs below (the vector may reallocate)

      pushDaughter(d.pdgB, d.dec.baryon, i, vx, vy, vz, vt, weight);
      pushDaughter(d.pdgM, d.dec.meson, i, vx, vy, vz, vt, weight);
      ++mEventDecayed;
      ++mNDecayed;
      if (d.chargedChannel) {
        ++mNDecayedPPi;
      }
    }
  }

  void pushDaughter(int pdg, const injpol::P4& p, int iMother, double vx, double vy, double vz, double vt, double weight)
  {
    TParticle d(pdg, 0, iMother, -1, -1, -1, p.p.x, p.p.y, p.p.z, p.e, vx, vy, vz, vt);
    o2::mcutils::MCGenHelper::encodeParticleStatusAndTracking(d, 1, 91, true); // 91: Pythia code for decay products
    d.SetUniqueID(kPDecay);                                                     // creation process -> secondary-like bookkeeping
    d.SetWeight(weight);
    mParticles.push_back(d);
  }

  void printSummary() const
  {
    if (mNAccepted == 0) {
      return;
    }
    LOG(info) << "[InjectedPol] summary: accepted events " << mNAccepted << " / Pythia events " << mTotalTries
              << " (efficiency " << double(mNAccepted) / mTotalTries << ")"
              << " | decayed (anti)Lambdas " << mNDecayed << " (p pi: " << mNDecayedPPi << ")"
              << " | |P|>1 clamped " << mNClamped;
  }

  injpol::Config mCfg;
  std::string mCfgPath;
  std::unique_ptr<Pythia8::SlowJet> mSlowJet;
  Pythia8::Event mWork;
  injpol::EventAxes mAxes;
  injpol::AxesInfo mAxesInfo;

  injpol::DecayConsts mDecayConsts;

  long mEventTries = 0, mTotalTries = 0, mNAccepted = 0;
  int mEventDecayed = 0;
  long mNDecayed = 0, mNDecayedPPi = 0, mNClamped = 0;
};

} // namespace eventgen
} // namespace o2

/// Entry point used by [GeneratorExternal] funcName. `cfgPath` may contain environment variables.
FairGenerator* generator_injectedPol(const char* cfgPath = "")
{
  return new o2::eventgen::GeneratorPythia8InjectedPol(cfgPath);
}
