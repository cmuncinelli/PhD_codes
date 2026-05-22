// ==========================================================================
// helicityEfficiencyToyModel.cxx
// ==========================================================================
//
// PURPOSE
// -------
// This standalone ROOT macro is a Monte Carlo toy model designed to isolate
// and diagnose the two main detector-induced fake-polarization effects that
// appear in Lambda (and AntiLambda) hyperon analyses:
//
//   (1) THE HELICITY EFFICIENCY EFFECT (forward-backward asymmetry):
//       A minimum-pT cut on daughter tracks preferentially removes the
//       decay configuration where the pion was emitted BACKWARD in the
//       Lambda rest frame. In that configuration the pion is soft in the lab;
//       the proton well-reconstructed in both geometries, but the soft pion
//       is not.  As a result, the isotropic decay angular distribution develops
//       a fake asymmetry in cos(theta*), where theta* is the angle between the proton
//       emission direction in the Lambda rest frame and the Lambda lab
//       momentum direction.  This is the "negative helicity" effect.
//
//   (2) THE DCA-INDUCED LEFT-RIGHT ASYMMETRY (phi* asymmetry):
//       A minimum-DCA cut requires each daughter track to be displaced from
//       the primary vertex (PV).  Because the proton and the pion carry
//       opposite charges, they curve in opposite directions in the magnetic
//       field.  Depending on the azimuthal angle phi* of the proton emission
//       direction around the Lambda axis, one of the two daughters may curve
//       TOWARD the PV (failing the DCA cut) while the other curves away.
//       This geometrical preference creates a fake asymmetry in phi*, the
//       "left-right" direction. Crucially, the sign of this asymmetry should
//       flip between Lambda pseudorapidity eta > 0 and eta < 0, because
//       flipping pz reverses the sense of the in-out coordinate (e3), which
//       mirrors the phi* distribution, and the effect should be exactly the
//       same for both eta > 0 and eta < 0, so just the reference direction flips.
//
// WHAT THIS TOY CAN DO THAT A FULL GEANT4 SIMULATION CANNOT:
// -----------------------------------------------------------
//   (a) Isolate each effect independently by switching cuts on/off.
//   (b) Use an ideal sharp pT threshold and a purely analytical DCA formula,
//       with no smearing from multiple scattering or detector material.
//   (c) Flip the magnetic field sign or the eta sign in a single line to
//       verify the sign of the asymmetry.
//   (d) Generate O(10^7) Lambdas in minutes, rather than hours.
//   (e) Guarantee that no polarization is input, so any non-zero result is
//       a pure acceptance/reconstruction artefact.
// What it CANNOT do: quantitative correction factors (which require GEANT4).
//
// COORDINATE SYSTEM  (defined per Lambda)
// ----------------------------------------
//   e1 = p_Lambda_unit                      (forward-backward axis)
//   e2 = (p_Lambda x z_hat) / |...|         (left-right axis; phi* = 0 here)
//   e3 = e1 x e2                            (in-out axis)
//
//   cos(theta*) = p*_proton_unit . e1
//   phi*        = atan2(p*_proton_unit . e3,
//                       p*_proton_unit . e2)
//
// Note: e2 = (py, -px, 0) / pT  and  e3 = (pz*px, pz*py, -pT^2) / (|p|*pT)
// When the Lambda is nearly collinear with z (pT -> 0), these axes become
// ill-defined; such events are skipped (see guard in main loop).
//
// RING OBSERVABLE PROXY
// ---------------------
// The ring observable uses the beam direction (+z) as the jet proxy:
//
//   R_proxy = (3 / alpha_Lambda) *
//             p*_proton_unit . (z_hat x lambda_unit)
//                            / |z_hat x lambda_unit|
//
// Algebraically this reduces to:
//   R_proxy = (3/alpha) * p*_proton_unit . (-py, px, 0) / pT
//           = -(3/alpha) * p*_proton_unit . e2
//           = -(3/alpha) * sin(theta*) * cos(phi*)
//
// Therefore R_proxy directly measures the LEFT-RIGHT (phi*) asymmetry.
//   - pT cut alone:   R_proxy stays near zero  (helicity effect is along e1,
//                     not e2, so it does not directly bias R_proxy).
//   - DCA cut alone:  R_proxy != 0, and it FLIPS SIGN between eta>0 and eta<0.
//   - Both cuts:      Combination of the above.
//
// EXPECTED RESULTS
// ----------------
//   No cuts  (bug-check):  (cos theta*, phi*) uniform.  R_proxy ~ 0.
//   pT cut only:           cos(theta*) suppressed at positive values.
//                          phi* still flat.  R_proxy ~ 0.
//   DCA cut only:          phi* shows a sin-like modulation.
//                          R_proxy != 0.  Sign FLIPS eta>0 vs eta<0.
//   Both cuts:             Both asymmetries superimposed.
//
// USAGE
// -----
// From the ROOT command line:
//
//   root -l -b -q 'helicityEfficiencyToyModel.cxx'
//       (all defaults: 1M Lambdas, output to helicityEffOutput.root)
//
//   root -l -b -q 'helicityEfficiencyToyModel.cxx(5000000,"myout.root")'
//       (5M Lambdas, custom output file)
//
//   root -l -b -q 'helicityEfficiencyToyModel.cxx(nL,"out.root",Bz,pTmin,pTmax,ymax,T,pTminP,pTminPi,dcaP,dcaPi,seed)'
//       (full parameter control -- see function signature below)
//
// OUTPUT STRUCTURE (in the ROOT file)
// ------------------------------------
// Two parallel families of cut scenarios are stored:
//
//   WithoutEtaGate/
//     NoCuts/      EtaPos/  EtaNeg/  All/
//     pTCutOnly/   EtaPos/  EtaNeg/  All/
//     DCACutOnly/  EtaPos/  EtaNeg/  All/
//     BothCuts/    EtaPos/  EtaNeg/  All/
//
//   WithEtaGate/   (additionally requires |eta_daughter| < etaMaxDetector)
//     NoCuts/      EtaPos/  EtaNeg/  All/
//     pTCutOnly/   EtaPos/  EtaNeg/  All/
//     DCACutOnly/  EtaPos/  EtaNeg/  All/
//     BothCuts/    EtaPos/  EtaNeg/  All/
//
//   Kinematics/
//     (Lambda pT spectrum, rapidity, decay radius, daughter pT before any cut)
//
// WithoutEtaGate is kept to demonstrate why omitting the daughter eta
// acceptance gate makes the physics cuts inconsistent.  WithEtaGate is
// the physically correct set and should be used for all conclusions.
//
// EtaPos / EtaNeg correspond to eta_Lambda >= 0 / < 0.
// All combines both halves.
//
// Each leaf directory (EtaPos, EtaNeg, All) contains:
//   h2d_cosTheta_phi  -- 2D: cos(theta*) vs phi*  [main diagnostic]
//   h1d_cosTheta      -- 1D projection of cos(theta*)
//   h1d_phi           -- 1D projection of phi*
//   h1d_ringProxy     -- distribution of R_proxy values
//   pRingProxy        -- TProfile: integrated <R_proxy>
//   pRingProxyVsEta   -- TProfile: <R_proxy> vs Lambda eta
//   pRingProxyVsPt    -- TProfile: <R_proxy> vs Lambda pT
//   h1d_ringProxyJet  -- Same as above, but using a random direction as if the jet proxy
//   pRingProxyJet     -- Same as above, but using a random direction as if the jet proxy
//   pRingProxyJetVsEta -- Same as above, but using a random direction as if the jet proxy
//   pRingProxyJetVsEtaJet -- TProfile: Sees the dependency with the random jet proxy's eta
//   pRingProxyJetVsPt -- Same as above, but using a random direction as if the jet proxy
//   h1d_decayRadius   -- transverse decay vertex radius [cm]
//   h1d_pT_proton     -- proton pT in lab frame
//   h1d_pT_pion       -- pion pT in lab frame
//   h1d_DCA_proton    -- proton DCA_xy to PV
//   h1d_DCA_pion      -- pion DCA_xy to PV
//   h1d_pT_lambda     -- Lambda pT
//   h1d_eta_lambda    -- Lambda eta
//
//   Reference analysis code: lambdaJetPolarizationIonsDerived.cxx (O2Physics)
// ==========================================================================

// --- Standard ROOT headers ---
// (these are only needed when code is compiled via g++, actually)
#include <TFile.h>
#include <TDirectory.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TLorentzVector.h> // Preserved solely for TGenPhaseSpace
#include <Math/Vector4D.h>
#include <Math/GenVector/VectorUtil.h>
// #include <TVector3.h>
#include <Math/Vector3D.h>
#include <TGenPhaseSpace.h>
#include <TRandom3.h>
#include <TMath.h>
#include <TString.h>
#include <TROOT.h>

// --- Standard C++ headers ---
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <iostream>
#include <exception>

using ROOT::Math::PxPyPzEVector;
using ROOT::Math::XYZVector; // From Vector3D

// ==========================================================================
// PHYSICS CONSTANTS
// All energies and masses in GeV, distances in cm, momenta in GeV/c.
// ==========================================================================

// Particle masses [GeV/c^2]
static const double kMassLambda = 1.115683;
static const double kMassProton = 0.938272;
static const double kMassPion   = 0.139570;

// Lambda mean proper decay length c*tau [cm]
static const double kCTauLambda = 7.89;

// Weak-decay asymmetry parameter (DPG 2025 update)
// alpha_Lambda = +0.749  for Lambda -> p + pi-
// alpha_Lambda = -0.758  for AntiLambda -> pbar + pi+
// Here we only simulate Lambda decays, so we use the positive value.
static const double kAlphaLambda  = 0.749;
static const double kPolPrefactor = 3.0 / kAlphaLambda; // (3/alpha) normalization

// Magnetic-field conversion: radius [cm] = pT [GeV/c] / (kBConv * B [T])
// This follows from R [m] = pT / (0.3 * B), then R [cm] = R [m] * 100.
// So kBConv = 0.3 / 100 = 0.003  (GeV/c per T per cm)
static const double kBConv = 0.003;

// Charge of the Lambda daughters (elementary-charge units)
// Proton: +1,  Pion:  -1  (for Lambda -> p pi-)
static const int kChargeProton = +1;
static const int kChargePion   = -1;

// Guard threshold: skip Lambdas whose transverse momentum fraction
// pT/|p| is smaller than this (Lambda nearly collinear with beam axis)
static const double kMinSinTheta = 1.e-5;

// ==========================================================================
// Mathematical constants and other useful constants for error bars:
// ==========================================================================
static const double Pi = TMath::Pi();
static const double TwoPi = 2.0*Pi;

// Number of independent chunks for the data-chunking error estimator.
// Each jet-group event is assigned to one chunk sequentially (round-robin).
// The final error is stddev(chunk_means) / sqrt(K_filled), which properly
// accounts for intra-event correlations.
// The whole problem in doing otherwise is that just taking an average of
// per-event averages did not seem to properly grasp the REALLY LARGE correlations
// that come from making cuts in eta, dca or min pT of the daughters!
// The no-cuts (and no eta cuts) histograms had stable error bars, but as soon
// as you start introducing cuts which create preferential regions of the phase space,
// we start seeing underestimated error bars: re-running the same code, with the same
// statistics, was making the mean <R> values fluctuate around zero with error bars
// that were not consistent of each other -- error was underestimated!!!
// The per-event average should already deal with the correlation/covariance that arises
// from (artificially) making 4 or so Lambdas share the same randomly sorted jet direction.
// (do notice that the per-event mean is not implemented as filling histograms with sumEvent/countsEvent,
//  because events with different amounts of Lambdas have different statistical weights that are counterbalanced
//  by weighting (multiplying) by countsEvent. The magic of per-event averaging comes in flushing the counters
//  once per event, which then makes the variance of this new set account for intra-event correlations.
//  See the comments in FlushEventMean if this is still not that clear here.)
// This should be a small correlation, but still a correlation. What we are trying to fix
// with this kChunks strategy is to get an error bar that grasps the correlation introduced
// from doing kinematic cuts: as I've shown (check the end of log 655), min pT cuts, DCA cuts
// and so on create regions of the phase space that were preferentially sampled, so there is
// no longer an IID constraint on the Lambdas and they have a covariance.
// At most 4 or so Lambdas per event is too little to grasp this phase-space covariance, so another "level"
// of averaging could possibly do the trick.
// At least, the goal of the kChunks is to do so. I might have missed something when thinking this through...
static const int kChunks = 200; // Increased from 50 to 200


// ==========================================================================
/**
 * @brief Computes the pT at the analytic maximum of the thermal Boltzmann mT
 *        spectrum dN/dpT ~ pT * exp(-mT/T), used as the rejection-sampling
 *        envelope in SampleLambdaPt().
 * @param T  Boltzmann temperature parameter [GeV].
 * @return   pT value at the spectrum peak [GeV/c].
 */
// ==========================================================================
double ComputeThermalPtMaximum(double T){
    double T2 = T * T;
    double m2 = kMassLambda * kMassLambda;
    double pt2 = 0.5 * (T2 + std::sqrt(T2 * T2 + 4.0 * T2 * m2));

    return std::sqrt(pt2);
}

/**
 * ==========================================================================
 * @brief Sample Lambda transverse momentum from a thermal Boltzmann mT spectrum.
 * --------------------------------------------------------------------------
 *
 * Samples pT according to dN / dpT ~ pT * exp( -mT / T ), where mT = sqrt( pT^2 + mLambda^2 )
 *
 * The sampling is performed with the accept-reject (rejection sampling) method:
 *   1. Generate a trial pT uniformly in [pTmin, pTmax]
 *   2. Evaluate the target function f(pT)
 *   3. Generate a second random number uniformly in [0, fmax]
 *   4. Accept the trial pT if the random number falls below f(pT)
 *   5. Otherwise repeat until acceptance
 *
 * The envelope maximum fmax must satisfy:
 *   fmax >= max( pT * exp(-mT/T) )
 * over the full sampling interval.
 *
 * IMPORTANT:
 *   fmax is precomputed outside this function for efficiency.
 *   Recomputing the envelope maximum for every sampled Lambda would be
 *   unnecessarily expensive.
 *
 * @param rng    TRandom3 random-number generator, passed by reference (if not passed by reference, the RNG list of already-used values would not be updated properly!)
 * @param T      Boltzmann temperature parameter [GeV]
 * @param pTmin  Lower bound of sampling interval [GeV/c]
 * @param pTmax  Upper bound of sampling interval [GeV/c]
 * @param fmax   Envelope maximum used in rejection sampling
 *
 * @return One accepted Lambda transverse momentum pT [GeV/c]
 * ==========================================================================
 */
static double SampleLambdaPt(TRandom3 &rng, double T, double pTmin, double pTmax, double fmax){
    // -- Repeat until one trial pT is accepted --
    while (true) {
        // -- Generate a trial pT uniformly inside the requested interval --
        double pt = pTmin + (pTmax - pTmin) * rng.Rndm();
        // -- Compute transverse mass corresponding to the trial pT --
        double mT = std::sqrt(pt * pt + kMassLambda * kMassLambda);
        // -- Evaluate the target Boltzmann weight at this pT --
        double f = pt * std::exp(-mT / T);

        // -- Accept the trial pT with probability f/fmax --
        //    Geometrically: accept if a random point inside the rectangle
        //    [pTmin,pTmax] x [0,fmax] falls below the target curve f(pT)
        if (rng.Rndm() * fmax < f)
            return pt;
    }
}


// ==========================================================================
/**
 * @brief Computes the analytical transverse distance of closest approach
 *        (DCA_xy) of a charged daughter track to the primary vertex at (0,0).
 *
 * A charged particle in a uniform field Bz [T] along +z, starting at
 * (xv, yv) [cm] with transverse momentum (px, py) [GeV/c], traces a circle
 * in the x-y plane.  The signed helix radius is:
 *   r_signed = pT / (q * kBConv * Bz)   [cm]
 * (q = +1 for proton, -1 for pion; Bz can be positive or negative)
 *
 * The center of the circle in the transverse plane is:
 *   Cx = xv + (py / pT) * r_signed
 *   Cy = yv - (px / pT) * r_signed
 *
 * The DCA_xy to the origin is:
 *   DCA_xy = | sqrt(Cx^2 + Cy^2) - |r_signed| |
 *
 * This expression follows from the fact that the particle moves on a circle
 * of radius |r_signed| centered at (Cx, Cy), and the closest point on that
 * circle to the origin is at distance |dist(center,origin) - |r_signed||.
 *
 * @note Only the transverse DCA is computed; the longitudinal DCA (z) is
 *       not used by this toy.
 *
 * @param xv      Decay vertex transverse x position [cm].
 * @param yv      Decay vertex transverse y position [cm].
 * @param px      Daughter transverse momentum x component [GeV/c].
 * @param py      Daughter transverse momentum y component [GeV/c].
 * @param pT      Daughter transverse momentum magnitude [GeV/c].
 * @param charge  Daughter charge in units of e (+1 proton, -1 pion).
 * @param Bz      Magnetic field z-component [T] (sign matters).
 * @return        DCA_xy [cm], or 1e9 if pT ~ 0.
 */
// ==========================================================================
static double ComputeDCAxy(double xv,  double yv,
                            double px,  double py,  double pT,
                            int    charge,
                            double Bz)
{
    // Guard against zero transverse momentum (straight-line track, infinite R)
    if (pT < 1.e-9) return 1.e9;

    // Signed radius of curvature [cm]
    double r_signed = pT / (charge * kBConv * Bz); // Positive charge in positive Bz bends clockwise in the x-y plane

    // Center of the helix circle in the transverse plane [cm]
    double Cx = xv + (py / pT) * r_signed;
    double Cy = yv - (px / pT) * r_signed;

    // Distance from the primary vertex (0, 0) to the circle center
    double dist_center = std::sqrt(Cx * Cx + Cy * Cy);

    // DCA = | dist_to_center - |radius| |
    double dca = std::fabs(dist_center - std::fabs(r_signed));
    return dca;
}


// ==========================================================================
/**
 * @brief Wraps an angle (in radians) into the interval [-pi, pi).
 * @note  Input is assumed to lie within [0, 2*pi); this allows a single
 *        conditional subtraction instead of a full fmod call.
 * @param phi  Input angle in radians, assumed in [0, 2*pi).
 * @return     Equivalent angle in [-pi, pi).
 */
// ==========================================================================
inline double wrapToPiFast(double phi){
    return (phi < Pi) ? phi : (phi - TwoPi);
}

// ==========================================================================
/**
 * @brief High-precision accumulator for integrated observables using
 *        compensated Kahan summation.
 *
 * This struct is designed as a lightweight, allocation-free replacement for
 * single-bin ROOT TProfiles inside extremely large event loops
 * (e.g. O(10^8) fills), where standard floating-point accumulation may lose
 * precision for small signals (~10^-6).
 *
 * The accumulator tracks:
 *   - N          : number of entries
 *   - sum_y      : compensated sum of values
 *   - sum_y2     : compensated sum of squared values
 *
 * Both sums use Kahan compensated summation to minimize catastrophic loss of
 * significance during long accumulation chains.
 *
 * The resulting mean and Standard Error of the Mean (SEM) reproduce the
 * statistical behavior of a 1-bin TProfile while avoiding:
 *   - ROOT virtual-function overhead in the hot loop
 *   - repeated histogram bin lookups
 *   - floating-point precision degradation from naive summation
 *
 * Intended usage:
 *   - Instantiate directly inside ScenarioHistos
 *   - Call Add(value) inside the event loop
 *   - Retrieve results via GetMean() and GetSEM()
 *   - Optionally export to a 1-bin TH1D after processing for ROOT I/O
 *
 * All methods are inline to maximize compiler optimization and vectorization
 * opportunities in performance-critical code paths.
 */
// ==========================================================================
struct KahanAccumulator {
    double sum_y  = 0.0;
    double c_y    = 0.0; // Kahan compensation for sum_y
    double sum_y2 = 0.0;
    double c_y2   = 0.0; // Kahan compensation for sum_y2
    long   N      = 0;

    // Fast inline add, meant to be used inside the hot loop
    inline void Add(double y) {
        // Kahan sum for the values
        double v_y = y - c_y;
        double t_y = sum_y + v_y;
        c_y = (t_y - sum_y) - v_y;
        sum_y = t_y;

        // Kahan sum for the squared values (for the SEM error bars)
        double y2 = y * y;
        double v_y2 = y2 - c_y2;
        double t_y2 = sum_y2 + v_y2;
        c_y2 = (t_y2 - sum_y2) - v_y2;
        sum_y2 = t_y2;

        N++;
    }

    /// @brief Returns the accumulated mean value.
    inline double GetMean() const {
        return (N > 0) ? (sum_y / (double)N) : 0.0;
    }

    /**
     * @brief Returns the Standard Error of the Mean (SEM).
     *
     * Computed using the unbiased sample variance:
     *   variance = (sum_y2 - N * mean^2) / (N - 1)
     * and:
     *   SEM = sqrt(variance / N)
     *
     * This implementation should match the statistical interpretation of a
     * ROOT TProfile bin error for integrated observables.
     */
    inline double GetSEM() const {
        if (N < 2) return 0.0;

        double mean = sum_y / (double)N;
            // Sample variance
        double variance = (sum_y2 - (double)N * mean * mean) / (double)(N - 1);
            // Numerical safeguard against tiny negative roundoff artifacts
        if (variance < 0.0) variance = 0.0;

        return std::sqrt(variance / (double)N);
    }

    /**
     * @brief Exports the accumulated statistics into a 1-bin ROOT histogram.
     *
     * The histogram bin content is filled with the accumulated mean, while
     * the bin error is set to the Standard Error of the Mean (SEM).
     *
     * @param h Pointer to the target histogram. If nullptr, the method
     *          returns immediately without performing any operation.
     */
    inline void FlushToTH1(TH1* h) const { // The method is declared as a const because it doesn't change the state of the KahanAccumulator
        if (!h) return; // Safety check
        h->SetBinContent(1, GetMean());
        h->SetBinError(1, GetSEM());
        h->SetEntries(N);
    }
};


// ==========================================================================
/**
 * @brief Holds all ROOT histogram and TProfile pointers for one combination
 *        of (cut scenario, eta half).
 *
 * Populated by BookScenario() and filled by FillScenario().
 * ROOT owns all objects once they are associated with a TDirectory via cd(),
 * so no manual deletion is required.
 */
// ==========================================================================
struct ScenarioHistos {
    // -- Main diagnostic: proton emission angles in Lambda rest frame --
    // These are the primary plots requested in the analysis note.
    TH2D*    h2d_cosTheta_phi;   // 2D: cos(theta*) vs phi*
    TH1D*    h1d_cosTheta;       // cos(theta*) projection
    TH1D*    h1d_phi;            // phi* projection

    // -- Ring observable proxy (using beam z-hat as jet direction) --
    // R_proxy = (3/alpha) * p*_D . (z x lambda_unit) / |z x lambda_unit|
    //         = -(3/alpha) * sin(theta*) * cos(phi*)
    // Directly measures the left-right asymmetry as a ring-like quantity.
    TH1D*     h1d_ringProxy;     // distribution of R_proxy (weird, but could be interesting!)
    TProfile* pRingProxy;        // integrated <R_proxy> (one bin at x=0)
    TProfile* pRingProxyVsEta;   // <R_proxy> vs Lambda pseudorapidity
    TProfile* pRingProxyVsPt;    // <R_proxy> vs Lambda pT

    // Same, but for the ring observable calculated with respect to a randomly sampled direction acting as the jet:
    // (AEE should persist even if the trigger direction is totally randomized. In a way, it depends only on the Lambdas,
    //  not on the jets' direction)
    TH1D*     h1d_ringProxyJet;     // distribution of R_proxyJet values (weird, but could be interesting!)
    TProfile* pRingProxyJet;        // integrated <R_proxyJet> (one bin at x=0)
    TProfile* pRingProxyJetVsEta;   // <R_proxyJet> vs Lambda pseudorapidity
    TProfile* pRingProxyJetVsEtaJet;// <R_proxyJet> vs Jet pseudorapidity
    TProfile* pRingProxyJetVsPt;    // <R_proxyJet> vs Lambda pT

    // Integrated <R_proxyJet> split by jet pseudorapidity sign.
    // Fills happen conditionally inside FillScenario based on eta_jet.
    // Allows the same integrated-ring study as MakeFig12 but binned by jet eta
    // instead of Lambda eta, directly revealing the eta_jet dependence.
    TProfile* pRingProxyJet_JetEtaPos; // <R_proxyJet> for eta_jet >= 0, all Lambda eta
    TProfile* pRingProxyJet_JetEtaNeg; // <R_proxyJet> for eta_jet <  0, all Lambda eta

    // Some histograms for a better error estimation than the SEM:
      // When we choose nLambdasSinceJetShuffle == 3, we are implicitly creating a correlation between 
      // four Lambdas, as if they were in the same event. That can make the SEM error bars underestimate
      // the error, which may make it look like we have non-zero spurious-signal even after integrating
      // everything in eta.
    // Per-event mean estimator for the integrated ring-jet observable:
    TH1D*  hEventMeanRingProxyJet; // one fill per event (jet group); mean +/- error = unbiased estimator
    TProfile* pEventMeanRingProxyJetIntegrated;
    double evtSumRpj  = 0.;  // running sum of ringProxyJet values in the current event
    int    evtCntRpj  = 0;   // number of entries accumulated so far in each event
    long   evtIdxRpj  = 0;   // total events seen; drives chunk assignment (evtIdxRpj % kChunks)

    // Data-chunking accumulators (one slot per chunk):
    // After the main loop, FinalizeChunks() computes mu_k = chunkSumRpj[k]/chunkCntRpj[k]
    // for every non-empty chunk and fills hChunkMeansRingProxyJet.
    // hChunkMeansRingProxyJet->GetMeanError() then gives stddev(mu_k)/sqrt(K_filled),
    // which is the (with luck!) bulletproof error estimate the chunking method was designed to produce.
    double chunkSumRpj[kChunks];    // sum of all accepted ringProxyJet values per chunk
    long   chunkCntRpj[kChunks];    // count of all accepted entries/Lambdas per chunk
    TH1D*  hChunkMeansRingProxyJet; // K fills (one per non-empty chunk); GetMeanError() = final error
    TProfile* pChunkMeansRingProxyJetIntegrated; // A single pre-integrated observable. Should be easier to read on the final .root.

    // -- Decay geometry and daughter kinematics (useful for diagnosis) --
    TH1D*    h1d_decayRadius;    // transverse decay vertex radius [cm]
    TH1D*    h1d_pT_proton;      // proton pT in lab frame [GeV/c]
    TH1D*    h1d_pT_pion;        // pion pT in lab frame [GeV/c]
    TH1D*    h1d_DCA_proton;     // proton DCA_xy to PV [cm]
    TH1D*    h1d_DCA_pion;       // pion DCA_xy to PV [cm]

    // -- Lambda kinematics --
    TH1D*    h1d_pT_lambda;      // Lambda pT [GeV/c]
    TH1D*    h1d_eta_lambda;     // Lambda pseudorapidity

    // -- Proton rest-frame direction vector field vs Lambda transverse momentum --
    // The three TProfile2Ds store <p*_x>, <p*_y>, <p*_z> in bins of (px_lam, py_lam).
    // Together they form a 2D vector field in transverse momentum space.
    // The expected ring-like structure from DCA cuts appears as arrows tangent to
    // circles of constant pT, pointing in the direction of the average polarisation
    // bias induced by the magnetic-field-dependent acceptance.
    // pPstarZ_vsPxPy also serves as the colormap channel for the TProfile2D plots (out-of-plane component).
    TProfile2D* pPstarX_vsPxPy;   // <p*_x> vs (px_lam, py_lam)
    TProfile2D* pPstarY_vsPxPy;   // <p*_y> vs (px_lam, py_lam)
    TProfile2D* pPstarZ_vsPxPy;   // <p*_z> vs (px_lam, py_lam)  [used as colormap -- cannot fit anything else in 2D!]

    // 1D compact version: <p*_x,y,z> vs azimuthal angle phi_lam.
    // Directly shows the sinusoidal modulation around the ring without
    // needing any fancy arrow-drawing code that we do in plotHelicityEfficiency.cxx
    TProfile*   pPstarX_vsPhiLam; // <p*_x> vs phi_lam
    TProfile*   pPstarY_vsPhiLam; // <p*_y> vs phi_lam
    TProfile*   pPstarZ_vsPhiLam; // <p*_z> vs phi_lam

    // -- Proton rest-frame direction vector field vs Lambda (pz, px) -- ZX plane --
    // (Wrote it in the ZX plane because that is the most usual projection of the detector!)
    // Captures the longitudinal structure: allows to visualize if the fake polarisation pattern
    // has a component correlated with both the transverse and beam directions
    // Arrow components: <p*_z> (horizontal, along beam) and <p*_x> (vertical).
    // Colormap:         <p*_y> (out-of-plane for the ZX cut).
    TProfile2D* pPstarX_vsPzPx;   // <p*_x> vs (pz_lam, px_lam)  [arrow y-component]
    TProfile2D* pPstarZ_vsPzPx;   // <p*_z> vs (pz_lam, px_lam)
    TProfile2D* pPstarY_vsPzPx;   // <p*_y> vs (pz_lam, px_lam)  [colormap]

    // 1D compact version vs Lambda pseudorapidity (natural 1D variable for ZX,
    // since pz = pT * sinh(eta) and the pattern is driven by longitudinal kinematics).
    TProfile*   pPstarX_vsEtaLam; // <p*_x> vs eta_lam
    TProfile*   pPstarY_vsEtaLam; // <p*_y> vs eta_lam
    TProfile*   pPstarZ_vsEtaLam; // <p*_z> vs eta_lam

    // Manual accumulators (for Kahan sums):
    KahanAccumulator kRingProxy;
    KahanAccumulator kRingProxyJet;
    KahanAccumulator kRingProxyJet_JetEtaPos;
    KahanAccumulator kRingProxyJet_JetEtaNeg;

    // Histograms to flush the Kahan sums after they have finished accumulating everything,
    // just to store in the final .root. The "precision computation" part should already be
    // dealt with at this point!
    // I am using TH1Ds instead of, say, TProfiles, because TH1D are essentially a dumb visual container here.
    // TProfile might want to do something else. I can fill them with a simple "FlushToTH1" call with 3 lines.
    TH1D* hRingProxy_Kahan;
    TH1D* hRingProxyJet_Kahan;
    TH1D* hRingProxyJet_JetEtaPos_Kahan;
    TH1D* hRingProxyJet_JetEtaNeg_Kahan;
};


// ==========================================================================
/**
 * @brief Holds all four cut-scenario groups (NoCuts, pTCutOnly, DCACutOnly,
 *        BothCuts) for one histogram family, each split into EtaPos, EtaNeg,
 *        and All sub-directories.
 *
 * Two instances of this struct are used in the main function: one for the
 * WithoutEtaGate family and one for the WithEtaGate family.
 * Populated by CreateFamily() and filled by FillFamily().
 */
// ==========================================================================
struct FamilyHistos {
    ScenarioHistos NC_Pos, NC_Neg, NC_All;  // NoCuts
    ScenarioHistos PT_Pos, PT_Neg, PT_All;  // pTCutOnly
    ScenarioHistos DC_Pos, DC_Neg, DC_All;  // DCACutOnly
    ScenarioHistos BC_Pos, BC_Neg, BC_All;  // BothCuts
};


// ==========================================================================
/**
 * @brief Allocates and initialises all histograms inside the given TDirectory.
 *
 * The directory must already exist.  Calls dir->cd() before creating objects
 * so ROOT associates them with the correct directory.
 *
 * @param dir  Pointer to an existing TDirectory in which histograms are booked.
 * @param etaMaxDetector  Maximum acceptance allowed for Toy Model lambda daughters (and jets)
 * @return     A filled ScenarioHistos struct with all histogram pointers set.
 */
// ==========================================================================
static ScenarioHistos BookScenario(TDirectory* dir, double etaMaxDetector)
{
    ScenarioHistos h;

    dir->cd();  // All subsequent "new TH1D/TH2D" calls attach to this dir

    // ---------- Main diagnostic: (cos theta*, phi*) space ----------
    // Expected: flat before cuts; asymmetric after cuts.
    h.h2d_cosTheta_phi = new TH2D("h2d_cosTheta_phi", "Proton emission in Lambda frame;"
        "cos(#theta*) [proton vs #Lambda dir];"
        "#phi* [proton azimuth around #Lambda] [rad]",
        50, -1., 1.,          // 50 bins in cos(theta*) in [-1, 1]
        64, -Pi, Pi); // 64 bins in phi* in [-pi, pi]
    h.h1d_cosTheta = new TH1D("h1d_cosTheta", "cos(#theta*) -- proton vs #Lambda direction; cos(#theta*);Counts", 50, -1., 1.);
    h.h1d_phi = new TH1D("h1d_phi", "#phi* -- proton azimuth around #Lambda axis; #phi* [rad];Counts",64, -Pi, Pi);

    // ---------- Ring observable proxy ----------
    // Range is set to [-3/alpha, +3/alpha] which covers the full prefactor.
    double rmax = kPolPrefactor * 1.05; // Slight margin beyond max value
    h.h1d_ringProxy = new TH1D("h1d_ringProxy", "Ring observable proxy (z-hat as jet direction); R_{proxy};Counts", 120, -rmax, rmax);
    h.pRingProxy = new TProfile("pRingProxy", "Integrated <R_{proxy}> (single-bin TProfile); bin;<R_{proxy}>", 1, -0.5, 0.5);
    h.pRingProxyVsEta = new TProfile("pRingProxyVsEta", "<R_{proxy}> vs #Lambda pseudorapidity; #eta_{#Lambda};<R_{proxy}>", 9, -etaMaxDetector, etaMaxDetector); // 18 bins, 0.1 wide each if etaMaxDetector = 0.9
    h.pRingProxyVsPt = new TProfile("pRingProxyVsPt", "<R_{proxy}> vs #Lambda p_{T}; p_{T}^{#Lambda} [GeV/c];<R_{proxy}>", 10, 0., 5.); // 20 bins of 0.25 GeV/c each

    // ---------- Ring observable proxy for randomly sampled jet direction ----------
    h.h1d_ringProxyJet = new TH1D("h1d_ringProxyJet", "Ring observable proxy (random jet direction); R_{proxyJet};Counts", 120, -rmax, rmax);
    h.pRingProxyJet = new TProfile("pRingProxyJet", "Integrated <R_{proxyJet}> (single-bin TProfile); bin;<R_{proxyJet}>", 1, -0.5, 0.5);
    h.pRingProxyJetVsEta = new TProfile("pRingProxyJetVsEta", "<R_{proxyJet}> vs #Lambda pseudorapidity; #eta_{#Lambda};<R_{proxyJet}>", 9, -etaMaxDetector, etaMaxDetector);
    h.pRingProxyJetVsEtaJet = new TProfile("pRingProxyJetVsEtaJet", "<R_{proxyJet}> vs Jet; #eta_{Jet};<R_{proxyJet}>", 9, -etaMaxDetector, etaMaxDetector);
    h.pRingProxyJetVsPt = new TProfile("pRingProxyJetVsPt", "<R_{proxyJet}> vs #Lambda p_{T}; p_{T}^{#Lambda} [GeV/c];<R_{proxyJet}>", 10, 0., 5.); // 20 bins of 0.25 GeV/c each

        // Adding a study for the dependency with jet eta:
    h.pRingProxyJet_JetEtaPos = new TProfile("pRingProxyJet_JetEtaPos", "Integrated <R_{proxyJet}> for #eta_{jet} #geq 0; bin; <R_{proxyJet}>", 1, -0.5, 0.5);
    h.pRingProxyJet_JetEtaNeg = new TProfile("pRingProxyJet_JetEtaNeg", "Integrated <R_{proxyJet}> for #eta_{jet} < 0; bin; <R_{proxyJet}>", 1, -0.5, 0.5);
    
    // For the per-event mean estimators for better error bars, considering correlations:
    double evRmax = kPolPrefactor * 1.05; // This proved to be an OK estimator
    h.hEventMeanRingProxyJet = new TH1D("hEventMeanRingProxyJet", "Event-mean <R_{proxyJet}> per jet group; #bar{R}_{proxyJet} per event; Events", 200, -evRmax, evRmax);
    h.pEventMeanRingProxyJetIntegrated = new TProfile("pEventMeanRingProxyJet", "Mean of all per-event means <R_{proxyJet}>; Dummy variable; Global mean <R_{proxyJet}>", 1, 0., 1.);
    h.evtSumRpj = 0.;
    h.evtCntRpj = 0;
    h.evtIdxRpj = 0;
    for (int k = 0; k < kChunks; ++k) {
        h.chunkSumRpj[k] = 0.; h.chunkCntRpj[k] = 0;
    }
    double evRmaxChunks = kPolPrefactor * 0.05; // Instead of being in the order of |<R>| = 4, the per-chunk averages should already be close to the 10^(-3) magnitude
    h.hChunkMeansRingProxyJet = new TH1D("hChunkMeansRingProxyJet", Form("Chunk means of <R_{proxyJet}> (%d chunks); #bar{R}_{proxyJet} per chunk; Chunks", kChunks), 100, -evRmaxChunks, evRmaxChunks);
    h.pChunkMeansRingProxyJetIntegrated = new TProfile("pChunkMeansRingProxyJetIntegrated", Form("Mean of chunks of <R_{proxyJet}> (%d chunks); Dummy variable; Chunks", kChunks), 1, 0., 1.);

    // ---------- Decay vertex and daughter kinematics ----------
    h.h1d_decayRadius = new TH1D("h1d_decayRadius", "Transverse decay vertex radius; r_{decay} = #sqrt{x_{v}^{2}+y_{v}^{2}} [cm];Counts", 180, 0., 180.); // 0 to 150 cm (would cover most of ITS/TPC inner field cage)
    h.h1d_pT_proton = new TH1D("h1d_pT_proton", "Proton p_{T} in lab frame; p_{T}^{proton} [GeV/c];Counts", 100, 0., 5.);
    h.h1d_pT_pion = new TH1D("h1d_pT_pion", "Pion p_{T} in lab frame; p_{T}^{#pi} [GeV/c];Counts", 100, 0., 5.);
    h.h1d_DCA_proton = new TH1D("h1d_DCA_proton", "Proton DCA_{xy} to primary vertex; DCA_{xy}^{proton} [cm];Counts", 200, 0., 10.); // Fine binning near 0 to see the cut clearly
    h.h1d_DCA_pion = new TH1D("h1d_DCA_pion", "Pion DCA_{xy} to primary vertex; DCA_{xy}^{#pi} [cm];Counts", 200, 0., 10.);
    // ---------- Lambda kinematics ----------
    h.h1d_pT_lambda = new TH1D("h1d_pT_lambda", "#Lambda p_{T}; p_{T}^{#Lambda} [GeV/c];Counts", 100, 0., 10.);
    h.h1d_eta_lambda = new TH1D("h1d_eta_lambda", "#Lambda pseudorapidity; #eta_{#Lambda};Counts", 36, -etaMaxDetector, etaMaxDetector);

    // Vector field: average proton rest-frame direction in (px, py) space.
    // Range [-3, 3] GeV/c covers >95% of the thermal spectrum (T~0.3 GeV).
    // 40x40 bins give cells of ~0.15 GeV/c width, fine enough to resolve the ring but coarse enough to have good statistics per cell at 1M+ Lambdas.
    // Notice that we can do a grid in (px,py) because that also relates finely with actual XY coordinates given a transverse travel time
    // for the Lambda before decaying! If we actually want to plot this in spatial coordinates, just take (px,py) and multiply by the sampled
    // Tau for each Lambda.
    const double pXYmax = 3.0;
    const int    pXYbin = 40;
    h.pPstarX_vsPxPy = new TProfile2D("pPstarX_vsPxPy", "<p*_{x}> vs (p_{x}^{#Lambda}, p_{y}^{#Lambda});" "p_{x}^{#Lambda} [GeV/c];p_{y}^{#Lambda} [GeV/c];<p*_{x}>", pXYbin, -pXYmax, pXYmax,  pXYbin, -pXYmax, pXYmax);
    h.pPstarY_vsPxPy = new TProfile2D("pPstarY_vsPxPy", "<p*_{y}> vs (p_{x}^{#Lambda}, p_{y}^{#Lambda});" "p_{x}^{#Lambda} [GeV/c];p_{y}^{#Lambda} [GeV/c];<p*_{y}>", pXYbin, -pXYmax, pXYmax,  pXYbin, -pXYmax, pXYmax);
    h.pPstarZ_vsPxPy = new TProfile2D("pPstarZ_vsPxPy", "<p*_{z}> vs (p_{x}^{#Lambda}, p_{y}^{#Lambda});" "p_{x}^{#Lambda} [GeV/c];p_{y}^{#Lambda} [GeV/c];<p*_{z}>", pXYbin, -pXYmax, pXYmax,  pXYbin, -pXYmax, pXYmax);

    // 1D phi profile: 32 bins of width ~0.2 rad, clean for a sinusoidal fit later if wanted
    h.pPstarX_vsPhiLam = new TProfile("pPstarX_vsPhiLam", "<p*_{x}> vs #phi_{#Lambda}; #phi_{#Lambda} [rad]; <p*_{x}>", 32, 0., TwoPi);
    h.pPstarY_vsPhiLam = new TProfile("pPstarY_vsPhiLam", "<p*_{y}> vs #phi_{#Lambda}; #phi_{#Lambda} [rad]; <p*_{y}>", 32, 0., TwoPi);
    h.pPstarZ_vsPhiLam = new TProfile("pPstarZ_vsPhiLam", "<p*_{z}> vs #phi_{#Lambda}; #phi_{#Lambda} [rad]; <p*_{z}>", 32, 0., TwoPi);

    // ZX vector field:
    // px range same as XY [-3,3]; pz range wider [-4,4] since pz = pT*sinh(eta)
    // and Lambda eta is larger than daughter eta.
    const double pXmax = 3.0, pZmax = 4.0;
    const int    pZXbin = 40;
    h.pPstarX_vsPzPx = new TProfile2D("pPstarX_vsPzPx", ";<p*_{x}> vs (p_{z}^{#Lambda}, p_{x}^{#Lambda});" "p_{z}^{#Lambda} [GeV/c];p_{x}^{#Lambda} [GeV/c];<p*_{x}>", pZXbin, -pZmax, pZmax, pZXbin, -pXmax, pXmax);
    h.pPstarZ_vsPzPx = new TProfile2D("pPstarZ_vsPzPx", ";<p*_{z}> vs (p_{z}^{#Lambda}, p_{x}^{#Lambda});" "p_{z}^{#Lambda} [GeV/c];p_{x}^{#Lambda} [GeV/c];<p*_{z}>", pZXbin, -pZmax, pZmax, pZXbin, -pXmax, pXmax);
    h.pPstarY_vsPzPx = new TProfile2D("pPstarY_vsPzPx", ";<p*_{y}> vs (p_{z}^{#Lambda}, p_{x}^{#Lambda});" "p_{z}^{#Lambda} [GeV/c];p_{x}^{#Lambda} [GeV/c];<p*_{y}>", pZXbin, -pZmax, pZmax, pZXbin, -pXmax, pXmax);

    // 1D eta profiles (18 bins, matching pRingProxyVsEta width)
    h.pPstarX_vsEtaLam = new TProfile("pPstarX_vsEtaLam", "<p*_{x}> vs #eta_{#Lambda}; #eta_{#Lambda}; <p*_{x}>", 18, -etaMaxDetector, etaMaxDetector);
    h.pPstarY_vsEtaLam = new TProfile("pPstarY_vsEtaLam", "<p*_{y}> vs #eta_{#Lambda}; #eta_{#Lambda}; <p*_{y}>", 18, -etaMaxDetector, etaMaxDetector);
    h.pPstarZ_vsEtaLam = new TProfile("pPstarZ_vsEtaLam", "<p*_{z}> vs #eta_{#Lambda}; #eta_{#Lambda}; <p*_{z}>", 18, -etaMaxDetector, etaMaxDetector);

    // Histograms to store the final Kahan accumulator results:
    h.hRingProxy_Kahan = new TH1D("hRingProxy_Kahan", "Integrated <R_{proxy}> (Kahan Acc.); bin;<R_{proxy}>", 1, -0.5, 0.5);
    h.hRingProxyJet_Kahan = new TH1D("hRingProxyJet_Kahan", "Integrated <R_{proxyJet}> (Kahan Acc.); bin;<R_{proxyJet}>", 1, -0.5, 0.5);
    h.hRingProxyJet_JetEtaPos_Kahan = new TH1D("hRingProxyJet_JetEtaPos_Kahan", "Integrated <R_{proxyJet}> (Kahan Acc.) for #eta_{jet} #geq 0; bin; <R_{proxyJet}>", 1, -0.5, 0.5);
    h.hRingProxyJet_JetEtaNeg_Kahan = new TH1D("hRingProxyJet_JetEtaNeg_Kahan", "Integrated <R_{proxyJet}> (Kahan Acc.) for #eta_{jet} < 0; bin; <R_{proxyJet}>", 1, -0.5, 0.5);

    return h;
}


// ==========================================================================
/**
 * @brief Fills all histograms in a ScenarioHistos struct for one Lambda
 *        candidate.  Called once per Lambda per cut scenario (up to 4 times).
 *
 * All arguments are computed in the main event loop; see variable names there.
 *
 * @param h           Reference to the ScenarioHistos struct to fill.
 * @param cosTheta    cos(theta*) = p*_proton_unit . e1.
 * @param phi         phi* = atan2(p*_proton.e3, p*_proton.e2) [rad].
 * @param ringProxy   R_proxy = (3/alpha) * p*_D . (z x lambda_unit)/||z x lambda_unit|| .
 * @param ringProxyJet Same as R_proxy, but the jet direction takes the place of the z axis
 * @param decayR      Transverse decay vertex radius [cm].
 * @param pT_proton   Proton pT in the lab frame [GeV/c].
 * @param pT_pion     Pion pT in the lab frame [GeV/c].
 * @param dca_proton  Proton DCA_xy to PV [cm].
 * @param dca_pion    Pion DCA_xy to PV [cm].
 * @param pT_lambda   Lambda pT [GeV/c].
 * @param eta_lambda  Lambda pseudorapidity.
 * @param eta_jet      Randomly sampled jet's pseudorapidity.
 */
// ==========================================================================
static void FillScenario(ScenarioHistos& h,
                          double cosTheta,
                          double phi,
                          double ringProxy,
                          double ringProxyJet,
                          double decayR,
                          double pT_proton,
                          double pT_pion,
                          double dca_proton,
                          double dca_pion,
                          double pT_lambda,
                          double eta_lambda,
                          double eta_jet,
                          double pstar_x,
                          double pstar_y,
                          double pstar_z,
                          double px_lam,
                          double py_lam,
                          double phi_lam,
                          double pz_lam)
{
    // Main diagnostic histograms
    h.h2d_cosTheta_phi->Fill(cosTheta, phi);
    h.h1d_cosTheta->Fill(cosTheta);
    h.h1d_phi->Fill(phi);

    // Ring proxy
    h.h1d_ringProxy->Fill(ringProxy);
    h.pRingProxy->Fill(0., ringProxy); // integrated average
    h.kRingProxy.Add(ringProxy); // Integrated average via Kahan summation
    h.pRingProxyVsEta->Fill(eta_lambda, ringProxy);
    h.pRingProxyVsPt->Fill(pT_lambda, ringProxy);

    h.h1d_ringProxyJet->Fill(ringProxyJet);
    h.pRingProxyJet->Fill(0., ringProxyJet); // integrated average
    h.kRingProxyJet.Add(ringProxyJet); // Integrated average via Kahan summation
    h.pRingProxyJetVsEta->Fill(eta_lambda, ringProxyJet);
    h.pRingProxyJetVsEtaJet->Fill(eta_jet, ringProxyJet);
    h.pRingProxyJetVsPt->Fill(pT_lambda, ringProxyJet);
    // Jet-eta-split integrated profiles
    if (eta_jet >= 0.) {
        h.pRingProxyJet_JetEtaPos->Fill(0., ringProxyJet);
        h.kRingProxyJet_JetEtaPos.Add(ringProxyJet);
    } else {
        h.pRingProxyJet_JetEtaNeg->Fill(0., ringProxyJet);
        h.kRingProxyJet_JetEtaNeg.Add(ringProxyJet);
    }

    // Accumulate for the event-mean estimator; flush happens at jet reshuffle:
    h.evtSumRpj += ringProxyJet;
    h.evtCntRpj++;

    // Kinematics
    h.h1d_decayRadius->Fill(decayR);
    h.h1d_pT_proton->Fill(pT_proton);
    h.h1d_pT_pion->Fill(pT_pion);
    h.h1d_DCA_proton->Fill(dca_proton);
    h.h1d_DCA_pion->Fill(dca_pion);
    h.h1d_pT_lambda->Fill(pT_lambda);
    h.h1d_eta_lambda->Fill(eta_lambda);

    // Vector field fills
    h.pPstarX_vsPxPy->Fill(px_lam, py_lam, pstar_x);
    h.pPstarY_vsPxPy->Fill(px_lam, py_lam, pstar_y);
    h.pPstarZ_vsPxPy->Fill(px_lam, py_lam, pstar_z);
    h.pPstarX_vsPhiLam->Fill(phi_lam, pstar_x);
    h.pPstarY_vsPhiLam->Fill(phi_lam, pstar_y);
    h.pPstarZ_vsPhiLam->Fill(phi_lam, pstar_z);

    // ZX vector field fills
    h.pPstarX_vsPzPx->Fill(pz_lam, px_lam, pstar_x);
    h.pPstarZ_vsPzPx->Fill(pz_lam, px_lam, pstar_z);
    h.pPstarY_vsPzPx->Fill(pz_lam, px_lam, pstar_y);
    h.pPstarX_vsEtaLam->Fill(eta_lambda, pstar_x);
    h.pPstarY_vsEtaLam->Fill(eta_lambda, pstar_y);
    h.pPstarZ_vsEtaLam->Fill(eta_lambda, pstar_z);
}


// ==========================================================================
/**
 * @brief Creates a three-level directory structure inside a parent TDirectory
 *        and books ScenarioHistos in each leaf directory.
 *
 * Directory layout created under @p parent:
 *   scenarioName/EtaPos/
 *   scenarioName/EtaNeg/
 *   scenarioName/All/
 *
 * Results are returned via the output references hPos, hNeg, and hAll.
 *
 * @note TFile inherits from TDirectory, so outFile can be passed directly
 *       as @p parent when booking at the file root level.
 *
 * @param parent        Pointer to the parent TDirectory under which the
 *                      scenario directory is created.
 * @param scenarioName  Name of the scenario sub-directory to create.
 * @param hPos          Output: ScenarioHistos booked for eta_Lambda >= 0.
 * @param hNeg          Output: ScenarioHistos booked for eta_Lambda < 0.
 * @param hAll          Output: ScenarioHistos booked for both eta halves combined.
 * @param etaMaxDetector  Maximum acceptance allowed for Toy Model lambda daughters (and jets)
 */
// ==========================================================================
static void CreateSubdirs(TDirectory*        parent,
                           const char*        scenarioName,
                           ScenarioHistos&    hPos,
                           ScenarioHistos&    hNeg,
                           ScenarioHistos&    hAll,
                           double etaMaxDetector)
{
    // Create the top-level scenario directory
    TDirectory* dirScen = parent->mkdir(scenarioName);

    // Create sub-directories for each eta half and the combined sample
    TDirectory* dirPos = dirScen->mkdir("EtaPos");
    TDirectory* dirNeg = dirScen->mkdir("EtaNeg");
    TDirectory* dirAll = dirScen->mkdir("All");

    // Book histograms inside each sub-directory
    hPos = BookScenario(dirPos, etaMaxDetector);
    hNeg = BookScenario(dirNeg, etaMaxDetector);
    hAll = BookScenario(dirAll, etaMaxDetector);

    // Return to parent directory so subsequent mkdir calls start from there
    parent->cd();
}


// ==========================================================================
/**
 * @brief Creates all four cut-scenario sub-directories (NoCuts, pTCutOnly,
 *        DCACutOnly, BothCuts) under @p parent and books all histograms into
 *        the provided FamilyHistos struct.
 *
 * @param parent  Pointer to the parent TDirectory (e.g. the WithEtaGate or WithoutEtaGate top-level directory).
 * @param f       Output: FamilyHistos struct whose 12 ScenarioHistos members will be populated.
 * @param etaMaxDetector  Maximum acceptance allowed for Toy Model lambda daughters (and jets)
 */
// ==========================================================================
static void CreateFamily(TDirectory* parent, FamilyHistos& f, double etaMaxDetector)
{
    CreateSubdirs(parent, "NoCuts",     f.NC_Pos, f.NC_Neg, f.NC_All, etaMaxDetector);
    CreateSubdirs(parent, "pTCutOnly",  f.PT_Pos, f.PT_Neg, f.PT_All, etaMaxDetector);
    CreateSubdirs(parent, "DCACutOnly", f.DC_Pos, f.DC_Neg, f.DC_All, etaMaxDetector);
    CreateSubdirs(parent, "BothCuts",   f.BC_Pos, f.BC_Neg, f.BC_All, etaMaxDetector);
}


// ==========================================================================
/**
 * @brief Fills all cut scenarios in a FamilyHistos struct for one Lambda
 *        candidate, applying the four cut combinations consistently.
 *
 * NoCuts is always filled.  The remaining three scenarios are gated on the
 * corresponding cut booleans.  Within each scenario the correct eta half
 * (EtaPos or EtaNeg) and the combined All directory are both filled.
 *
 * @param f            Reference to the FamilyHistos struct to fill.
 * @param passPt       True if both daughters pass the minimum-pT cut.
 * @param passDca      True if both daughters pass the minimum-DCA cut.
 * @param etaPos       True if eta_Lambda >= 0 (routes fill to EtaPos vs EtaNeg).
 * @param cosTheta     cos(theta*) = p*_proton_unit . e1.
 * @param phi_star     phi* = atan2(p*_proton.e3, p*_proton.e2) [rad].
 * @param ringProxy    R_proxy = (3/alpha) * p*_D . (z x lambda_unit)/|...| .
 * @param ringProxyJet Same as R_proxy, but the jet direction takes the place of the z axis
 * @param decayR       Transverse decay vertex radius [cm].
 * @param pT_p         Proton pT in the lab frame [GeV/c].
 * @param pT_pi        Pion pT in the lab frame [GeV/c].
 * @param dca_proton   Proton DCA_xy to the primary vertex [cm].
 * @param dca_pion     Pion DCA_xy to the primary vertex [cm].
 * @param pT_lam       Lambda pT [GeV/c].
 * @param eta_lam      Lambda pseudorapidity.
 * @param eta_jet      Randomly sampled jet's pseudorapidity.
 */
// ==========================================================================
static void FillFamily(FamilyHistos& f,
                        bool   passPt,
                        bool   passDca,
                        bool   etaPos,
                        double cosTheta,
                        double phi_star,
                        double ringProxy,
                        double ringProxyJet,
                        double decayR,
                        double pT_p,
                        double pT_pi,
                        double dca_proton,
                        double dca_pion,
                        double pT_lam,
                        double eta_lam,
                        double eta_jet,
                        double pstar_x,
                        double pstar_y,
                        double pstar_z,
                        double px_lam,
                        double py_lam,
                        double phi_lam,
                        double pz_lam)
{
    // Local lambda that fills one scenario's eta-half and All subdirectory.
    // Captures all per-event quantities by reference so we only pass the
    // three ScenarioHistos pointers that change between scenarios.
    auto fill = [&](ScenarioHistos& hPos, ScenarioHistos& hNeg, ScenarioHistos& hAll) {
        ScenarioHistos& hEta = etaPos ? hPos : hNeg;
        FillScenario(hEta,  cosTheta, phi_star, ringProxy, ringProxyJet, decayR, pT_p, pT_pi, dca_proton, dca_pion, pT_lam, eta_lam, eta_jet, pstar_x, pstar_y, pstar_z, px_lam, py_lam, phi_lam, pz_lam);
        FillScenario(hAll,  cosTheta, phi_star, ringProxy, ringProxyJet, decayR, pT_p, pT_pi, dca_proton, dca_pion, pT_lam, eta_lam, eta_jet, pstar_x, pstar_y, pstar_z, px_lam, py_lam, phi_lam, pz_lam);
    };

    // Scenario 1: No cuts (always filled -- serves as a flat-distribution
    //             bug-check; any asymmetry here is a code error)
    fill(f.NC_Pos, f.NC_Neg, f.NC_All);

    // Scenario 2: pT cut only
    if (passPt)
        fill(f.PT_Pos, f.PT_Neg, f.PT_All);

    // Scenario 3: DCA cut only
    if (passDca)
        fill(f.DC_Pos, f.DC_Neg, f.DC_All);

    // Scenario 4: Both cuts
    if (passPt && passDca)
        fill(f.BC_Pos, f.BC_Neg, f.BC_All);
}

// ==========================================================================
/**
 * @brief Closes the current implicit event in one ScenarioHistos accumulator
 *        and records its mean ringProxyJet value.
 *
 * @details
 * Lambdas are grouped into implicit "events" by the jet-reshuffle logic in
 * the main loop: all Lambdas that share the same random jet direction belong
 * to the same event.  While those Lambdas are being processed, their
 * ringProxyJet values are accumulated in @p h.evtSumRpj / @p h.evtCntRpj
 * rather than being flushed immediately.
 *
 * When the jet direction is about to change (or after the main loop ends),
 * this function computes the event mean R_e = evtSumRpj / evtCntRpj and
 * fills it into hEventMeanRingProxyJet.  The accumulators are then zeroed
 * so the next event starts clean.
 *
 * The spread of hEventMeanRingProxyJet across many events provides an
 * uncertainty estimate that correctly accounts for (ARTIFICIAL!) intra-event
 * correlations, unlike the plain SEM returned by TProfile which assumes all
 * entries are independent:
 *   sigma(R_bar) = hEventMeanRingProxyJet->GetMeanError() = StdDev(R_e) / sqrt(N_events)
 *
 * Nothing is performed if no entries were accumulated (e.g. all Lambdas in the event were
 * rejected by cuts).
 *
 * @param h  Reference to the ScenarioHistos whose accumulator is to be flushed.
 */
// ==========================================================================
static void FlushEventMean(ScenarioHistos& h) {
    // Chunk index is determined BEFORE the guard, so it advances for every
    // jet-group event regardless of how many Lambdas survived cuts in this
    // scenario. This keeps the chunk assignment in sync across all twelve
    // ScenarioHistos that are flushed simultaneously by FlushEventMeansFamily.
    // Possibly not the best choice, but this makes it so we could handle everything
    // with just a single index!
    int chunkIdx = (int)(h.evtIdxRpj % (long)kChunks);

    if (h.evtCntRpj > 0) { // Only fill the scenarios where at least 1 Lambda passed the cuts
        // Trying to implement weights based on the number of accepted particles per event:
            // The no-cuts (with NO ETA CUT as well! That is a cut!) distribution of averages should look Gaussian.
            // The distributions after each cut may be weirder combinations that are non-Gaussian (leptokurtic distributions, for instance).
            // What we are trying to do is make this mean of in-event means properly grasp the covariance of the ring observable and then
            // have an SEM error that is not underestimated as it was before.
        double eventMean = h.evtSumRpj / h.evtCntRpj;
        double weight = h.evtCntRpj; // More particles = higher statistical weight. We use the number of particles instead of, say, sqrt(counts)
                                     // as the weight because this converges to the grand mean over all Lambdas (instead of a mean in each event, then
                                     // a mean over events) and is an unbiased estimator.
        h.hEventMeanRingProxyJet->Fill(eventMean, weight); // In the grand mean, this has the same contribution as if we considered all 4 lambdas of
                                                           // each event separately, but notice that the variance of these means is different: it may
                                                           // be larger than the variance that would appear if we filled per-particle instead of per-event
                                                           // because we are now considering the covariance/"correlation" of having the same jet direction
                                                           // for all particles in an "event".
                                                           // In other words, the distribution of event means generally has a different variance than the
                                                           // per-particle distribution because particles within the same event share common
                                                           // event-level quantities (such as the jet direction), introducing covariance terms between
                                                           // particles in the same event.
        h.pEventMeanRingProxyJetIntegrated->Fill(0., eventMean, weight); // Fill the dummy bin with eventMean as the value and give it a weight based
                                                                        // on counts. This should be the exact same as the TH1 if it indeed does report
                                                                        // a running average for GetMean and GetMeanError instead of providing us a binned
                                                                        // average and a binned average error, which would be bad here where the binning
                                                                        // is much coarser than the true average. This is a DEBUG plot, just to make sure
                                                                        // everything is working as it should in regular TH1s.
                                                                        // Indeed, this does seem to be the case from the TStatBoxes after running.

        // Accumulate raw (unweighted) sums into the chunk:
          // (This also accumulates on a per-event average basis! The event is the smallest statistical entity here: all Lambdas in the same event are
          //  guaranteed to have the same jet direction, and that is correlation. The kChunks strategy is then meant to grab another type of covariance/
          //  correlation that was not being grasped yet: the correlation that arises from making kinematic cuts that select preferential phase space regions!)
          // Using raw sums avoids the weighting paradox: mu_k will be computed as a plain ratio (total R in chunk) / (total entries in chunk),
          // which is exactly the arithmetic mean over all Lambdas in that chunk.
        h.chunkSumRpj[chunkIdx] += h.evtSumRpj;
        h.chunkCntRpj[chunkIdx] += h.evtCntRpj;

        h.evtSumRpj = 0.;
        h.evtCntRpj = 0; // This variable already knows how many Lambdas there are for each cut, be it in eta, DCA or pT min.
                         // In that regard, it could possibly be used as a weight so as to not consider events where 1 Lambda
                         // passed in the cuts as statistically equivalent to events where all 4 Lambdas passed the cuts.
    }
    h.evtIdxRpj++; // Always advance: keeps chunk index consistent across all scenarios
}


// ==========================================================================
/**
 * @brief Flushes the event accumulators for all twelve ScenarioHistos inside
 *        a FamilyHistos struct.
 *
 * @details
 * Calls FlushEventMean() on every (cut scenario) x (eta half) combination
 * held by @p f.  This covers the full 4 x 3 = 12 ScenarioHistos members:
 *
 *   NoCuts     x { EtaPos, EtaNeg, All }
 *   pTCutOnly  x { EtaPos, EtaNeg, All }
 *   DCACutOnly x { EtaPos, EtaNeg, All }
 *   BothCuts   x { EtaPos, EtaNeg, All }
 *
 * Must be called:
 *   (a) inside the jet-reshuffle block, BEFORE drawing the new jet direction,
 *       to close the event that just ended; and
 *   (b) once after the main loop, to flush the last partial event that may
 *       not have triggered a reshuffle.
 *
 * @param f  Reference to the FamilyHistos whose twelve accumulators are flushed.
 */
// ==========================================================================
static void FlushEventMeansFamily(FamilyHistos& f) {
    FlushEventMean(f.NC_Pos); FlushEventMean(f.NC_Neg); FlushEventMean(f.NC_All);
    FlushEventMean(f.PT_Pos); FlushEventMean(f.PT_Neg); FlushEventMean(f.PT_All);
    FlushEventMean(f.DC_Pos); FlushEventMean(f.DC_Neg); FlushEventMean(f.DC_All);
    FlushEventMean(f.BC_Pos); FlushEventMean(f.BC_Neg); FlushEventMean(f.BC_All);
}

// ==========================================================================
/**
 * @brief Computes the per-chunk means for one ScenarioHistos and fills
 *        hChunkMeansRingProxyJet with one entry per non-empty chunk.
 *
 * @details
 * Must be called ONCE, after the last FlushEventMean() call for this
 * ScenarioHistos (i.e. after both FlushEventMeansFamily calls that follow
 * the main loop).
 *
 * For each chunk k with at least one accumulated entry:
 *   mu_k = chunkSumRpj[k] / chunkCntRpj[k]
 * is filled into hChunkMeansRingProxyJet (unweighted).
 *
 * After all K_filled non-empty chunks are filled:
 *   hChunkMeansRingProxyJet->GetMean()      = grand mean  (agrees with TProfile)
 *   hChunkMeansRingProxyJet->GetMeanError() = stddev(mu_k) / sqrt(K_filled)
 *                                           = the bulletproof final error
 *
 * @param h  Reference to the ScenarioHistos to finalise.
 */
// ==========================================================================
static void FinalizeChunks(ScenarioHistos& h) {
    for (int k = 0; k < kChunks; ++k) {
        if (h.chunkCntRpj[k] > 0) {
            double mu_k = h.chunkSumRpj[k] / (double)h.chunkCntRpj[k];
            h.hChunkMeansRingProxyJet->Fill(mu_k); // unweighted: ROOT SEM = stddev/sqrt(K_filled)
            h.pChunkMeansRingProxyJetIntegrated->Fill(0., mu_k);
        }
    }
}


// ==========================================================================
/**
 * @brief Calls FinalizeChunks() on all twelve ScenarioHistos in a FamilyHistos.
 *
 * @details
 * Must be called once, after the two final FlushEventMeansFamily() calls
 * that drain the last partial event at the end of the main loop.
 *
 * @param f  Reference to the FamilyHistos to finalise.
 */
// ==========================================================================
static void FinalizeChunksFamily(FamilyHistos& f) {
    FinalizeChunks(f.NC_Pos); FinalizeChunks(f.NC_Neg); FinalizeChunks(f.NC_All);
    FinalizeChunks(f.PT_Pos); FinalizeChunks(f.PT_Neg); FinalizeChunks(f.PT_All);
    FinalizeChunks(f.DC_Pos); FinalizeChunks(f.DC_Neg); FinalizeChunks(f.DC_All);
    FinalizeChunks(f.BC_Pos); FinalizeChunks(f.BC_Neg); FinalizeChunks(f.BC_All);
}

// ==========================================================================
/**
 * @brief Flushes the Kahan accumulators into their corresponding TH1D targets
 * for a single ScenarioHistos.
 *
 * @param h  Reference to the ScenarioHistos containing the accumulators and TH1Ds.
 */
// ==========================================================================
static void FlushKahan(ScenarioHistos& h) {
    h.kRingProxy.FlushToTH1(h.hRingProxy_Kahan);
    h.kRingProxyJet.FlushToTH1(h.hRingProxyJet_Kahan);
    h.kRingProxyJet_JetEtaPos.FlushToTH1(h.hRingProxyJet_JetEtaPos_Kahan);
    h.kRingProxyJet_JetEtaNeg.FlushToTH1(h.hRingProxyJet_JetEtaNeg_Kahan);
}

// ==========================================================================
/**
 * @brief Calls FlushKahan() on all twelve ScenarioHistos in a FamilyHistos.
 *
 * @details
 * Must be called once, after the main event loop is completely finished,
 * to lock the high-precision sums into the TH1D objects before writing 
 * to the output ROOT file.
 *
 * @param f  Reference to the FamilyHistos to finalize.
 */
// ==========================================================================
static void FlushKahanFamily(FamilyHistos& f) {
    FlushKahan(f.NC_Pos); FlushKahan(f.NC_Neg); FlushKahan(f.NC_All);
    FlushKahan(f.PT_Pos); FlushKahan(f.PT_Neg); FlushKahan(f.PT_All);
    FlushKahan(f.DC_Pos); FlushKahan(f.DC_Neg); FlushKahan(f.DC_All);
    FlushKahan(f.BC_Pos); FlushKahan(f.BC_Neg); FlushKahan(f.BC_All);
}


// ==========================================================================
/**
 * @brief Generates unpolarized Lambda decays and studies acceptance-induced
 *        fake-polarization effects by varying kinematic and topological cuts.
 *
 * @details
 * Histograms are stored in two parallel families:
 *   - WithoutEtaGate: no requirement on daughter pseudorapidity.  Kept as a
 *     reference to study the inconsistency that arises when Lambdas
 *     whose daughters fall outside the detector acceptance are included when
 *     the Lambdas themselves were only generated within acceptance (approximately!
 *     The Lambdas are generated in a rapidity cut, not pseudorapidity, but
 *     still, not applying a pseudorap cut would be inconsistent if they are
       generated in a small region but may be detected in a large eta region).
 *   - WithEtaGate: additionally requires |eta_daughter| < etaMaxDetector for
 *     BOTH daughters before any other cut is applied.  This is the physically
 *     consistent set and should be used for all physics conclusions.
 *
 * Within each family, four cut scenarios are filled: NoCuts, pTCutOnly,
 * DCACutOnly, and BothCuts.  Within each scenario, histograms are split into
 * EtaPos (eta_Lambda >= 0), EtaNeg (eta_Lambda < 0), and All sub-directories
 * (see file header for the full output structure).
 *
 * @param nLambdas       Number of Lambdas to simulate (default: 10,000,000).
 * @param outputPath     Output ROOT file path (default: helicityEffOutput.root).
 * @param Bz_Tesla       Magnetic field along z [T], positive = ALICE default
 *                       (default: +0.5 T).
 * @param pTmin_Lambda   Minimum Lambda pT for generation [GeV/c] (default: 0.0).
 * @param pTmax_Lambda   Maximum Lambda pT for generation [GeV/c] (default: 10.0).
 * @param rapMax_Lambda  Maximum |rapidity| for GENERATING Lambdas (default: 5.0).
 *                       Set wide enough to capture all Lambdas whose daughters
 *                       can fall inside |eta| < etaMaxDetector.
 * @param etaMaxDetector Maximum |eta| for daughter tracks (default: 0.9,
 *                       matching the ALICE inner barrel). Applied as a gate
 *                       to select the WithEtaGate family; not applied to the
 *                       Lambda itself. Detectors only see charged particles,
 *                       and only within a select acceptance in |eta|.
 * @param T_thermal      Boltzmann temperature for the mT spectrum [GeV] (default: 0.3).
 * @param pTmin_proton   Minimum proton pT cut [GeV/c] (default: 0.0).
 * @param pTmin_pion     Minimum pion pT cut [GeV/c] (default: 0.0).
 * @param dcaMin_proton  Minimum proton DCA_xy to the primary vertex [cm] (default: 0.0).
 * @param dcaMin_pion    Minimum pion DCA_xy to the primary vertex [cm] (default: 0.0).
 * @param seed           TRandom3 seed; 0 = time-based random seed (default: 0).
 */
// ==========================================================================
void helicityEfficiencyToyModel(
    long        nLambdas       = 10000000,
    const char* outputPath     = "helicityEffOutput.root",
    double      Bz_Tesla       = 0.5,
    double      pTmin_Lambda   = 0.0,
    double      pTmax_Lambda   = 10.0,
    double      rapMax_Lambda  = 5.0,
    double      etaMaxDetector = 0.9,
    double      T_thermal      = 0.300,
    double      pTmin_proton   = 0.0,
    double      pTmin_pion     = 0.0,
    double      dcaMin_proton  = 0.00,
    double      dcaMin_pion    = 0.00,
    int         seed           = 0)
{
    // -----------------------------------------------------------------------
    // 0) Print run parameters to stdout for a clear record in the log
    // -----------------------------------------------------------------------
    printf("\n");
    printf("========================================================\n");
    printf("  helicityEfficiencyToyModel -- run parameters\n");
    printf("========================================================\n");
    printf("  N Lambdas to simulate : %ld\n",  nLambdas);
    printf("  Output ROOT file      : %s\n",   outputPath);
    printf("  Magnetic field Bz     : %.3f T\n",  Bz_Tesla);
    printf("  Lambda pT range       : [%.2f, %.2f] GeV/c\n", pTmin_Lambda, pTmax_Lambda);
    printf("  Lambda |rapidity| max : %.2f\n", rapMax_Lambda);
    printf("  Detector |eta| max    : %.2f (daughter eta gate)\n", etaMaxDetector);
    printf("  Thermal temperature T : %.3f GeV\n", T_thermal);
    printf("  Min proton pT cut     : %.3f GeV/c\n", pTmin_proton);
    printf("  Min pion pT cut       : %.3f GeV/c\n", pTmin_pion);
    printf("  Min proton DCA_xy     : %.4f cm\n", dcaMin_proton);
    printf("  Min pion DCA_xy       : %.4f cm\n", dcaMin_pion);
    printf("  RNG seed              : %d\n",   seed);
    printf("========================================================\n\n");


    // -----------------------------------------------------------------------
    // 1) Open output ROOT file and set up directory structure
    // -----------------------------------------------------------------------
    TFile* outFile = TFile::Open(outputPath, "RECREATE");
    if (!outFile || outFile->IsZombie()) {
        printf("ERROR: Cannot open output file '%s'. Aborting.\n", outputPath);
        return;
    }

    // Two top-level directories, one per histogram family.
    // WithoutEtaGate: no daughter eta requirement (kept for comparison only;
    //                 see file header for the physics motivation).
    // WithEtaGate:    requires |eta_daughter| < etaMaxDetector (consistent set).
    TDirectory* dirNG = outFile->mkdir("WithoutEtaGate");
    TDirectory* dirEG = outFile->mkdir("WithEtaGate");

    // Create the four cut-scenario sub-directories in each family and book
    // all histograms.  CreateFamily calls CreateSubdirs four times internally.
    FamilyHistos famNG; // Without eta gate
    FamilyHistos famEG; // With eta gate
    CreateFamily(dirNG, famNG, etaMaxDetector);
    CreateFamily(dirEG, famEG, etaMaxDetector);

    // Extra kinematics directory: Lambda distribution BEFORE any decay cuts,
    // only subject to the generation acceptance (pT and rapidity windows).
    TDirectory* dirKin = outFile->mkdir("Kinematics");
    dirKin->cd();
    TH1D* hKin_pT_lambda  = new TH1D("hKin_pT_lambda", "Generated #Lambda p_{T};p_{T}^{#Lambda} [GeV/c];Counts", 100, 0., 10.);
    TH1D* hKin_rap_lambda = new TH1D("hKin_rap_lambda", "Generated #Lambda rapidity;y_{#Lambda};Counts", 36, -rapMax_Lambda, rapMax_Lambda);
    TH1D* hKin_eta_lambda = new TH1D("hKin_eta_lambda", "Generated #Lambda pseudorapidity;#eta_{#Lambda};Counts", 36, -etaMaxDetector, etaMaxDetector);
    TH1D* hKin_phi_lambda = new TH1D("hKin_phi_lambda", "Generated #Lambda azimuth;#phi_{#Lambda} [rad];Counts", 64, -Pi, Pi);
    TH1D* hKin_decayR     = new TH1D("hKin_decayR", "All transverse decay radii;r_{decay} [cm];Counts", 150, 0., 150.);
    TH1D* hKin_pT_proton  = new TH1D("hKin_pT_proton", "All proton p_{T} (pre-cut);p_{T}^{p} [GeV/c];Counts", 100, 0., 5.);
    TH1D* hKin_pT_pion    = new TH1D("hKin_pT_pion", "All pion p_{T} (pre-cut);p_{T}^{#pi} [GeV/c];Counts", 100, 0., 5.);
    TH1D* hKin_DCA_proton = new TH1D("hKin_DCA_proton", "All proton DCA_{xy} (pre-cut);DCA_{xy}^{p} [cm];Counts", 200, 0., 10.);
    TH1D* hKin_DCA_pion   = new TH1D("hKin_DCA_pion", "All pion DCA_{xy} (pre-cut);DCA_{xy}^{#pi} [cm];Counts", 200, 0., 10.);
    TH1D* hKin_eta_jet    = new TH1D("hKin_eta_jet", "Jet pseudorapidity; #eta_{Jet};Counts", 36, -etaMaxDetector, etaMaxDetector);
    outFile->cd(); // Back to root of output file

    // -----------------------------------------------------------------------
    // 2) Initialise the random number generator and phase-space decay engine
    // -----------------------------------------------------------------------
    TRandom3 rng(seed);

    // TGenPhaseSpace: decays a parent 4-vector into N daughters according to
    // Lorentz-invariant phase space (flat in LIPS = unpolarized).
    TGenPhaseSpace decayGen; // Its gRandom should be zero-initialized so uses a random seed.
                             // This means the "seed" parameter does not exactly fix all outputs,
                             // but it fixes jet direction sampling.

    // Daughter masses array for Lambda -> p + pi-
    Double_t daughterMasses[2] = { kMassProton, kMassPion };

    // -----------------------------------------------------------------------
    // 3) Main event loop
    // -----------------------------------------------------------------------
    printf("Starting generation of %ld Lambdas...\n", nLambdas);

    long nGenerated = 0; // Lambdas for which decay was successfully set up
    long nCollinear = 0; // Lambdas skipped because pT/|p| < kMinSinTheta
    int  nLambdasSinceJetShuffle = 0; // Lambdas since the jet direction last shuffled

    // Progress report every 10%
    long nProgress = nLambdas / 10;
    if (nProgress < 1) nProgress = 1;

    // Calculating useful variables for the SampleLambdaPt function:
    // We sample a Boltzmann function, which has an analytic maximum that can be
    // calculated exactly for the rejection-sampling envelope below.
    //
    // For best rejection-sampling efficiency, we must evaluate the maximum
    // INSIDE the sampled interval [pTmin_Lambda, pTmax_Lambda] to avoid 
    // unnecessary rejection inefficiency, keeping all checks outside the hot event loop.
    double ptPeak = ComputeThermalPtMaximum(T_thermal);

    if (ptPeak < pTmin_Lambda) {
        ptPeak = pTmin_Lambda;
        printf("SampleLambdaPt: thermal spectrum peak lies below sampled pT range; using pTmin_Lambda for envelope maximum.\n");
    } else if (ptPeak > pTmax_Lambda) {
        ptPeak = pTmax_Lambda;
        printf("SampleLambdaPt: thermal spectrum peak lies above sampled pT range; using pTmax_Lambda for envelope maximum.\n");
    }

    double mTpeak = std::sqrt(ptPeak * ptPeak + kMassLambda * kMassLambda);
    double fmax = 1.10 * ptPeak * std::exp(-mTpeak / T_thermal);

    // Creating a jet proxy facing a random direction
        // Preliminary analysis on data shows that AEE effects indeed persist even with a random jet direction.
        // This makes sense, as the AEE is an effect somewhat "intrinsic" to the Lambda daughters interacting
        // with the magnetic field, not so much of the ring observable's proxy direction
    double randPhiJet = rng.Uniform(0., TwoPi);
    double sinPhiJet  = std::sin(randPhiJet);
    double cosPhiJet  = std::cos(randPhiJet);
    // Sampling cos(theta) only within the detector range, also for consistency:
    double cosThetaMax = std::tanh(etaMaxDetector);
    double cosThetaJet = rng.Uniform(-cosThetaMax, cosThetaMax); // Uniformly sampled in solid angle, yet within the detectable region
    double sinThetaJet = std::sqrt(1. - cosThetaJet * cosThetaJet);
    
    // Calculating the jet component directions (could use TVector3, but this formulation is a bit faster in the hot loop)
    // TVector3 jetProxyDirection(sinThetaJet * cosPhiJet, sinThetaJet * sinPhiJet, cosThetaJet);
    double jx = sinThetaJet * cosPhiJet;
    double jy = sinThetaJet * sinPhiJet;
    double jz = cosThetaJet;
    double eta_jet = std::atanh(cosThetaJet);
    hKin_eta_jet->Fill(eta_jet); // Colecting the kinematics of the first sampled jet as well

    for (long iLam = 0; iLam < nLambdas; ++iLam) {
        // ---- Progress printout ----
        if ((iLam + 1) % nProgress == 0)
            std::cout << "  " << (iLam + 1) << " / " << nLambdas << " (" << 100. * (iLam + 1) / nLambdas << "%)" << std::endl;

        // ==================================================================
        // 3.1  Generate Lambda 4-momentum
        // ==================================================================
        // Sample pT from thermal mT spectrum [GeV/c]
        double pT_lam = SampleLambdaPt(rng, T_thermal, pTmin_Lambda, pTmax_Lambda, fmax);

        // Sample Lambda azimuthal angle uniformly in [0, 2pi)
        double phi_lam = rng.Uniform(0., TwoPi);

        // Sample Lambda rapidity uniformly in [-rapMax, +rapMax]
        // This makes the distribution NON UNIFORM in cos(theta) (i.e., a non-spherical distribution!)
        // The idea is to have something Lorentz-invariant, like what we see in HI collisions.
        double rap_lam = rng.Uniform(-rapMax_Lambda, rapMax_Lambda);

        // Construct the 4-momentum from (pT, phi, y):
        //   mT = sqrt(pT^2 + m^2)
        //   E  = mT * cosh(y)
        //   pz = mT * sinh(y)
        //   px = pT * cos(phi)
        //   py = pT * sin(phi)
        double mT_lam  = std::sqrt(pT_lam * pT_lam + kMassLambda * kMassLambda);
        double E_lam   = mT_lam * std::cosh(rap_lam);
        double pz_lam  = mT_lam * std::sinh(rap_lam);
        double px_lam  = pT_lam * std::cos(phi_lam);
        double py_lam  = pT_lam * std::sin(phi_lam);
        double p_lam = std::sqrt(pT_lam*pT_lam + pz_lam*pz_lam);

        // Lambda pseudorapidity (for eta-split histograms)
        double eta_lam = std::asinh(pz_lam / pT_lam);

        // Build TLorentzVector for Lambda (used by TGenPhaseSpace, so can't really swap this for Vector4D...)
        TLorentzVector lv_lam(px_lam, py_lam, pz_lam, E_lam);
            // Converting the Lambda 4-vec to the latest structure (should yield faster computation still!):
        PxPyPzEVector lambda4vec(px_lam, py_lam, pz_lam, E_lam);

        // Fill pre-cut kinematics
        hKin_pT_lambda->Fill(pT_lam);
        hKin_rap_lambda->Fill(rap_lam);
        hKin_eta_lambda->Fill(eta_lam);
        hKin_phi_lambda->Fill(wrapToPiFast(phi_lam));

        // ==================================================================
        // 3.2  Build Lambda frame axes
        // ==================================================================
        // e1 = Lambda unit momentum vector (forward-backward axis)
        XYZVector e1(px_lam / p_lam, py_lam / p_lam, pz_lam / p_lam);

        // e2 = (p_Lambda x z_hat) / |p_Lambda x z_hat| (transverse, left-right axis)
            // Explicitly:  e2 = (py, -px, 0) / pT
        // Reject nearly beam-collinear Lambdas (pT/|p| is too small):
        if (pT_lam / p_lam < kMinSinTheta) {
            ++nCollinear;
            continue;
        }
        XYZVector e2(py_lam / pT_lam, -px_lam / pT_lam, 0.);

        // e3 = e1 x e2  (in-out axis)
        // Explicitly:  e3 = (pz*px, pz*py, -pT^2) / (|p|*pT)
        XYZVector e3 = e1.Cross(e2);  // Unit by construction since e1, e2 are unit orthonormal

        // ==================================================================
        // 3.3  Sample proper decay length and compute decay vertex
        // ==================================================================
        // Sample proper decay length from exponential distribution:
        //   l_proper = -ctau * ln(U),  U ~ Uniform(0,1)
        // Then the lab-frame decay length is:
        //   r_lab = (p_vec / m_Lambda) * l_proper
        // since:
        //   l_lab = beta*gamma*c*tau = (|p|/m_Lambda) * l_proper
        // and:
        //   r_hat = p_vec / |p|
        // giving:
        //   r_lab = r_hat * l_lab = (p_vec/m_Lambda) * l_proper
        // Decay vertex position:
        //   v = (px, py, pz) / |p| * l_lab
        double l_proper = -kCTauLambda * std::log(rng.Rndm()); // Sampling exponentially decaying distribution

        // -- Lab-frame decay vertex coordinates [cm] --
        double xv = (px_lam / kMassLambda) * l_proper;
        double yv = (py_lam / kMassLambda) * l_proper;
        double zv = (pz_lam / kMassLambda) * l_proper;

        // -- Transverse decay radius [cm] --
        double decayR = std::sqrt(xv * xv + yv * yv);
        hKin_decayR->Fill(decayR);

        // ==================================================================
        // 3.4  Generate unpolarized Lambda decay via TGenPhaseSpace
        // ==================================================================
        // TGenPhaseSpace::SetDecay() sets up the kinematics.
        // TGenPhaseSpace::Generate() produces one random phase-space point.
        // Daughters are indexed as: 0 = proton, 1 = pion.
        decayGen.SetDecay(lv_lam, 2, daughterMasses);
        decayGen.Generate();  // Isotropic, unpolarized decay

        // Retrieve daughter 4-vectors in the LAB frame
        // (TODO: could calculate this by hand, as done in PythiaGenMin's MC by the way! Much more efficient than TGenPhaseSpace, I guess...)
        // (Just make sure to generate homogeneously in cosTheta instead of Theta, for proper 3D isotropic generation of protons and pions!)
        TLorentzVector tlv_proton = *(decayGen.GetDecay(0));
        TLorentzVector tlv_pion   = *(decayGen.GetDecay(1));

        // Lab-frame transverse momenta of daughters [GeV/c]
        double px_p = tlv_proton.Px();
        double py_p = tlv_proton.Py();
        double pz_p = tlv_proton.Pz();
        double E_p  = tlv_proton.E();

        double px_pi = tlv_pion.Px();
        double py_pi = tlv_pion.Py();
        double pz_pi = tlv_pion.Pz();
        double E_pi  = tlv_pion.E();

        // Now convert to the most recent ROOT toolset, faster and optimization-friendly:
        PxPyPzEVector proton4vec(px_p, py_p, pz_p, E_p);
        PxPyPzEVector pion4vec(px_pi, py_pi, pz_pi, E_pi);

        double pT_p  = proton4vec.Pt();
        double pT_pi = pion4vec.Pt();

        // Fill pre-cut daughter kinematics
        hKin_pT_proton->Fill(pT_p);
        hKin_pT_pion->Fill(pT_pi);

        // -- Daughter pseudorapidities (used for the eta gate below) --
        // Lambdas are generated over a wide rapidity window (rapMax_Lambda = 5.0
        // by default) so that ALL Lambdas whose daughters fall inside the inner
        // barrel acceptance (|eta| < etaMaxDetector = 0.9) should have been included.
        // The WithEtaGate family requires BOTH daughters to satisfy this cut,
        // making the acceptance consistent with that of a real detector.
        // The WithoutEtaGate family does NOT apply this requirement and is
        // kept only to demonstrate the inconsistency of omitting it.
        double eta_p  = proton4vec.Eta();
        double eta_pi = pion4vec.Eta();

        // ==================================================================
        // 3.5  Boost proton to Lambda rest frame
        // ==================================================================
        // Boost the proton 4-vector from the lab frame into the Lambda rest frame.
        // In TLorentzVector, BoostVector() returns beta = p/E (rest --> lab).
        // To go from lab --> rest frame, we must use the opposite boost: -p/E.
        // TLorentzVector lv_proton_star = lv_proton;
        // lv_proton_star.Boost(-lv_lam.BoostVector());

        // New boost style:
        auto proton4vec_star = ROOT::Math::VectorUtil::boost(proton4vec, lambda4vec.BoostToCM()); // BoostToCM gives a trivector that goes from laboratory
                                                                                                  // frame to Lambda's rest frame (convenient new function,
                                                                                                  // different from TLorentzVector's BoostVector())
                                                                                                  // This already contains the "minus sign" for the boost
        // Unit vector of proton momentum in Lambda rest frame
        XYZVector p_star_unit = proton4vec_star.Vect().Unit();

        // ==================================================================
        // 3.6  Compute the three key angular observables
        // ==================================================================
        // -- cos(theta*): proton emission vs Lambda momentum direction (e1) --
        // This is the forward-backward angle.
        double cosTheta = p_star_unit.Dot(e1);

        // -- phi*: proton azimuth around the Lambda axis --
        // phi* = 0 is in the "left" direction (e2 = p_Lambda x z_hat).
        // phi* = pi/2 is in the "in-out" direction (e3 = e1 x e2).
        double phi_star = std::atan2(p_star_unit.Dot(e3), p_star_unit.Dot(e2));

        // ==================================================================
        // 3.7  Compute the ring observable proxy (z_hat as if jet direction)
        // ==================================================================
        // t_hat = z_hat = (0, 0, 1)
        // t_hat x lambda_momentum = z_hat x (lx, ly, lz) = (-ly, lx, 0)
        // |z_hat x lambda_momentum| = sqrt(lx^2 + ly^2) = pT / |p| = sinTheta_lam
        //
        // Normalized cross product = (-ly, lx, 0) / sinTheta_lam
        //                          = (-py, px, 0) / pT (after scaling)
        // ring_proxy = (3/alpha) * p_star_unit . (-py, px, 0) / pT
        double ringProxy = kPolPrefactor * (-p_star_unit.X() * py_lam + p_star_unit.Y() * px_lam) / pT_lam;

        // ==================================================================
        // 3.7.1  Compute the ring observable proxy using randomly sampled jet
        // ==================================================================
        if (nLambdasSinceJetShuffle == 3){
            // Close the current event before drawing a new jet:
            FlushEventMeansFamily(famNG);
            FlushEventMeansFamily(famEG);

            // Every fourth Lambda that used the same jet direction will trigger a reshuffle of the jet direction
            // Physically, this would mean saying that there are about 4 Lambdas in the same event. One can test
            // if this makes much of an impact due to introducing a correlation between the Lambdas.
            // Realistically, data shows about 2 Lambda *candidates* on average, for events that do have a Lambda.
            // And data also shows about 2 jets (on average) per event, in events that do have a jet, so this may
            // not be so far off. Probably would be interesting to know the correlation between N_Lambda and N_Jets
            // per event, to estimate this number a little bit better.
            randPhiJet = rng.Uniform(0., TwoPi);
            cosThetaJet = rng.Uniform(-cosThetaMax, cosThetaMax);
            sinPhiJet = std::sin(randPhiJet);
            cosPhiJet = std::cos(randPhiJet);
            sinThetaJet = std::sqrt(1. - cosThetaJet * cosThetaJet);
            
            // Re-evaluate jet direction components:
            jx = sinThetaJet * cosPhiJet;
            jy = sinThetaJet * sinPhiJet;
            jz = cosThetaJet;
            eta_jet = std::atanh(cosThetaJet);
            hKin_eta_jet->Fill(eta_jet); // This is filled inside the if block, as we only want to map this once for each jet

            nLambdasSinceJetShuffle = 0;
        }
        else nLambdasSinceJetShuffle++;

        // Calculating the cross-product by hand (just a tiny bit faster, probably):
        // 1) Cross product c = jet_hat x p_Lambda
        double cx = jy * pz_lam - jz * py_lam;
        double cy = jz * px_lam - jx * pz_lam;
        double cz = jx * py_lam - jy * px_lam;

        // 2) Normalization factor for the cross product
        double inv_c_norm = 1.0 / std::sqrt(cx * cx + cy * cy + cz * cz);

        // 3) Dot product with p_star_unit
        double ringProxyJet = kPolPrefactor * (p_star_unit.X() * cx + p_star_unit.Y() * cy + p_star_unit.Z() * cz) * inv_c_norm;

        // ==================================================================
        // 3.8  Compute DCA_xy of each daughter to the primary vertex
        // ==================================================================
        // The daughter tracks originate at the decay vertex (xv, yv).
        // We compute the analytical transverse DCA using the helix formula.
        double dca_proton = ComputeDCAxy(xv, yv, px_p, py_p, pT_p, kChargeProton, Bz_Tesla);
        double dca_pion   = ComputeDCAxy(xv, yv, px_pi, py_pi, pT_pi, kChargePion, Bz_Tesla);

        // Fill pre-cut DCA histograms
        hKin_DCA_proton->Fill(dca_proton);
        hKin_DCA_pion->Fill(dca_pion);

        // ==================================================================
        // 3.9  Apply selection cuts
        // ==================================================================
        // -- Eta gate: both daughters must lie inside the detector acceptance --
        // This is the acceptance pre-condition for the WithEtaGate family.
        bool passEtaGate = (std::fabs(eta_p)  < etaMaxDetector) && (std::fabs(eta_pi) < etaMaxDetector);

        // -- pT cut: both daughters must exceed minimum pT --
        bool passPtCut = (pT_p  >= pTmin_proton) && (pT_pi >= pTmin_pion);

        // -- DCA cut: both daughters must be sufficiently displaced from PV --
        // (This selects genuine secondary tracks in data, i.e. rejects primaries.
        //  We lack actual primary-secondary distinction here, so it just emulates it)
        bool passDcaCut = (dca_proton >= dcaMin_proton) && (dca_pion >= dcaMin_pion);

        // -- Eta > 0 flag (routes fill into EtaPos vs EtaNeg subdirectory) --
        bool etaPos = (eta_lam >= 0.);

        // ==================================================================
        // 3.10  Fill histogram families
        // ==================================================================
        // WithoutEtaGate: fill for all events regardless of daughter eta.
        // Kept as a comparison to show the inconsistency of omitting the gate.
        FillFamily(famNG, passPtCut, passDcaCut, etaPos,
                   cosTheta, phi_star, ringProxy, ringProxyJet, decayR,
                   pT_p, pT_pi, dca_proton, dca_pion, pT_lam, eta_lam, eta_jet,
                   p_star_unit.X(), p_star_unit.Y(), p_star_unit.Z(), px_lam, py_lam, phi_lam, pz_lam);

        // WithEtaGate: only fill when both daughters are inside the acceptance.
        // This is the physically consistent set.
        if (passEtaGate)
            FillFamily(famEG, passPtCut, passDcaCut, etaPos,
                       cosTheta, phi_star, ringProxy, ringProxyJet, decayR,
                       pT_p, pT_pi, dca_proton, dca_pion, pT_lam, eta_lam, eta_jet,
                       p_star_unit.X(), p_star_unit.Y(), p_star_unit.Z(), px_lam, py_lam, phi_lam, pz_lam);

        ++nGenerated;

    } // end main loop over Lambdas

    // Flush whatever partial event is left in the accumulators:
    FlushEventMeansFamily(famNG);
    FlushEventMeansFamily(famEG);

    // Compute chunk means now that all events have been flushed:
    FinalizeChunksFamily(famNG);
    FinalizeChunksFamily(famEG);

    // Flush Kahan accumulators to their TH1D containers:
    FlushKahanFamily(famNG);
    FlushKahanFamily(famEG);

    // -----------------------------------------------------------------------
    // 4) Print summary statistics
    // -----------------------------------------------------------------------
    printf("\n");
    printf("========================================================\n");
    printf("  Generation complete.\n");
    printf("  Total Lambdas requested : %ld\n",   nLambdas);
    printf("  Successfully generated  : %ld\n",   nGenerated);
    printf("  Skipped (collinear)     : %ld\n",   nCollinear);
    printf("========================================================\n\n");

    // Helper lambda to print cut statistics for one family
    auto PrintCutStats = [](const char* familyLabel, const FamilyHistos& f) {
        double nNC = (double)f.NC_All.h1d_pT_lambda->GetEntries();
        double nPT = (double)f.PT_All.h1d_pT_lambda->GetEntries();
        double nDC = (double)f.DC_All.h1d_pT_lambda->GetEntries();
        double nBC = (double)f.BC_All.h1d_pT_lambda->GetEntries();
        printf("  [%s] cut statistics (fraction of NoCuts in this family):\n", familyLabel);
        printf("    NoCuts    : %.0f  (100%%)\n",       nNC);
        printf("    pTCutOnly : %.0f  (%.1f%%)\n",     nPT, 100.*nPT/nNC);
        printf("    DCACutOnly: %.0f  (%.1f%%)\n",     nDC, 100.*nDC/nNC);
        printf("    BothCuts  : %.0f  (%.1f%%)\n",     nBC, 100.*nBC/nNC);
        printf("\n");
    };

    PrintCutStats("WithoutEtaGate", famNG);
    PrintCutStats("WithEtaGate",    famEG);

    // Helper lambda to print the jet-proxy ring observable for one family,
    // showing both the naive TProfile SEM and the event-mean-based uncertainty
    // side by side so they can be compared directly.
    auto PrintRingJetFamily = [](const char* familyLabel, const FamilyHistos& f) {

        // Prints one scenario row: TProfile SEM vs event-mean sigma/sqrt(N_events)
        auto PrintRingJet = [](const char* label,
                               const ScenarioHistos& hA,
                               const ScenarioHistos& hP,
                               const ScenarioHistos& hN) {

            // --- TProfile SEM (assumes independent entries) ---
            double rAll_prof = hA.pRingProxyJet->GetBinContent(1);
            double eAll_prof = hA.pRingProxyJet->GetBinError(1);
            double rPos_prof = hP.pRingProxyJet->GetBinContent(1);
            double ePos_prof = hP.pRingProxyJet->GetBinError(1);
            double rNeg_prof = hN.pRingProxyJet->GetBinContent(1);
            double eNeg_prof = hN.pRingProxyJet->GetBinError(1);

            // --- Kahan TH1D-flushed values (much more precise than simple TProfiles!) ---
                // Could do this directly with the KahanAccumulator struct logic, but chose to standardize with the TH1 logic from above:
            // double rAll_Kah = hA.kRingProxyJet.GetMean();
            // double eAll_Kah = hA.kRingProxyJet.GetSEM();
            double rAll_Kah  = hA.hRingProxyJet_Kahan->GetBinContent(1);
            double eAll_Kah  = hA.hRingProxyJet_Kahan->GetBinError(1);
            double rPos_Kah  = hP.hRingProxyJet_Kahan->GetBinContent(1);
            double ePos_Kah  = hP.hRingProxyJet_Kahan->GetBinError(1);
            double rNeg_Kah  = hN.hRingProxyJet_Kahan->GetBinContent(1);
            double eNeg_Kah  = hN.hRingProxyJet_Kahan->GetBinError(1);

            // --- Event-mean estimator: sigma(R_e) / sqrt(N_events) ---
            double rAll_evt = hA.hEventMeanRingProxyJet->GetMean();
            double eAll_evt = hA.hEventMeanRingProxyJet->GetMeanError();
            double rPos_evt = hP.hEventMeanRingProxyJet->GetMean();
            double ePos_evt = hP.hEventMeanRingProxyJet->GetMeanError();
            double rNeg_evt = hN.hEventMeanRingProxyJet->GetMean();
            double eNeg_evt = hN.hEventMeanRingProxyJet->GetMeanError();

            // Number of events (jet groups) that contributed at least one entry
            double nEvtAll = (double)hA.hEventMeanRingProxyJet->GetEntries();
            double nEvtPos = (double)hP.hEventMeanRingProxyJet->GetEntries();
            double nEvtNeg = (double)hN.hEventMeanRingProxyJet->GetEntries();
            
            // --- Chunking estimator (stddev(mu_k) / sqrt(K_filled), unweighted) ---
            double rAll_chk = hA.hChunkMeansRingProxyJet->GetMean();
            double eAll_chk = hA.hChunkMeansRingProxyJet->GetMeanError();
            double rPos_chk = hP.hChunkMeansRingProxyJet->GetMean();
            double ePos_chk = hP.hChunkMeansRingProxyJet->GetMeanError();
            double rNeg_chk = hN.hChunkMeansRingProxyJet->GetMean();
            double eNeg_chk = hN.hChunkMeansRingProxyJet->GetMeanError();

            double nChkAll = (double)hA.hChunkMeansRingProxyJet->GetEntries();
            double nChkPos = (double)hP.hChunkMeansRingProxyJet->GetEntries();
            double nChkNeg = (double)hN.hChunkMeansRingProxyJet->GetEntries();

            printf("    %-12s\n", label);
            printf("      TProfile SEM    All: %+.4e +/- %.4e  EtaPos: %+.4e +/- %.4e  EtaNeg: %+.4e +/- %.4e\n",
                   rAll_prof, eAll_prof, rPos_prof, ePos_prof, rNeg_prof, eNeg_prof);
            printf("      Kahan SEM    All: %+.4e +/- %.4e  EtaPos: %+.4e +/- %.4e  EtaNeg: %+.4e +/- %.4e\n",
                   rAll_Kah, eAll_Kah, rPos_Kah, ePos_Kah, rNeg_Kah, eNeg_Kah);
            printf("      Event-mean      All: %+.4e +/- %.4e  EtaPos: %+.4e +/- %.4e  EtaNeg: %+.4e +/- %.4e\n",
                   rAll_evt,  eAll_evt,  rPos_evt,  ePos_evt,  rNeg_evt,  eNeg_evt);
            printf("      Chunking (%2.0f chk) All: %+.4e +/- %.4e  EtaPos: %+.4e +/- %.4e  EtaNeg: %+.4e +/- %.4e\n",
                   nChkAll, rAll_chk, eAll_chk, rPos_chk, ePos_chk, rNeg_chk, eNeg_chk);
            printf("      N_events   All: %.0f (evt) / %.0f (chk)   EtaPos: %.0f / %.0f   EtaNeg: %.0f / %.0f\n",
                   nEvtAll, nChkAll, nEvtPos, nChkPos, nEvtNeg, nChkNeg);

            double ratioAll_ec = (eAll_prof > 0.) ? eAll_evt / eAll_prof : 0.;
            double ratioAll_cc = (eAll_prof > 0.) ? eAll_chk / eAll_prof : 0.;
            double ratioAll_kc = (eAll_prof > 0.) ? eAll_Kah / eAll_prof : 0.;
            printf("      err ratio vs TProfile SEM -- evt/prof: %.3f   chk/prof: %.3f   kah/prof: %.3f"
                   "  (>1 => SEM underestimated)\n", ratioAll_ec, ratioAll_cc, ratioAll_kc);
        };

        printf("  [%s] <R_proxyJet> summary -- TProfile SEM vs event-mean estimator:\n", familyLabel);
        PrintRingJet("NoCuts",     f.NC_All, f.NC_Pos, f.NC_Neg);
        PrintRingJet("pTCutOnly",  f.PT_All, f.PT_Pos, f.PT_Neg);
        PrintRingJet("DCACutOnly", f.DC_All, f.DC_Pos, f.DC_Neg);
        PrintRingJet("BothCuts",   f.BC_All, f.BC_Pos, f.BC_Neg);
        printf("\n");
    };

    PrintRingJetFamily("WithoutEtaGate", famNG);
    PrintRingJetFamily("WithEtaGate",    famEG);

    // -----------------------------------------------------------------------
    // 5) Write all histograms to disk and close
    // -----------------------------------------------------------------------
    outFile->Write("", TObject::kOverwrite); // This should write to all registered objects from memory.
                                             // Much cleaner than what I used to do before of writing every single histogram explicitly!
    outFile->Close();

    printf("Output written to: %s\n", outputPath);
    printf("Done.\n\n");
}


// Defining a CLING block such that, even if we run outside of ROOT's "root -l -b -q helicityEfficiencyToyModel.cxx" syntax,
// the program is able to be compiled via g++. Notice that ROOT expects to see a function with the name of helicityEfficiencyToyModel,
// to match the file name and allow CLING to compile Just-In-Time (JIT), yet g++ expects to see a function named "main" for it to
// properly compile. Thus, if we don't see CLING (which is the case of compiling via g++), then we define the "main" function below:
#ifndef __CLING__
int main(int argc, char** argv) {
    // 1. Initialize default values identical to the function signature
    long        nLambdas       = 10000000; // This can't go over much more than 2 billion Lambdas! Be careful!
    const char* outputPath     = "helicityEffOutput.root";
    double      Bz_Tesla       = 0.5;
    double      pTmin_Lambda   = 0.0;
    double      pTmax_Lambda   = 10.0;
    double      rapMax_Lambda  = 5.0;
    double      etaMaxDetector = 0.9;
    double      T_thermal      = 0.300;
    double      pTmin_proton   = 0.0;
    double      pTmin_pion     = 0.0;
    double      dcaMin_proton  = 0.00;
    double      dcaMin_pion    = 0.00;
    int         seed           = 0;

    // 2. Safely parse command-line positional arguments sequentially
    try {
        if (argc > 1)  nLambdas       = std::stol(argv[1]); // If needed, go to stoll for string to long-long instead of just to long (stol)
        if (argc > 2)  outputPath     = argv[2];
        if (argc > 3)  Bz_Tesla       = std::stod(argv[3]);
        if (argc > 4)  pTmin_Lambda   = std::stod(argv[4]);
        if (argc > 5)  pTmax_Lambda   = std::stod(argv[5]);
        if (argc > 6)  rapMax_Lambda  = std::stod(argv[6]);
        if (argc > 7)  etaMaxDetector = std::stod(argv[7]);
        if (argc > 8)  T_thermal      = std::stod(argv[8]);
        if (argc > 9)  pTmin_proton   = std::stod(argv[9]);
        if (argc > 10) pTmin_pion     = std::stod(argv[10]);
        if (argc > 11) dcaMin_proton  = std::stod(argv[11]);
        if (argc > 12) dcaMin_pion    = std::stod(argv[12]);
        if (argc > 13) seed           = std::stoi(argv[13]);
    } 
    catch (const std::exception& e) {
        std::cerr << "\nERROR: Command-line argument conversion failed! " << e.what() << std::endl;
        std::cerr << "Usage: " << argv[0] << " [nLambdas] [outputPath] [Bz_Tesla] [pTmin_Lambda] ... [seed]\n" << std::endl;
        return 1;
    }

    // 3. Execute the native macro logic
    helicityEfficiencyToyModel(
        nLambdas, outputPath, Bz_Tesla, pTmin_Lambda, pTmax_Lambda,
        rapMax_Lambda, etaMaxDetector, T_thermal, pTmin_proton,
        pTmin_pion, dcaMin_proton, dcaMin_pion, seed
    );

    return 0;
}
#endif




// end of helicityEfficiencyToyModel.cxx