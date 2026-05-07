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
//   NoCuts/
//     EtaPos/   -- all histograms for eta_Lambda > 0
//     EtaNeg/   -- all histograms for eta_Lambda < 0
//     All/      -- both eta halves combined
//   pTCutOnly/
//     EtaPos/ EtaNeg/ All/
//   DCACutOnly/
//     EtaPos/ EtaNeg/ All/
//   BothCuts/
//     EtaPos/ EtaNeg/ All/
//   Kinematics/
//     (Lambda pT spectrum, rapidity, decay radius, daughter pT before any cut)
//
// Each leaf directory contains:
//   h2d_cosTheta_phi  -- 2D: cos(theta*) vs phi*  [main diagnostic]
//   h1d_cosTheta      -- 1D projection of cos(theta*)
//   h1d_phi           -- 1D projection of phi*
//   h1d_ringProxy     -- distribution of R_proxy values
//   pRingProxy        -- TProfile: integrated <R_proxy>
//   pRingProxyVsEta   -- TProfile: <R_proxy> vs Lambda eta
//   pRingProxyVsPt    -- TProfile: <R_proxy> vs Lambda pT
//   h1d_decayRadius   -- transverse decay vertex radius [cm]
//   h1d_pT_proton     -- proton pT in lab frame
//   h1d_pT_pion       -- pion pT in lab frame
//   h1d_DCA_proton    -- proton DCA_xy to PV
//   h1d_DCA_pion      -- pion DCA_xy to PV
//   h1d_pT_lambda     -- Lambda pT
//   h1d_eta_lambda    -- Lambda eta
//
//   Reference analysis code: lambdaJetPolarizationIonsDerived.cxx (ALICE O2)
// ==========================================================================

// --- Standard ROOT headers ---
#include <TFile.h>
#include <TDirectory.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TGenPhaseSpace.h>
#include <TRandom3.h>
#include <TMath.h>
#include <TString.h>
#include <TROOT.h>

// --- Standard C++ headers ---
#include <cmath>
#include <cstdio>
#include <cstring>


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
static const double kMinSinTheta = 1.e-4;

constexpr double TwoPi = TMath::TwoPI();
constexpr double Pi = TMath::Pi();


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
 * @param rng    Pointer to TRandom3 random-number generator (if not passed by reference, the RNG list of already-used values would not be updated properly!)
 * @param T      Boltzmann temperature parameter [GeV]
 * @param pTmin  Lower bound of sampling interval [GeV/c]
 * @param pTmax  Upper bound of sampling interval [GeV/c]
 * @param fmax   Envelope maximum used in rejection sampling
 *
 * @return One accepted Lambda transverse momentum pT [GeV/c]
 * ==========================================================================
 */
static double SampleLambdaPt(TRandom3* rng, double T, double pTmin, double pTmax, double fmax){
    // -- Repeat until one trial pT is accepted --
    while (true) {
        // -- Generate a trial pT uniformly inside the requested interval --
        double pt = pTmin + (pTmax - pTmin) * rng->Rndm();
        // -- Compute transverse mass corresponding to the trial pT --
        double mT = std::sqrt(pt * pt + kMassLambda * kMassLambda);
        // -- Evaluate the target Boltzmann weight at this pT --
        double f = pt * std::exp(-mT / T);

        // -- Accept the trial pT with probability f/fmax --
        //    Geometrically: accept if a random point inside the rectangle
        //    [pTmin,pTmax] x [0,fmax] falls below the target curve f(pT)
        if (rng->Rndm() * fmax < f)
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
    TH1D*    h1d_ringProxy;      // distribution of R_proxy
    TProfile* pRingProxy;        // integrated <R_proxy> (one bin at x=0)
    TProfile* pRingProxyVsEta;   // <R_proxy> vs Lambda pseudorapidity
    TProfile* pRingProxyVsPt;    // <R_proxy> vs Lambda pT

    // -- Decay geometry and daughter kinematics (useful for diagnosis) --
    TH1D*    h1d_decayRadius;    // transverse decay vertex radius [cm]
    TH1D*    h1d_pT_proton;      // proton pT in lab frame [GeV/c]
    TH1D*    h1d_pT_pion;        // pion pT in lab frame [GeV/c]
    TH1D*    h1d_DCA_proton;     // proton DCA_xy to PV [cm]
    TH1D*    h1d_DCA_pion;       // pion DCA_xy to PV [cm]

    // -- Lambda kinematics --
    TH1D*    h1d_pT_lambda;      // Lambda pT [GeV/c]
    TH1D*    h1d_eta_lambda;     // Lambda pseudorapidity
};


// ==========================================================================
/**
 * @brief Allocates and initialises all histograms inside the given TDirectory.
 *
 * The directory must already exist.  Calls dir->cd() before creating objects
 * so ROOT associates them with the correct directory.
 *
 * @param dir  Pointer to an existing TDirectory in which histograms are booked.
 * @return     A filled ScenarioHistos struct with all histogram pointers set.
 */
// ==========================================================================
static ScenarioHistos BookScenario(TDirectory* dir)
{
    ScenarioHistos h;

    dir->cd();  // All subsequent "new TH1D/TH2D" calls attach to this dir

    // ---------- Main diagnostic: (cos theta*, phi*) space ----------
    // Expected: flat before cuts; asymmetric after cuts.
    h.h2d_cosTheta_phi = new TH2D("h2d_cosTheta_phi", "Proton emission in Lambda frame;"
        "cos(#theta*) [proton vs #Lambda dir];"
        "#phi* [proton azimuth around #Lambda] [rad]",
        50, -1., 1.,          // 50 bins in cos(theta*) in [-1, 1]
        64, -TMath::Pi(), TMath::Pi()); // 64 bins in phi* in [-pi, pi]
    h.h1d_cosTheta = new TH1D("h1d_cosTheta", "cos(#theta*) -- proton vs #Lambda direction; cos(#theta*);Counts", 50, -1., 1.);
    h.h1d_phi = new TH1D("h1d_phi", "#phi* -- proton azimuth around #Lambda axis; #phi* [rad];Counts",64, -TMath::Pi(), TMath::Pi());

    // ---------- Ring observable proxy ----------
    // Range is set to [-3/alpha, +3/alpha] which covers the full prefactor.
    double rmax = kPolPrefactor * 1.05; // Slight margin beyond max value
    h.h1d_ringProxy = new TH1D("h1d_ringProxy", "Ring observable proxy (z-hat as jet direction); R_{proxy};Counts", 120, -rmax, rmax);
    h.pRingProxy = new TProfile("pRingProxy", "Integrated <R_{proxy}> (single-bin TProfile); bin;<R_{proxy}>", 1, -0.5, 0.5);
    h.pRingProxyVsEta = new TProfile("pRingProxyVsEta", "<R_{proxy}> vs #Lambda pseudorapidity; #eta_{#Lambda};<R_{proxy}>", 18, -0.9, 0.9); // 18 bins, 0.1 wide each
    h.pRingProxyVsPt = new TProfile("pRingProxyVsPt", "<R_{proxy}> vs #Lambda p_{T}; p_{T}^{#Lambda} [GeV/c];<R_{proxy}>", 20, 0., 5.); // 20 bins of 0.25 GeV/c each

    // ---------- Decay vertex and daughter kinematics ----------
    h.h1d_decayRadius = new TH1D("h1d_decayRadius", "Transverse decay vertex radius; r_{decay} = #sqrt{x_{v}^{2}+y_{v}^{2}} [cm];Counts", 150, 0., 150.);     // 0 to 150 cm (covers most of ITS/TPC inner field cage)
    h.h1d_pT_proton = new TH1D("h1d_pT_proton", "Proton p_{T} in lab frame; p_{T}^{proton} [GeV/c];Counts", 100, 0., 5.);
    h.h1d_pT_pion = new TH1D("h1d_pT_pion", "Pion p_{T} in lab frame; p_{T}^{#pi} [GeV/c];Counts", 100, 0., 5.);
    h.h1d_DCA_proton = new TH1D("h1d_DCA_proton", "Proton DCA_{xy} to primary vertex; DCA_{xy}^{proton} [cm];Counts", 200, 0., 10.); // Fine binning near 0 to see the cut clearly
    h.h1d_DCA_pion = new TH1D("h1d_DCA_pion", "Pion DCA_{xy} to primary vertex; DCA_{xy}^{#pi} [cm];Counts", 200, 0., 10.);
    // ---------- Lambda kinematics ----------
    h.h1d_pT_lambda = new TH1D("h1d_pT_lambda", "#Lambda p_{T}; p_{T}^{#Lambda} [GeV/c];Counts", 100, 0., 10.);
    h.h1d_eta_lambda = new TH1D("h1d_eta_lambda", "#Lambda pseudorapidity; #eta_{#Lambda};Counts", 36, -0.9, 0.9); // 0.05 per bin, matching typical ALICE acceptance

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
 * @param decayR      Transverse decay vertex radius [cm].
 * @param pT_proton   Proton pT in the lab frame [GeV/c].
 * @param pT_pion     Pion pT in the lab frame [GeV/c].
 * @param dca_proton  Proton DCA_xy to PV [cm].
 * @param dca_pion    Pion DCA_xy to PV [cm].
 * @param pT_lambda   Lambda pT [GeV/c].
 * @param eta_lambda  Lambda pseudorapidity.
 */
// ==========================================================================
static void FillScenario(ScenarioHistos& h,
                          double cosTheta,
                          double phi,
                          double ringProxy,
                          double decayR,
                          double pT_proton,
                          double pT_pion,
                          double dca_proton,
                          double dca_pion,
                          double pT_lambda,
                          double eta_lambda)
{
    // Main diagnostic histograms
    h.h2d_cosTheta_phi->Fill(cosTheta, phi);
    h.h1d_cosTheta->Fill(cosTheta);
    h.h1d_phi->Fill(phi);

    // Ring proxy
    h.h1d_ringProxy->Fill(ringProxy);
    h.pRingProxy->Fill(0., ringProxy); // integrated average
    h.pRingProxyVsEta->Fill(eta_lambda, ringProxy);
    h.pRingProxyVsPt->Fill(pT_lambda, ringProxy);

    // Kinematics
    h.h1d_decayRadius->Fill(decayR);
    h.h1d_pT_proton->Fill(pT_proton);
    h.h1d_pT_pion->Fill(pT_pion);
    h.h1d_DCA_proton->Fill(dca_proton);
    h.h1d_DCA_pion->Fill(dca_pion);
    h.h1d_pT_lambda->Fill(pT_lambda);
    h.h1d_eta_lambda->Fill(eta_lambda);
}


// ==========================================================================
/**
 * @brief Creates a three-level directory structure inside the output file and
 *        books ScenarioHistos in each leaf directory.
 *
 * Directory layout created under the file root:
 *   scenarioName/EtaPos/
 *   scenarioName/EtaNeg/
 *   scenarioName/All/
 *
 * Results are returned via the output references hPos, hNeg, and hAll.
 *
 * @param outFile       Pointer to the open output TFile.
 * @param scenarioName  Name of the top-level scenario directory to create.
 * @param hPos          Output: ScenarioHistos booked for eta_Lambda >= 0.
 * @param hNeg          Output: ScenarioHistos booked for eta_Lambda < 0.
 * @param hAll          Output: ScenarioHistos booked for both eta halves combined.
 */
// ==========================================================================
static void CreateSubdirs(TFile*             outFile,
                           const char*        scenarioName,
                           ScenarioHistos&    hPos,
                           ScenarioHistos&    hNeg,
                           ScenarioHistos&    hAll)
{
    // Create the top-level scenario directory
    TDirectory* dirScen = outFile->mkdir(scenarioName);

    // Create sub-directories for each eta half and the combined sample
    TDirectory* dirPos = dirScen->mkdir("EtaPos");
    TDirectory* dirNeg = dirScen->mkdir("EtaNeg");
    TDirectory* dirAll = dirScen->mkdir("All");

    // Book histograms inside each sub-directory
    hPos = BookScenario(dirPos);
    hNeg = BookScenario(dirNeg);
    hAll = BookScenario(dirAll);

    // Return to root directory of file so subsequent mkdir calls start there
    outFile->cd();
}


// ==========================================================================
/**
 * @brief Generates unpolarized Lambda decays and studies acceptance-induced
 *        fake-polarization effects by varying kinematic and topological cuts.
 *
 * @details
 * Four cut scenarios are filled per run: NoCuts, pTCutOnly, DCACutOnly, and
 * BothCuts.  Within each scenario, histograms are split into EtaPos, EtaNeg,
 * and All sub-directories (see file header for the full output structure).
 *
 * @param nLambdas       Number of Lambdas to simulate (default: 10,000,000).
 * @param outputPath     Output ROOT file path (default: helicityEffOutput.root).
 * @param Bz_Tesla       Magnetic field along z [T], positive = ALICE default
 *                       (default: +0.5 T).
 * @param pTmin_Lambda   Minimum Lambda pT for generation [GeV/c] (default: 0.0).
 * @param pTmax_Lambda   Maximum Lambda pT for generation [GeV/c] (default: 10.0).
 * @param rapMax_Lambda  Maximum |rapidity| for GENERATING Lambdas (default: 5.0).
 * @param etaMaxDetector Maximum |eta| for detecting Lambda daughters (default: 0.9,
 *                       matching the ALICE inner barrel). Not applied to the Lambda
 *                       itself; cutting on Lambda eta would bias the distribution.
 *                       Detectors only see charged particles after all.
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
    printf("  Lambda pT range       : [%.2f, %.2f] GeV/c\n",
           pTmin_Lambda, pTmax_Lambda);
    printf("  Lambda |rapidity| max : %.2f\n", rapMax_Lambda);
    printf("  Detector |eta| max : %.2f\n", etaMaxDetector);
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

    // Create histogram containers for each (scenario, eta) combination.
    // Naming convention: hXxx_Pos = eta > 0, hXxx_Neg = eta < 0, hXxx_All = both.
    ScenarioHistos hNC_Pos, hNC_Neg, hNC_All;         // No cuts
    ScenarioHistos hPT_Pos, hPT_Neg, hPT_All;         // pT cut only
    ScenarioHistos hDC_Pos, hDC_Neg, hDC_All;         // DCA cut only
    ScenarioHistos hBC_Pos, hBC_Neg, hBC_All;         // Both cuts

    // CreateSubdirs makes the TDirectory hierarchy and books all histograms:
    CreateSubdirs(outFile, "NoCuts",    hNC_Pos, hNC_Neg, hNC_All);
    CreateSubdirs(outFile, "pTCutOnly", hPT_Pos, hPT_Neg, hPT_All);
    CreateSubdirs(outFile, "DCACutOnly",hDC_Pos, hDC_Neg, hDC_All);
    CreateSubdirs(outFile, "BothCuts",  hBC_Pos, hBC_Neg, hBC_All);

    // Extra kinematics directory: Lambda distribution BEFORE any decay cuts,
    // only subject to the generation acceptance (pT and rapidity windows).
    TDirectory* dirKin = outFile->mkdir("Kinematics");
    dirKin->cd();
    TH1D* hKin_pT_lambda  = new TH1D("hKin_pT_lambda", "Generated #Lambda p_{T};p_{T}^{#Lambda} [GeV/c];Counts", 100, 0., 10.);
    TH1D* hKin_rap_lambda = new TH1D("hKin_rap_lambda", "Generated #Lambda rapidity;y_{#Lambda};Counts", 36, -0.9, 0.9);
    TH1D* hKin_eta_lambda = new TH1D("hKin_eta_lambda", "Generated #Lambda pseudorapidity;#eta_{#Lambda};Counts", 36, -0.9, 0.9);
    TH1D* hKin_phi_lambda = new TH1D("hKin_phi_lambda", "Generated #Lambda azimuth;#phi_{#Lambda} [rad];Counts", 64, -TMath::Pi(), TMath::Pi());
    TH1D* hKin_decayR     = new TH1D("hKin_decayR", "All transverse decay radii;r_{decay} [cm];Counts", 150, 0., 150.);
    TH1D* hKin_pT_proton  = new TH1D("hKin_pT_proton", "All proton p_{T} (pre-cut);p_{T}^{p} [GeV/c];Counts", 100, 0., 5.);
    TH1D* hKin_pT_pion    = new TH1D("hKin_pT_pion", "All pion p_{T} (pre-cut);p_{T}^{#pi} [GeV/c];Counts", 100, 0., 5.);
    TH1D* hKin_DCA_proton = new TH1D("hKin_DCA_proton", "All proton DCA_{xy} (pre-cut);DCA_{xy}^{p} [cm];Counts", 200, 0., 10.);
    TH1D* hKin_DCA_pion   = new TH1D("hKin_DCA_pion", "All pion DCA_{xy} (pre-cut);DCA_{xy}^{#pi} [cm];Counts", 200, 0., 10.);
    outFile->cd(); // Back to root of output file

    // -----------------------------------------------------------------------
    // 2) Initialise the random number generator and phase-space decay engine
    // -----------------------------------------------------------------------
    TRandom3 rng(seed);

    // TGenPhaseSpace: decays a parent 4-vector into N daughters according to
    // Lorentz-invariant phase space (flat in LIPS = unpolarized).
    TGenPhaseSpace decayGen;

    // Daughter masses array for Lambda -> p + pi-
    Double_t daughterMasses[2] = { kMassProton, kMassPion };

    // -----------------------------------------------------------------------
    // 3) Main event loop
    // -----------------------------------------------------------------------
    printf("Starting generation of %ld Lambdas...\n", nLambdas);

    long nGenerated   = 0; // Lambdas for which decay was successfully set up
    long nCollinear   = 0; // Lambdas skipped because pT/|p| < kMinSinTheta

    // Progress report every 10%
    long nProgress = nLambdas / 10;
    if (nProgress < 1) nProgress = 1;

    // Calculating useful variables for the SampleLambdaPt function:
    // We only sample a Boltzmann function, which has an analytic maximum that can be calculated for the rejection sampling below
    double ptPeak = ComputeThermalPtMaximum(T_thermal);
    double mTpeak = std::sqrt(ptPeak * ptPeak + kMassLambda * kMassLambda);
    double fmax = 1.10 * ptPeak * std::exp(-mTpeak / T_thermal);

    for (long iLam = 0; iLam < nLambdas; ++iLam) {
        // ---- Progress printout ----
        if ((iLam + 1) % nProgress == 0)
            printf("  %ld / %ld (%.0f%%)\n",iLam + 1, nLambdas, 100. * (iLam + 1) / nLambdas);

        // ==================================================================
        // 3.1  Generate Lambda 4-momentum
        // ==================================================================
        // Sample pT from thermal mT spectrum [GeV/c]
        double pT_lam = SampleLambdaPt(&rng, T_thermal, pTmin_Lambda, pTmax_Lambda, fmax);

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

        // Build TLorentzVector for Lambda (used by TGenPhaseSpace)
        TLorentzVector lv_lam(px_lam, py_lam, pz_lam, E_lam);

        // Fill pre-cut kinematics
        hKin_pT_lambda->Fill(pT_lam);
        hKin_rap_lambda->Fill(rap_lam);
        hKin_eta_lambda->Fill(eta_lam);
        hKin_phi_lambda->Fill(wrapToPiFast(phi_lam));

        // ==================================================================
        // 3.2  Build Lambda frame axes
        // ==================================================================
        // e1 = Lambda unit momentum vector (forward-backward axis)
        TVector3 e1(px_lam / p_lam, py_lam / p_lam, pz_lam / p_lam);

        // e2 = (p_Lambda x z_hat) / |p_Lambda x z_hat| (transverse, left-right axis)
            // Explicitly:  e2 = (py, -px, 0) / pT
        // Reject nearly beam-collinear Lambdas (pT/|p| is too small):
        if (pT_lam / p_lam < kMinSinTheta) {
            ++nCollinear;
            continue;
        }
        TVector3 e2(py_lam / pT_lam, -px_lam / pT_lam, 0.);

        // e3 = e1 x e2  (in-out axis)
        // Explicitly:  e3 = (pz*px, pz*py, -pT^2) / (|p|*pT)
        TVector3 e3 = e1.Cross(e2);  // Unit by construction since e1, e2 are unit orthonormal

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
        double l_proper = -kCTauLambda * std::log(rng.Rndm());

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
        TLorentzVector lv_proton = *(decayGen.GetDecay(0));
        TLorentzVector lv_pion   = *(decayGen.GetDecay(1));

        // Lab-frame transverse momenta of daughters [GeV/c]
        double px_p  = lv_proton.Px();
        double py_p  = lv_proton.Py();
        double pT_p  = lv_proton.Pt();

        double px_pi = lv_pion.Px();
        double py_pi = lv_pion.Py();
        double pT_pi = lv_pion.Pt();

        // Fill pre-cut daughter kinematics
        hKin_pT_proton->Fill(pT_p);
        hKin_pT_pion->Fill(pT_pi);

        // Cutting in daughter eta:
        // (TODO: actually apply this cut with etaMaxDetector!
        //  All other results that previously existed in this code make no sense if you don't implement this!
        //  Could maybe create a whole other family of cuts such as the noCuts, pTCuts, both cuts and so on TDirectory families, but now for the Eta cut?)
        //  Or even better: we should create a folder that contains all the other families without this eta cut, and then clone all existing families
        //  to properly include this eta cut.
        //  All of the previous histograms are somewhat inconsistent: we generate Lambdas in a small rapidity range (was 0.9 by default. Now 5.0 to 
        //  guarantee we include forward Lambdas that could still be detected in the Inner Barrel acceptance of |eta|<0.9), yet we
        //  allow for their daughters to be created with all possible rapidities and don't consider that the original generated
        //  rapidity distribution (|y_Lambda|<0.9) was smaller than it should be. We were cutting in rapidity but forgot to cut
        //  in eta to make everything consistent!
        double eta_p = lv_proton.Eta();
        double eta_pi = lv_pion.Eta();

        // ==================================================================
        // 3.5  Boost proton to Lambda rest frame
        // ==================================================================
        // Boost the proton 4-vector from the lab frame into the Lambda rest frame.
        // In TLorentzVector, BoostVector() returns beta = p/E (rest --> lab).
        // To go from lab --> rest frame, we must use the opposite boost: -p/E.
        TLorentzVector lv_proton_star = lv_proton;
        lv_proton_star.Boost(-lv_lam.BoostVector());

        // Unit vector of proton momentum in Lambda rest frame
        TVector3 p_star_unit = lv_proton_star.Vect().Unit();

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
        // t_hat x lambda_unit = z_hat x (lx, ly, lz) = (-ly, lx, 0)
        // |z_hat x lambda_unit| = sqrt(lx^2 + ly^2) = pT / |p| = sinTheta_lam
        //
        // Normalized cross product = (-ly, lx, 0) / sinTheta_lam
        //                          = (-py, px, 0) / pT     (after scaling)
        //
        // ring_proxy = (3/alpha) * p_star_unit . (-py, px, 0) / pT
        //
        // Note: this equals -(3/alpha) * p_star_unit . e2 = -(3/alpha) * sin(theta*) * cos(phi*)
        TVector3 crossVec(-py_lam / pT_lam, px_lam / pT_lam, 0.); // crossVec is already unit-normalized (|(-py, px, 0)| / pT = 1)
        double ringProxy = kPolPrefactor * p_star_unit.Dot(crossVec);

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
        // -- pT cut: both daughters must exceed minimum pT --
        bool passPtCut = (pT_p  >= pTmin_proton) && (pT_pi >= pTmin_pion);

        // -- DCA cut: both daughters must be sufficiently displaced from PV --
        // (This selects genuine secondary tracks in data, i.e. rejects primaries.
        //  We lack actual primary-secondary distinction here, so it just emulates it)
        bool passDcaCut = (dca_proton >= dcaMin_proton) && (dca_pion >= dcaMin_pion);

        // -- Eta > 0 flag --
        bool etaPos = (eta_lam >= 0.);

        // ==================================================================
        // 3.10  Fill histogram scenarios
        // ==================================================================
        // We pack the computed quantities into one call per scenario per event.
        // Helper macro to avoid repeating the 12-argument call four times:
        #define FILL_BOTH_ETA(HPOS, HNEG, HALL) \
            do { \
                ScenarioHistos& heta = etaPos ? (HPOS) : (HNEG); \
                FillScenario(heta,  cosTheta, phi_star, ringProxy, decayR, pT_p, pT_pi, dca_proton, dca_pion, pT_lam, eta_lam); \
                FillScenario((HALL), cosTheta, phi_star, ringProxy, decayR, pT_p, pT_pi, dca_proton, dca_pion, pT_lam, eta_lam); \
            } while (0)

        // Scenario 1: No cuts (always filled -- serves as a flat-distribution
        //             bug-check; any asymmetry here is a code error)
        FILL_BOTH_ETA(hNC_Pos, hNC_Neg, hNC_All);

        // Scenario 2: pT cut only
        if (passPtCut)
            FILL_BOTH_ETA(hPT_Pos, hPT_Neg, hPT_All);

        // Scenario 3: DCA cut only
        if (passDcaCut)
            FILL_BOTH_ETA(hDC_Pos, hDC_Neg, hDC_All);

        // Scenario 4: Both cuts
        if (passPtCut && passDcaCut)
            FILL_BOTH_ETA(hBC_Pos, hBC_Neg, hBC_All);

        // Remove local macro to avoid leakage outside the loop
        #undef FILL_BOTH_ETA

        ++nGenerated;

    } // end main loop over Lambdas

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

    printf("  Cut statistics (fraction of NoCuts sample):\n");
    double nNC = (double)hNC_All.h1d_pT_lambda->GetEntries();
    double nPT = (double)hPT_All.h1d_pT_lambda->GetEntries();
    double nDC = (double)hDC_All.h1d_pT_lambda->GetEntries();
    double nBC = (double)hBC_All.h1d_pT_lambda->GetEntries();
    printf("  NoCuts    : %.0f  (100%%)\n",       nNC);
    printf("  pTCutOnly : %.0f  (%.1f%%)\n",     nPT, 100.*nPT/nNC);
    printf("  DCACutOnly: %.0f  (%.1f%%)\n",     nDC, 100.*nDC/nNC);
    printf("  BothCuts  : %.0f  (%.1f%%)\n",     nBC, 100.*nBC/nNC);
    printf("\n");

    // Print integrated ring proxy values for quick sanity check:
    printf("  Integrated <R_proxy> summary  (expected: 0 for NoCuts):\n");
    auto PrintRing = [](const char* label, ScenarioHistos& hA,
                        ScenarioHistos& hP, ScenarioHistos& hN) {
        double rAll  = hA.pRingProxy->GetBinContent(1);
        double eAll  = hA.pRingProxy->GetBinError(1);
        double rPos  = hP.pRingProxy->GetBinContent(1);
        double ePos  = hP.pRingProxy->GetBinError(1);
        double rNeg  = hN.pRingProxy->GetBinContent(1);
        double eNeg  = hN.pRingProxy->GetBinError(1);
        printf("  %-12s  All: %+.4e +/- %.4e  EtaPos: %+.4e +/- %.4e  EtaNeg: %+.4e +/- %.4e\n", label, rAll, eAll, rPos, ePos, rNeg, eNeg);
    };
    PrintRing("NoCuts",     hNC_All, hNC_Pos, hNC_Neg);
    PrintRing("pTCutOnly",  hPT_All, hPT_Pos, hPT_Neg);
    PrintRing("DCACutOnly", hDC_All, hDC_Pos, hDC_Neg);
    PrintRing("BothCuts",   hBC_All, hBC_Pos, hBC_Neg);
    printf("\n");

    // -----------------------------------------------------------------------
    // 5) Write all histograms to disk and close
    // -----------------------------------------------------------------------
    outFile->Write("", TObject::kOverwrite);
    outFile->Close();

    printf("Output written to: %s\n", outputPath);
    printf("Done.\n\n");
}
// end of helicityEfficiencyToyModel.cxx