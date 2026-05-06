// ============================================================
//  dataBasedModelForAEE_Combined.cxx
//
//  Computes the spurious polarization bias from the
//  "Azimuthal Efficiency Effect" (AEE) using two methods:
//
//  PART 1 -- Analytical model:
//    Weights for eta_L, phi_L, and Delta_phi are parametric
//    functions (Gaussians + polynomials) fitted to data.
//
//  PART 2 -- Data-driven model:
//    Weights are taken directly from binned histograms
//    filled with the actual measured distributions.
//
//  Both curves are drawn on the same canvas with a legend.
//
//  Run with: root -l dataBasedModelForAEE_Combined.cxx
// ============================================================

#include <iostream>
#include <vector>
#include <cmath>
#include "TGraph.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TAxis.h"
#include "TLegend.h"

// ============================================================
// PART 1 -- Analytical weight functions
// ============================================================

// Weight for phi_Lambda: two Gaussians fitted to data.
// Convention: Delta_phi = phi_Lambda - phi_Jet.
double weight_phi_Lambda(double phi_L) {
    // Fit parameters from data: mean and width of each Gaussian.
    // sigma = [2]/sqrt(2) because the fit used ((x-[1])/[2])^2.
    double mean1 = 1.64205, sigma1 = 2.92509 / std::sqrt(2.0);
    double mean2 = 4.97346, sigma2 = 2.93542 / std::sqrt(2.0);

    double g1 = exp(-0.5 * pow((phi_L - mean1) / sigma1, 2));
    double g2 = exp(-0.5 * pow((phi_L - mean2) / sigma2, 2));

    return g1 + g2;
}

// Weight for eta_Lambda: twin-peaked, slightly skewed to negative eta.
// Modelled as a flat base plus two Gaussians.
double weight_eta_Lambda(double eta_L) {
    double base       = 24.0;
    double left_peak  = 8.0 * exp(-0.5 * pow((eta_L + 0.7) / 0.3, 2)); // peak ~32 at eta=-0.7
    double right_peak = 6.0 * exp(-0.5 * pow((eta_L - 0.7) / 0.3, 2)); // peak ~30 at eta=+0.7
    return base + left_peak + right_peak;
}

// Weight for Delta_phi: near-side peak at 0, away-side peak at +/-pi.
double weight_Delta_phi(double delta_phi) {
    double base      = 1.0;
    double near_side = 2.0 * exp(-0.5 * pow(delta_phi / 0.4, 2));

    // Away-side peak is periodic at +/-pi; fold distance around pi.
    double dphi_away = fabs(delta_phi) - M_PI;
    double away_side = 1.0 * exp(-0.5 * pow(dphi_away / 0.6, 2));

    return base + near_side + away_side;
}

// ============================================================
// PART 2 -- Raw-histogram weight functions
// ============================================================

// Global histogram pointers (filled in init_histograms()).
TH1D* h_dphi  = nullptr;
TH1D* h_eta_L = nullptr;
TH1D* h_phi_L = nullptr;

// Fill the three histograms from the measured bin contents.
void init_histograms() {

    // --- Delta Phi: 40 bins over [-pi, pi] ---
    // Index 0 and 41 are dummy edge entries (set to 0).
    std::vector<double> dphi_data = {
        0,
        10020202, 9969643,  9862473,  9708729,  9536856,
        9343049,  9145533,  8959756,  8798170,  8672703,
        8578961,  8565351,  8625574,  8757253,  8996051,
        9357137,  9775244,  10161318, 12715456, 18791967,
        18868823, 12893765, 10132769, 9703422,  9341379,
        9002815,  8778303,  8632089,  8578076,  8592720,
        8667630,  8797434,  8954797,  9140359,  9333072,
        9528306,  9697959,  9840104,  9950018,  10020955,
        0
    };
    h_dphi = new TH1D("h_dphi", "", 40, -M_PI, M_PI);
    for (int i = 1; i <= 40; i++) h_dphi->SetBinContent(i, dphi_data[i]);

    // --- Eta Lambda: 50 bins over [-1, 1] ---
    std::vector<double> eta_data = {
        0,
        13792829, 16485943, 19159780, 21708131, 24056965,
        26081141, 27691808, 28806967, 29420026, 29628383,
        29490395, 29130425, 28677729, 28228442, 27799236,
        27378181, 26989194, 26611215, 26264608, 25953587,
        25670748, 25401841, 25026810, 24580036, 24194138,
        24047982, 24209714, 24482498, 24727947, 24923335,
        25184205, 25478833, 25844317, 26239088, 26674797,
        27135384, 27613330, 28115860, 28580151, 28903225,
        28979713, 28658472, 27872154, 26590142, 24783752,
        22569504, 20078673, 17428179, 14733852, 12078547,
        24  // last entry is a residual count; kept for fidelity to source data
    };
    h_eta_L = new TH1D("h_eta_L", "", 50, -1.0, 1.0);
    for (int i = 1; i <= 50; i++) h_eta_L->SetBinContent(i, eta_data[i]);

    // --- Phi Lambda: 50 bins over [0, 2*pi] ---
    std::vector<double> phi_data = {
        0,
        21745660, 23532625, 24358459, 24699952, 25374244,
        25417191, 26725585, 27871285, 26929914, 26796465,
        26594863, 26753481, 27207171, 27019165, 26631242,
        26261256, 24719276, 25391664, 27100858, 25174719,
        24944675, 25660121, 24371851, 24540276, 24219110,
        22984932, 23624190, 21917737, 21803268, 23946859,
        23854239, 24383135, 25939298, 25655165, 26353369,
        25809298, 23832359, 23979304, 25093117, 25778153,
        26665234, 26059304, 25545173, 25314411, 24760420,
        25068612, 23669303, 22887946, 24246152, 24950138,
        12  // last entry is a residual count; kept for fidelity to source data
    };
    h_phi_L = new TH1D("h_phi_L", "", 50, 0.0, 2.0 * M_PI);
    for (int i = 1; i <= 50; i++) h_phi_L->SetBinContent(i, phi_data[i]);
}

// ============================================================
// Main function
// ============================================================
void dataBasedModelForAEE() {

    // ----------------------------------------------------------
    // Shared settings
    // ----------------------------------------------------------
    int    n_eta_t_points = 40;      // number of eta_t scan points
    double eta_t_min      = -0.5;
    double eta_t_max      =  0.5;

    // ----------------------------------------------------------
    // PART 1 -- Analytical integration
    // ----------------------------------------------------------
    int n_grid = 50;  // grid points per dimension (50^3 per eta_t point)

    double eta_L_min = -0.9, eta_L_max = 0.9;
    double phi_L_min =  0.0, phi_L_max = 2.0 * M_PI;
    double dphi_min  = -M_PI, dphi_max =  M_PI;

    double d_eta_L = (eta_L_max - eta_L_min) / n_grid;
    double d_phi_L = (phi_L_max - phi_L_min) / n_grid;
    double d_dphi  = (dphi_max  - dphi_min)  / n_grid;

    TGraph* gr_analytical = new TGraph(n_eta_t_points);

    for (int i = 0; i < n_eta_t_points; i++) {

        double eta_t = eta_t_min + i * (eta_t_max - eta_t_min) / (n_eta_t_points - 1);

        double integral_sum = 0.0;
        double weight_sum   = 0.0;

        // 3D midpoint-rule integration over (eta_L, phi_L, Delta_phi)
        for (int j = 0; j < n_grid; j++) {
            double eta_L = eta_L_min + (j + 0.5) * d_eta_L;

            for (int k = 0; k < n_grid; k++) {
                double phi_L = phi_L_min + (k + 0.5) * d_phi_L;

                for (int m = 0; m < n_grid; m++) {
                    double dphi = dphi_min + (m + 0.5) * d_dphi;

                    // Combined weight from all three distributions
                    double w = weight_eta_Lambda(eta_L)
                             * weight_phi_Lambda(phi_L)
                             * weight_Delta_phi(dphi);

                    // AEE numerator:
                    // Convention: Delta_phi = phi_Lambda - phi_Jet
                    double num = -sinh(eta_L) * cos(2.0 * phi_L - dphi)
                                + sinh(eta_t) * cos(2.0 * phi_L);

                    // AEE denominator (epsilon avoids division by zero
                    // when eta_t == eta_L and dphi == 0 simultaneously)
                    double den_sq = pow(sin(dphi),      2)
                                  + pow(sinh(eta_t),    2)
                                  + pow(sinh(eta_L),    2)
                                  - 2.0 * sinh(eta_t) * sinh(eta_L) * cos(dphi);
                    double den = sqrt(den_sq + 1e-6);

                    double R_AEE = num / den;

                    integral_sum += R_AEE * w;
                    weight_sum   += w;
                }
            }
        }

        double expectation_value = integral_sum / weight_sum;
        gr_analytical->SetPoint(i, eta_t, expectation_value);
    }

    // ----------------------------------------------------------
    // PART 2 -- Data-driven integration (histogram bins as weights)
    // ----------------------------------------------------------
    init_histograms();  // populate h_dphi, h_eta_L, h_phi_L

    TGraph* gr_datadriven = new TGraph(n_eta_t_points);

    for (int i = 0; i < n_eta_t_points; i++) {

        double eta_t = eta_t_min + i * (eta_t_max - eta_t_min) / (n_eta_t_points - 1);

        double integral_sum = 0.0;
        double weight_sum   = 0.0;

        // Loop directly over histogram bins; bin centers replace the grid points.
        for (int j = 1; j <= h_eta_L->GetNbinsX(); j++) {
            double eta_L  = h_eta_L->GetXaxis()->GetBinCenter(j);
            double w_eta  = h_eta_L->GetBinContent(j);

            for (int k = 1; k <= h_phi_L->GetNbinsX(); k++) {
                double phi_L  = h_phi_L->GetXaxis()->GetBinCenter(k);
                double w_phi  = h_phi_L->GetBinContent(k);

                for (int m = 1; m <= h_dphi->GetNbinsX(); m++) {
                    double dphi   = h_dphi->GetXaxis()->GetBinCenter(m);
                    double w_dphi = h_dphi->GetBinContent(m);

                    // Product of the three histogram weights
                    double w = w_eta * w_phi * w_dphi;

                    // Same AEE formula as Part 1
                    double num = -sinh(eta_L) * cos(2.0 * phi_L - dphi)
                                + sinh(eta_t) * cos(2.0 * phi_L);

                    double den_sq = pow(sin(dphi),      2)
                                  + pow(sinh(eta_t),    2)
                                  + pow(sinh(eta_L),    2)
                                  - 2.0 * sinh(eta_t) * sinh(eta_L) * cos(dphi);
                    double den = sqrt(den_sq + 1e-6);

                    double R_AEE = num / den;

                    integral_sum += R_AEE * w;
                    weight_sum   += w;
                }
            }
        }

        double expectation_value = integral_sum / weight_sum;
        gr_datadriven->SetPoint(i, eta_t, expectation_value);
    }

    // ----------------------------------------------------------
    // Plotting -- both curves on the same canvas
    // ----------------------------------------------------------
    TCanvas* c1 = new TCanvas("c1", "AEE Polarization Bias -- Comparison", 900, 600);
    c1->SetGrid();

    // Style for the analytical curve (blue)
    gr_analytical->SetTitle("AEE Spurious Bias;#eta_{Jet}; <R^{AEE}>");
    gr_analytical->SetLineColor(kBlue);
    gr_analytical->SetLineWidth(3);
    gr_analytical->SetLineStyle(1);  // solid

    // Style for the data-driven curve (red, dashed)
    gr_datadriven->SetLineColor(kRed);
    gr_datadriven->SetLineWidth(3);
    gr_datadriven->SetLineStyle(2);  // dashed

    // Draw analytical first (sets the axis range), then overlay data-driven
    gr_analytical->Draw("AL");
    gr_datadriven->Draw("L same");

    // Legend
    // TLegend* leg = new TLegend(0.15, 0.70, 0.55, 0.88);
    TLegend* leg = new TLegend(0.55, 0.75, 0.88, 0.88);
    leg->SetBorderSize(1);
    leg->SetFillColor(0);
    leg->AddEntry(gr_analytical, "Analytical model (Part 1)", "l");
    leg->AddEntry(gr_datadriven, "Data-driven model (Part 2)", "l");
    leg->Draw();

    c1->Update();
}