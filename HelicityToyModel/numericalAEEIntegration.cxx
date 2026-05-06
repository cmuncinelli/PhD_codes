// Save this file as: numericalAEEIntegration.cxx
// Run in terminal with: root -l numericalAEEIntegration.cxx

#include <iostream>
#include <cmath>
#include "TGraph.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TMath.h"

// 1. Parameterization for the Jet correlation Delta_phi distribution
double weight_Delta_phi(double delta_phi) {
    double base = 1.0;
    double near_side = 2.0 * exp(-0.5 * pow(delta_phi / 0.4, 2));
    double dphi_away = fabs(delta_phi) - M_PI;
    double away_side = 1.0 * exp(-0.5 * pow(dphi_away / 0.6, 2));
    return base + near_side + away_side;
}

// 2. Parameterization for the eta_Lambda distribution
double weight_eta_Lambda(double eta_L) {
    double base = 24.0;
    double left_peak  = 8.0 * exp(-0.5 * pow((eta_L + 0.7)/0.3, 2)); 
    double right_peak = 6.0 * exp(-0.5 * pow((eta_L - 0.7)/0.3, 2)); 
    return base + left_peak + right_peak;
}

// Define the empirical Delta Eta distribution (42 bins, from -1.5 to 1.5)
std::vector<double> deta_data = {
    10593, 77913, 403417, 1208377, 2445989, 3863489, 5351153, 6852882,
    8319063, 9730288, 11099972, 12402580, 13612437, 14565267, 15028795,
    14964161, 14631424, 14616621, 15635345, 17377582, 19014211, 18976490,
    17256221, 15469078, 14402665, 14363705, 14629830, 14608905, 14094588,
    13123680, 11946987, 10678902, 9351462, 7962911, 6512273, 5022804,
    3552551, 2194710, 1044481, 326157, 58664, 7598
};

// Interpolation function for Delta Eta (Delta Eta = eta_Lambda - eta_Jet)
double weight_Delta_eta(double delta_eta) {
    if (delta_eta <= -1.5 || delta_eta >= 1.5) return 0.0;
    
    int n_bins = 42;
    double bin_width = 3.0 / n_bins;
    
    // Find the exact floating point bin index
    double exact_bin = (delta_eta + 1.5) / bin_width;
    int bin_idx = std::floor(exact_bin);
    
    // Simple linear interpolation between bin centers
    if (bin_idx < 0) return deta_data[0];
    if (bin_idx >= n_bins - 1) return deta_data[n_bins - 1];
    
    double t = exact_bin - bin_idx - 0.5; // Offset to bin center
    if (t < 0 && bin_idx > 0) {
        bin_idx--;
        t += 1.0;
    } else if (t > 0 && bin_idx == n_bins - 1) {
        t = 0.0; 
    }
    
    return deta_data[bin_idx] * (1.0 - t) + deta_data[bin_idx + 1] * t;
}

void numericalAEEIntegration() {
    int n_grid = 200; 
    
    double eta_L_min = -0.9, eta_L_max = 0.9;
    double dphi_min  = -M_PI, dphi_max = M_PI;
    double eta_t_min = -0.5, eta_t_max = 0.5;
    
    double d_eta_L = (eta_L_max - eta_L_min) / n_grid;
    double d_dphi  = (dphi_max - dphi_min) / n_grid;
    double d_eta_t = (eta_t_max - eta_t_min) / n_grid;
    
    int n_scan_points = 50;

    TGraph* gr_flat  = new TGraph(); 
    TGraph* gr_param = new TGraph(); 
    TGraph* gr_flat_div  = new TGraph();
    TGraph* gr_param_div = new TGraph();
    TGraph* gr_marginal_flat  = new TGraph(n_scan_points);
    TGraph* gr_marginal_param = new TGraph(n_scan_points);
    
    // =========================================================================
    // LOOP 1: Evaluate R as a function of eta_t
    // =========================================================================
    int div_point_idx = 0;
    for (int i = 0; i < n_scan_points; i++) {
        double eta_t = eta_t_min + i * ((eta_t_max - eta_t_min) / (n_scan_points - 1));
        
        double integral_flat = 0.0, weight_sum_flat = 0.0;
        double integral_param = 0.0, weight_sum_param = 0.0;
        
        for (int j = 0; j < n_grid; j++) {
            double eta_L = eta_L_min + (j + 0.5) * d_eta_L;
            for (int m = 0; m < n_grid; m++) {
                double dphi = dphi_min + (m + 0.5) * d_dphi;
                
                double num = - sinh(eta_L) * cos(dphi) + sinh(eta_t);
                double den_sq = pow(sin(dphi), 2) + pow(sinh(eta_t), 2) + pow(sinh(eta_L), 2) 
                                - 2.0*sinh(eta_t)*sinh(eta_L)*cos(dphi);
                double den = sqrt(den_sq + 1e-6); 
                double R_kinematic = num / den;
                
                double w_flat = 1.0; 
                integral_flat += R_kinematic * w_flat;
                weight_sum_flat += w_flat;
                
                double w_param = weight_eta_Lambda(eta_L) * weight_Delta_phi(dphi);
                integral_param += R_kinematic * w_param;
                weight_sum_param += w_param;
            }
        }
        
        double y_flat = integral_flat / weight_sum_flat;
        double y_param = integral_param / weight_sum_param;
        
        gr_flat->SetPoint(i, eta_t, y_flat);
        gr_param->SetPoint(i, eta_t, y_param);
        
        if (fabs(eta_t) > 1e-4) {
            gr_flat_div->SetPoint(div_point_idx, eta_t, y_flat / eta_t);
            gr_param_div->SetPoint(div_point_idx, eta_t, y_param / eta_t);
            div_point_idx++;
        }
    }
    
    // =========================================================================
    // LOOP 2: Evaluate Marginal R as a function of eta_Lambda
    // =========================================================================
    for (int i = 0; i < n_scan_points; i++) {
        double eta_L = eta_L_min + i * ((eta_L_max - eta_L_min) / (n_scan_points - 1));
        
        double integral_flat = 0.0, weight_sum_flat = 0.0;
        double integral_param = 0.0, weight_sum_param = 0.0;
        
        for (int j = 0; j < n_grid; j++) {
            double eta_t = eta_t_min + (j + 0.5) * d_eta_t;
            for (int m = 0; m < n_grid; m++) {
                double dphi = dphi_min + (m + 0.5) * d_dphi;
                
                double num = - sinh(eta_L) * cos(dphi) + sinh(eta_t);
                double den_sq = pow(sin(dphi), 2) + pow(sinh(eta_t), 2) + pow(sinh(eta_L), 2) 
                                - 2.0*sinh(eta_t)*sinh(eta_L)*cos(dphi);
                double den = sqrt(den_sq + 1e-6); 
                double R_kinematic = num / den;
                
                double w_flat = 1.0; 
                integral_flat += R_kinematic * w_flat;
                weight_sum_flat += w_flat;
                
                double w_param = weight_Delta_phi(dphi); 
                integral_param += R_kinematic * w_param;
                weight_sum_param += w_param;
            }
        }
        
        gr_marginal_flat->SetPoint(i, eta_L, integral_flat / weight_sum_flat);
        gr_marginal_param->SetPoint(i, eta_L, integral_param / weight_sum_param);
    }
    
    // =========================================================================
    // PLOTTING
    // =========================================================================
    
    auto style_graph = [](TGraph* g_flat, TGraph* g_param, const char* title) {
        g_flat->SetTitle(title);
        g_flat->SetLineColor(kRed);
        g_flat->SetLineStyle(2); 
        g_flat->SetLineWidth(3);
        
        g_param->SetLineColor(kBlue);
        g_param->SetLineWidth(3);
    };

    // A robust function to auto-scale axes without crashing ROOT
    auto adjust_axes = [](TCanvas* c, TGraph* g1, TGraph* g2) {
        c->Update(); // Force ROOT to build the TAxis first
        double min1 = TMath::MinElement(g1->GetN(), g1->GetY());
        double max1 = TMath::MaxElement(g1->GetN(), g1->GetY());
        double min2 = TMath::MinElement(g2->GetN(), g2->GetY());
        double max2 = TMath::MaxElement(g2->GetN(), g2->GetY());
        
        double ymin = std::min(min1, min2);
        double ymax = std::max(max1, max2);
        
        // Use absolute difference for the margin!
        double margin = 0.2 * fabs(ymax - ymin);
        if (margin < 1e-5) margin = 0.1 * fabs(ymin);
        if (margin < 1e-5) margin = 0.1; // Fallback if everything is exactly 0
        
        g1->GetYaxis()->SetRangeUser(ymin - margin, ymax + margin);
        c->Modified();
        c->Update();
    };

    style_graph(gr_flat, gr_param, "Kinematic Bias vs #eta_{Jet}; #eta_{Jet}; <#mathcal{R}>");
    style_graph(gr_flat_div, gr_param_div, "Ratio Check: <#mathcal{R}> / #eta_{Jet}; #eta_{Jet}; <#mathcal{R}> / #eta_{Jet}");
    style_graph(gr_marginal_flat, gr_marginal_param, "Marginal Bias vs #eta_{#Lambda}; #eta_{#Lambda}; <#mathcal{R}>");

    // Canvas 1: R vs eta_t
    TCanvas* c1 = new TCanvas("c1", "R vs eta_t", 900, 600);
    c1->SetGrid();
    gr_flat->Draw("AL"); 
    gr_param->Draw("L SAME"); 
    adjust_axes(c1, gr_flat, gr_param);
    TLegend* leg1 = new TLegend(0.15, 0.75, 0.45, 0.88);
    leg1->AddEntry(gr_flat, "Flat Phase Space", "l");
    leg1->AddEntry(gr_param, "Realistic Weights", "l");
    leg1->Draw();

    // Canvas 2: R/eta_t vs eta_t
    TCanvas* c2 = new TCanvas("c2", "R/eta_t vs eta_t", 900, 600);
    c2->SetGrid();
    gr_flat_div->Draw("AL"); 
    gr_param_div->Draw("L SAME"); 
    adjust_axes(c2, gr_flat_div, gr_param_div);
    TLegend* leg2 = new TLegend(0.15, 0.75, 0.45, 0.88);
    leg2->AddEntry(gr_flat, "Flat Phase Space", "l");
    leg2->AddEntry(gr_param, "Realistic Weights", "l");
    leg2->Draw();

    // Canvas 3: Marginal R vs eta_Lambda
    TCanvas* c3 = new TCanvas("c3", "Marginal R vs eta_L", 900, 600);
    c3->SetGrid();
    gr_marginal_flat->Draw("AL"); 
    gr_marginal_param->Draw("L SAME"); 
    adjust_axes(c3, gr_marginal_flat, gr_marginal_param);
    TLegend* leg3 = new TLegend(0.55, 0.75, 0.88, 0.88);
    leg3->AddEntry(gr_flat, "Flat Phase Space", "l");
    leg3->AddEntry(gr_param, "Realistic Weights", "l");
    leg3->Draw();
}