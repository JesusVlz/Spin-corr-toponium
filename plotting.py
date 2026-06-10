#!/usr/bin/env python3
import ROOT
import os
import argparse
import math

def setup_style():
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetLegendBorderSize(0)
    ROOT.gStyle.SetLegendFillColor(0)
    ROOT.gStyle.SetLegendFont(42)
    ROOT.gStyle.SetLabelFont(42, "XYZ")
    ROOT.gStyle.SetLabelSize(0.04, "XYZ")
    ROOT.gStyle.SetTitleFont(42, "XYZ")
    ROOT.gStyle.SetTitleSize(0.05, "XYZ")
    ROOT.gStyle.SetTickLength(0.02, "XYZ")
    ROOT.gStyle.SetNdivisions(510, "XYZ")
    ROOT.TH1.AddDirectory(False)

def add_text_labels(pad, title="", info="#bf{Parton Level}"):
    pad.cd()
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextFont(42)
    latex.SetTextSize(0.05)
    latex.DrawLatex(0.18, 0.92, title)
    latex.SetTextAlign(31) 
    latex.DrawLatex(0.345, 0.92, info)

def create_ratio_graph(h_mc, h_data_graph):
    ratio_graph = ROOT.TGraphErrors()
    n_points = h_data_graph.GetN()
    point_idx = 0
    for i in range(n_points):
        x_data = h_data_graph.GetPointX(i)
        y_data = h_data_graph.GetPointY(i)
        err_data = (h_data_graph.GetErrorYhigh(i) + h_data_graph.GetErrorYlow(i)) / 2.0
        bin_mc = h_mc.FindBin(x_data)
        y_mc = h_mc.GetBinContent(bin_mc)
        err_mc = h_mc.GetBinError(bin_mc)
        if y_mc > 0 and y_data > 0:
            ratio = y_data / y_mc
            err_ratio = ratio * math.sqrt(pow(err_data/y_data, 2) + pow(err_mc/y_mc, 2))
            ratio_graph.SetPoint(point_idx, x_data, ratio)
            ratio_graph.SetPointError(point_idx, 0, err_ratio)
            point_idx += 1
    return ratio_graph

# COMMENT OUT mc4_file IN THE ARGUMENTS IF USING ONLY 3 MCs. ALSO COMMENT OUT mc3_file IF USING ONLY 2.
def plot_comparison(mc1_file, mc2_file, mc3_file, config, labels, output_dir, y_limits, ratio_limits):
    mc_name, data_file_name, data_path, xtitle, ytitle = config
    
    h1 = mc1_file.Get(mc_name)
    h2 = mc2_file.Get(mc_name)
    h3 = mc3_file.Get(mc_name) # COMMENT OUT IF USING ONLY 2 MCs
    #h4 = mc4_file.Get(mc_name) # COMMENT OUT IF USING ONLY 2 OR 3 MCs
    
    f_data = ROOT.TFile.Open(data_file_name)
    h_data_graph = f_data.Get(data_path).Clone()
    f_data.Close()

    n_bins = h1.GetNbinsX()
    x_min = h1.GetXaxis().GetXmin()
    x_max = h1.GetXaxis().GetXmax()

    c = ROOT.TCanvas(f"c_{mc_name}", "", 800, 800)
    left_margin = 0.16
    
    pad1 = ROOT.TPad("p1", "", 0, 0.3, 1, 1)
    pad1.SetBottomMargin(0.025)
    pad1.SetLeftMargin(left_margin)
    pad1.SetRightMargin(0.05)
    pad1.Draw()
    
    pad2 = ROOT.TPad("p2", "", 0, 0, 1, 0.3)
    pad2.SetTopMargin(0.02)
    pad2.SetBottomMargin(0.35)
    pad2.SetLeftMargin(left_margin)
    pad2.SetRightMargin(0.05)
    pad2.Draw()

    # --- Top Pad ---
    pad1.cd()
    h_frame = ROOT.TH1F(f"frame_top_{mc_name}", "", n_bins, x_min, x_max)
    h_frame.SetMinimum(y_limits[0])
    h_frame.SetMaximum(y_limits[1])
    h_frame.GetYaxis().SetTitle(ytitle)
    h_frame.GetYaxis().SetTitleSize(0.055)
    h_frame.GetYaxis().SetLabelSize(0.055)
    h_frame.GetYaxis().SetTitleOffset(1.4)
    h_frame.GetXaxis().SetLabelSize(0)
    h_frame.Draw("AXIS")

    h1.SetLineColor(ROOT.kGreen); h1.SetLineWidth(2)
    h2.SetLineColor(ROOT.kBlue); h2.SetLineWidth(2); h2.SetLineStyle(2);
    
    # MC3 CONFIGURATION (COMMENT OUT IF USING ONLY 2 MCs)
    h3.SetLineColor(ROOT.kRed); h3.SetLineWidth(2) 
    
    # MC4 CONFIGURATION (COMMENT OUT IF USING ONLY 2 OR 3 MCs)
    #h4.SetLineColor(ROOT.kMagenta); h4.SetLineWidth(2) 

    h1.Draw("HIST E SAME")
    h2.Draw("HIST E SAME")
    h3.Draw("HIST E SAME") # COMMENT OUT IF USING ONLY 2 MCs
    #h4.Draw("HIST E SAME") # COMMENT OUT IF USING ONLY 2 OR 3 MCs
    
    h_data_graph.SetMarkerStyle(20); h_data_graph.SetLineColor(ROOT.kBlack); h_data_graph.SetMarkerSize(1.2)
    h_data_graph.Draw("PE SAME")

    # Legend (adjusted to fit 4 MCs + Data)
    leg = ROOT.TLegend(0.22, 0.68, 0.45, 0.88) 
    leg.AddEntry(h_data_graph, labels['data'], "ep")
    leg.AddEntry(h1, f"{labels['mc1']} ", "l")
    leg.AddEntry(h2, f"{labels['mc2']} ", "l")
    leg.AddEntry(h3, f"{labels['mc3']} ", "l") # COMMENT OUT IF USING ONLY 2 MCs
    #leg.AddEntry(h4, f"{labels['mc4']} (#mu={h4.GetMean():.2f})", "l") # COMMENT OUT IF USING ONLY 2 OR 3 MCs
    leg.Draw()
    add_text_labels(pad1)

    # --- Bottom Pad (Ratio) ---
    pad2.cd()
    h_ratio_frame = ROOT.TH1F(f"frame_ratio_{mc_name}", "", n_bins, x_min, x_max)
    h_ratio_frame.SetMinimum(ratio_limits[0])
    h_ratio_frame.SetMaximum(ratio_limits[1])
    h_ratio_frame.GetYaxis().SetTitle("Data/MC")
    h_ratio_frame.GetYaxis().CenterTitle()
    h_ratio_frame.GetYaxis().SetTitleSize(0.12)
    h_ratio_frame.GetYaxis().SetTitleOffset(0.6)
    h_ratio_frame.GetYaxis().SetLabelSize(0.12)
    h_ratio_frame.GetYaxis().SetNdivisions(505)
    h_ratio_frame.GetXaxis().SetTitle(xtitle)
    h_ratio_frame.GetXaxis().SetTitleSize(0.14)
    h_ratio_frame.GetXaxis().SetTitleOffset(1.1)
    h_ratio_frame.GetXaxis().SetLabelSize(0.12)
    h_ratio_frame.Draw("AXIS")

    ratio1 = create_ratio_graph(h1, h_data_graph)
    ratio2 = create_ratio_graph(h2, h_data_graph)
    ratio3 = create_ratio_graph(h3, h_data_graph) # COMMENT OUT IF USING ONLY 2 MCs
    #ratio4 = create_ratio_graph(h4, h_data_graph) # COMMENT OUT IF USING ONLY 2 OR 3 MCs

    ratio1.SetMarkerColor(ROOT.kGreen); ratio1.SetLineColor(ROOT.kGreen); ratio1.SetMarkerStyle(25); ratio1.SetMarkerSize(1.5) 
    ratio2.SetMarkerColor(ROOT.kBlue); ratio2.SetLineColor(ROOT.kBlue); ratio2.SetMarkerStyle(26); ratio2.SetMarkerSize(1.5)
    
    # MC3 RATIO (COMMENT OUT IF USING ONLY 2 MCs)
    ratio3.SetMarkerColor(ROOT.kRed); ratio3.SetLineColor(ROOT.kRed); ratio3.SetMarkerStyle(24); ratio3.SetMarkerSize(1.5)

    # MC4 RATIO (COMMENT OUT IF USING ONLY 2 OR 3 MCs)
    #ratio4.SetMarkerColor(ROOT.kMagenta); ratio4.SetLineColor(ROOT.kMagenta); ratio4.SetMarkerStyle(23); ratio4.SetMarkerSize(1.5) # Filled diamond

    ratio1.Draw("PE SAME")
    ratio2.Draw("PE E SAME")
    ratio3.Draw("PE E SAME") # COMMENT OUT IF USING ONLY 2 MCs
    #ratio4.Draw("HIST E SAME") # COMMENT OUT IF USING ONLY 2 OR 3 MCs
    
    line = ROOT.TLine(x_min, 1, x_max, 1)
    line.SetLineColor(ROOT.kBlack); line.SetLineStyle(2); line.Draw()

    c.SaveAs(os.path.join(output_dir, f"{mc_name}_data.pdf"))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--f1', default='/home/jesusavc/phd/topponium/outputs/spin_corr_topponium_res_500k.root')
    parser.add_argument('--f2', default='/home/jesusavc/phd/topponium/outputs/spin_corr_topponium_full_500k.root')
    parser.add_argument('--f3', default='/home/jesusavc/phd/topponium/outputs/spin_corr_ttbar_500k.root') # COMMENT OUT IF USING ONLY 2 MCs
   # parser.add_argument('--f4', default='/home/jesusavc/hep-software/Delphes-3.5.0/spin_corr_ttbar_madspin_test.root') # COMMENT OUT IF USING ONLY 2 OR 3 MCs
    parser.add_argument('--out', default='figures')
    args = parser.parse_args()

    setup_style()
    if not os.path.exists(args.out): os.makedirs(args.out)

    f1 = ROOT.TFile.Open(args.f1)
    f2 = ROOT.TFile.Open(args.f2)
    f3 = ROOT.TFile.Open(args.f3) # COMMENT OUT IF USING ONLY 2 MCs
    #f4 = ROOT.TFile.Open(args.f4) # COMMENT OUT IF USING ONLY 2 OR 3 MCs

    y_limit_map = {
        "h_cos_opening_angle": (0.0, 1.0),
        "h_cos_opening_angle_lab": (0.0, 1.8),
        "h_dphi": (0.18, 0.7),
        "h_c1c2": (0.0, 1.5),
        "h_lep_costheta_minus": (0.44, 0.59),
        "h_lep_costheta_plus": (0.44, 0.59)
    }
    
    ratio_limit_map = {
        "h_cos_opening_angle": (0.4, 3.8), 
        "h_cos_opening_angle_lab": (0.4, 3.7), 
        "h_dphi": (0.4, 3.7),
        "h_c1c2": (0.45,3.0),
        "h_lep_costheta_minus": (0.90, 1.09),
        "h_lep_costheta_plus": (0.90, 1.09)
    }

    configs = [
        ("h_cos_opening_angle", "data/HEPData-ins1742786-v1-cosvarphi.root", "cosvarphi/Graph1D_y1", "cos#phi", "1/#sigma d#sigma/dcos#phi"),
        ("h_cos_opening_angle_lab", "data/HEPData-ins1742786-v1-cosvarphi_{mathrm{lab}}.root", "cosvarphi_{mathrm{lab}}/Graph1D_y1", "cos#phi_{lab}", "1/#sigma d#sigma/dcos#phi_{lab}"),
        ("h_dphi", "data/HEPData-ins1742786-v1-_Deltaphi_{ellell}_.root", "|Deltaphi_{ellell}|/Graph1D_y1", "|#Delta#phi_{ll}|", "1/#sigma d#sigma/d#Delta#phi"),
        ("h_c_kk", "data/HEPData-ins1742786-v1-costheta_{1}^{k}costheta_{2}^{k}.root", "costheta_{1}^{k}costheta_{2}^{k}/Graph1D_y1", "cos#theta_{1}^{k}cos#theta_{2}^{k}", "1/#sigma d#sigma/dcos#theta_{1}^{k}cos#theta_{2}^{k}"),
        ("h_lep_costheta_minus", "data/HEPData-ins1742786-v1-costheta_{2}^{k}.root", "costheta_{2}^{k}/Graph1D_y1", "cos(#theta_{2}^{k})", "1/#sigma d#sigma/dcos(#theta_{2}^{k})"),
        ("h_lep_costheta_plus", "data/HEPData-ins1742786-v1-costheta_{1}^{k}.root", "costheta_{1}^{k}/Graph1D_y1", "cos(#theta_{1}^{k})", "1/#sigma d#sigma/dcos(#theta_{1}^{k})")
    ]

    # Added labels for MC3 and MC4 (COMMENT OUT THE ONES YOU DON'T USE)
    labels = {
        'mc1': '#bf{Simplified #eta_{t}-restrict}', 
        'mc2': '#bf{Simplified #eta_{t}-full}', 
        'mc3': '#bf{pp #rightarrow t#bar{t} (SM)}',      # COMMENT OUT IF USING ONLY 2 MCs
        'data': '#bf{CMS Data}'
    }

    for config in configs:
        name = config[0]
        y_lim = y_limit_map.get(name, (0.0, 1.8))
        r_lim = ratio_limit_map.get(name, (0.4, 1.6))
        
        # REMOVE f4 FROM HERE IF USING ONLY 3. REMOVE f3 AND f4 IF USING ONLY 2.
        plot_comparison(f1, f2, f3, config, labels, args.out, y_lim, r_lim) 

    f1.Close(); f2.Close()
    f3.Close() # COMMENT OUT IF USING ONLY 2 MCs
    #f4.Close() # COMMENT OUT IF USING ONLY 2 OR 3 MCs

if __name__ == "__main__":
    main()