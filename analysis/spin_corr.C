/*
 * TTbar Spin Correlation Analysis for Delphes
 * Converted from Rivet CMS_2016_I1413748
 *
 * Modified version:
 * - Removed h_dphidressedleptons (particle-level)
 * - Removed h_lep_costheta_CPV
 * - Added h_cos_opening_angle_lab (opening angle in lab frame before boost)
 * - Added full k, n, r spin quantization basis (C_kk, C_nn, C_rr)
 * - Corrected n and r axes to match exact CMS definitions (arXiv:1907.03729)
 * - Added individual 1D cosines for n and r axes
 * - Added prints for mean spin correlation coefficients and uncertainties
 *
 * Usage:
 * root -l examples/spin_corr.C'("ttbar.root","output.root")'
 */

#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#else
class ExRootTreeReader;
class ExRootResult;
#endif

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>

#include "TFile.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include "TSystem.h"
#include "TVector3.h"

using namespace std;

// Constants
const double PI = TMath::Pi();

const double DphiMax =  PI;
const double DphiMin = 0.0;

// Histogram binning (same as Rivet analysis)
const int nbins_dphi = 6;
const double bins_dphi[] = {0., 1.*PI/6., 1.*PI/3., 1.*PI/2., 2.*PI/3., 5.*PI/6., PI};

const int nbins_costheta = 6;
const double bins_costheta[] = {-1., -2./3., -1./3., 0., 1./3., 2./3., 1.};

const int nbins_c1c2 = 6;
const double bins_c1c2[] = {-1., -2./3., -1./3., 0., 1./3., 2./3., 1.};

const int nbins_cos_opening = 6;
const double bins_cos_opening[] = {-1., -2./3., -1./3., 0., 1./3., 2./3., 1.};

//------------------------------------------------------------------------------

// Helper function: Check if particle is prompt (not from hadron decay)
bool isPrompt(GenParticle* particle, TClonesArray* branchParticle) {
    if (!particle) return false;

    // Status 1 particles are final state
    if (particle->Status == 1) return true;

    // Check parent
    if (particle->M1 >= 0 && particle->M1 < branchParticle->GetEntriesFast()) {
        GenParticle* mother = (GenParticle*) branchParticle->At(particle->M1);
        int absPID = abs(mother->PID);

        // Not from hadron decay (reject if parent is hadron)
        if (absPID > 100) return false;
    }

    return true;
}

//------------------------------------------------------------------------------

// Helper function: Check if lepton has photon parent (reject leptons from photon conversion)
bool hasPhotonParent(GenParticle* particle, TClonesArray* branchParticle) {
    if (!particle) return false;

    GenParticle* current = particle;
    while (current && current->M1 >= 0 && current->M1 < branchParticle->GetEntriesFast()) {
        GenParticle* mother = (GenParticle*) branchParticle->At(current->M1);
        if (mother->PID == 22) return true; // Found photon parent
        // Stop if we reach a top quark or W boson
        if (abs(mother->PID) == 6 || abs(mother->PID) == 24) break;
        current = mother;
    }
    return false;
}

//------------------------------------------------------------------------------

// Helper function: Check if particle comes from a specific parent (by PID)
bool hasParentWithPID(GenParticle* particle, int parentPID, TClonesArray* branchParticle) {
    if (!particle) return false;

    GenParticle* current = particle;
    while (current && current->M1 >= 0 && current->M1 < branchParticle->GetEntriesFast()) {
        GenParticle* mother = (GenParticle*) branchParticle->At(current->M1);
        if (mother->PID == parentPID) return true;
        current = mother;
    }
    return false;
}

//------------------------------------------------------------------------------

// Helper function: Fill histogram with under/overflow protection
void fillWithUFOF(TH1D* h, double x) {
    double xMin = h->GetXaxis()->GetXmin();
    double xMax = h->GetXaxis()->GetXmax();
    double xFill = max(min(x, xMax - 1e-9), xMin + 1e-9);
    h->Fill(xFill);
}

//------------------------------------------------------------------------------

void spin_corr(const char* inputFile, const char* outputFile = "ttbar_spin_output.root") {

    gSystem->Load("libDelphes");

    cout << "========================================" << endl;
    cout << "TTbar Spin Correlation Analysis (Full k, n, r Basis)" << endl;
    cout << "========================================" << endl;
    cout << "Input file: " << inputFile << endl;
    cout << "Output file: " << outputFile << endl;
    cout << endl;

    // Create chain of root trees
    TChain chain("Delphes");
    chain.Add(inputFile);

    // Create object of class ExRootTreeReader
    ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
    Long64_t numberOfEntries = treeReader->GetEntries();

    cout << "Total events to process: " << numberOfEntries << endl;
    cout << endl;

    // Get pointers to branches used in this analysis
    TClonesArray *branchParticle = treeReader->UseBranch("Particle");

    // Create output histograms
    TH1D* h_dphi = new TH1D("h_dphi",
        "#Delta#phi(l^{+},l^{-}) - Parton Level;#Delta#phi;1/#sigma d#sigma/d#Delta#phi",
        nbins_dphi, bins_dphi);

    TH1D* h_cos_opening_angle_lab = new TH1D("h_cos_opening_angle_lab",
        "cos(#theta_{opening}) - Lab Frame;cos(#theta_{opening});1/#sigma d#sigma/dcos(#theta_{opening})",
        nbins_cos_opening, bins_cos_opening);

    TH1D* h_cos_opening_angle = new TH1D("h_cos_opening_angle",
        "cos(#theta_{opening}) - Top Rest Frame;cos(#theta_{opening});1/#sigma d#sigma/dcos(#theta_{opening})",
        nbins_cos_opening, bins_cos_opening);

    TH1D* h_c_kk = new TH1D("h_c_kk",
        "C_{kk};cos(#theta_{1}^{k})cos(#theta_{2}^{k});1/#sigma d#sigma/d(C_{kk})",
        nbins_c1c2, bins_c1c2);

    TH1D* h_c_nn = new TH1D("h_c_nn",
        "C_{nn};cos(#theta_{1}^{n})cos(#theta_{2}^{n});1/#sigma d#sigma/d(C_{nn})",
        nbins_c1c2, bins_c1c2);

    TH1D* h_c_rr = new TH1D("h_c_rr",
        "C_{rr};cos(#theta_{1}^{r})cos(#theta_{2}^{r});1/#sigma d#sigma/d(C_{rr})",
        nbins_c1c2, bins_c1c2);

    // 1D distributions for k-axis
    TH1D* h_lep_costheta_plus = new TH1D("h_lep_costheta_plus",
        "cos(#theta_{1}^{k});cos(#theta_{1}^{k});1/#sigma d#sigma/dcos(#theta_{1}^{k})",
        nbins_costheta, bins_costheta);
        
    TH1D* h_lep_costheta_minus = new TH1D("h_lep_costheta_minus",
        "cos(#theta_{2}^{k});cos(#theta_{2}^{k});1/#sigma d#sigma/dcos(#theta_{2}^{k})",
        nbins_costheta, bins_costheta);

    // 1D distributions for n-axis
    TH1D* h_lep_costheta_n_plus = new TH1D("h_lep_costheta_n_plus",
        "cos(#theta_{1}^{n});cos(#theta_{1}^{n});1/#sigma d#sigma/dcos(#theta_{1}^{n})",
        nbins_costheta, bins_costheta);
        
    TH1D* h_lep_costheta_n_minus = new TH1D("h_lep_costheta_n_minus",
        "cos(#theta_{2}^{n});cos(#theta_{2}^{n});1/#sigma d#sigma/dcos(#theta_{2}^{n})",
        nbins_costheta, bins_costheta);

    // 1D distributions for r-axis
    TH1D* h_lep_costheta_r_plus = new TH1D("h_lep_costheta_r_plus",
        "cos(#theta_{1}^{r});cos(#theta_{1}^{r});1/#sigma d#sigma/dcos(#theta_{1}^{r})",
        nbins_costheta, bins_costheta);
        
    TH1D* h_lep_costheta_r_minus = new TH1D("h_lep_costheta_r_minus",
        "cos(#theta_{2}^{r});cos(#theta_{2}^{r});1/#sigma d#sigma/dcos(#theta_{2}^{r})",
        nbins_costheta, bins_costheta);

    // Set sum of weights squared for proper normalization and exact mean tracking
    h_dphi->Sumw2();
    h_cos_opening_angle_lab->Sumw2();
    h_cos_opening_angle->Sumw2();
    h_c_kk->Sumw2();
    h_c_nn->Sumw2();
    h_c_rr->Sumw2();
    h_lep_costheta_plus->Sumw2();
    h_lep_costheta_minus->Sumw2();
    h_lep_costheta_n_plus->Sumw2();
    h_lep_costheta_n_minus->Sumw2();
    h_lep_costheta_r_plus->Sumw2();
    h_lep_costheta_r_minus->Sumw2();

    // Counters for debugging
    int nEventsPartonLevel = 0;
    int nRejectedPhotonParent = 0;
    int nRejectedChargeInconsistent = 0;
    int nRejectedMultipleLeptons = 0;

    // Event loop
    for(Int_t entry = 0; entry < numberOfEntries; ++entry) {

        // Load selected branches with data from specified event
        treeReader->ReadEntry(entry);

        if (entry % 1000 == 0) {
            cout << "Processing event " << entry << " / " << numberOfEntries << endl;
        }

        // Find top quarks
        GenParticle* topPlus = nullptr;
        GenParticle* topMinus = nullptr;

        for (int i = 0; i < branchParticle->GetEntriesFast(); ++i) {
            GenParticle* particle = (GenParticle*) branchParticle->At(i);
            int pid = particle->PID;
            int absPID = abs(pid);

            // Find top quarks (last in decay chain before decay)
            if (absPID == 6) {
                // Check if this is the last top before decay
                bool isLastTop = true;
                if (particle->D1 >= 0) {
                    for (int j = particle->D1; j <= particle->D2 && j < branchParticle->GetEntriesFast(); ++j) {
                        GenParticle* daughter = (GenParticle*) branchParticle->At(j);
                        if (abs(daughter->PID) == 6) {
                            isLastTop = false;
                            break;
                        }
                    }
                }

                if (isLastTop) {
                    if (pid > 0) topPlus = particle;
                    else topMinus = particle;
                }
            }
        }

        // ===== PARTON-LEVEL ANALYSIS =====
        if (topPlus && topMinus) {

            vector<GenParticle*> chargedLeptons;
            int nTrueLeptonicTops = 0;

            // Loop over both tops
            vector<GenParticle*> tops = {topPlus, topMinus};

            for (size_t k = 0; k < tops.size(); ++k) {
                GenParticle* lepTop = tops[k];

                // Find all prompt charged lepton candidates from this top
                vector<GenParticle*> lepton_candidates;

                for (int i = 0; i < branchParticle->GetEntriesFast(); ++i) {
                    GenParticle* particle = (GenParticle*) branchParticle->At(i);
                    int absPID = abs(particle->PID);

                    // Look for prompt charged leptons (e or mu)
                    if ((absPID == 11 || absPID == 13)) {

                        // Check if this lepton comes from this specific top
                        if (hasParentWithPID(particle, lepTop->PID, branchParticle)) {
                            lepton_candidates.push_back(particle);
                        }
                    }
                }

                // Process lepton candidates following Rivet logic
                bool isTrueLeptonicTop = false;

                for (size_t i = 0; i < lepton_candidates.size(); ++i) {
                    GenParticle* lepton_candidate = lepton_candidates[i];

                    // Reject leptons from photon conversion (like in Rivet)
                    if (hasPhotonParent(lepton_candidate, branchParticle)) {
                        nRejectedPhotonParent++;
                        continue;
                    }

                    // Check charge consistency: lepton should have same sign as top
                    bool chargeConsistent = (lepton_candidate->Charge * lepTop->PID > 0);

                    if (!isTrueLeptonicTop && chargeConsistent) {
                        chargedLeptons.push_back(lepton_candidate);
                        isTrueLeptonicTop = true;
                    }
                    else if (chargeConsistent) {
                        // Found extra lepton with correct charge - reject event
                        nRejectedMultipleLeptons++;
                    }
                }

                if (isTrueLeptonicTop) nTrueLeptonicTops++;
            }

            // Require exactly 2 true leptonic tops with opposite sign leptons
            bool oppositesign = false;
            if (nTrueLeptonicTops == 2 && chargedLeptons.size() == 2) {
                oppositesign = (chargedLeptons[0]->Charge * chargedLeptons[1]->Charge < 0);
                if (!oppositesign) nRejectedChargeInconsistent++;
            }

            if (nTrueLeptonicTops == 2 && oppositesign) {

                // Build 4-vectors
                TLorentzVector topPlus_p4, topMinus_p4;
                topPlus_p4.SetPxPyPzE(topPlus->Px, topPlus->Py, topPlus->Pz, topPlus->E);
                topMinus_p4.SetPxPyPzE(topMinus->Px, topMinus->Py, topMinus->Pz, topMinus->E);

                TLorentzVector lepPlus_p4, lepMinus_p4;
                if (chargedLeptons[0]->Charge > 0) {
                    lepPlus_p4.SetPxPyPzE(chargedLeptons[0]->Px, chargedLeptons[0]->Py,
                                         chargedLeptons[0]->Pz, chargedLeptons[0]->E);
                    lepMinus_p4.SetPxPyPzE(chargedLeptons[1]->Px, chargedLeptons[1]->Py,
                                          chargedLeptons[1]->Pz, chargedLeptons[1]->E);
                } else {
                    lepPlus_p4.SetPxPyPzE(chargedLeptons[1]->Px, chargedLeptons[1]->Py,
                                         chargedLeptons[1]->Pz, chargedLeptons[1]->E);
                    lepMinus_p4.SetPxPyPzE(chargedLeptons[0]->Px, chargedLeptons[0]->Py,
                                          chargedLeptons[0]->Pz, chargedLeptons[0]->E);
                }

                // Calculate dphi in lab frame (before boost)
                double dphi_parton = fabs(lepPlus_p4.DeltaPhi(lepMinus_p4));

                // Calculate opening angle in LAB FRAME (before boost)
                TVector3 lepPlus_vec_lab = lepPlus_p4.Vect();
                TVector3 lepMinus_vec_lab = lepMinus_p4.Vect();
                double cos_opening_angle_lab = lepPlus_vec_lab.Dot(lepMinus_vec_lab) /
                                               (lepPlus_vec_lab.Mag() * lepMinus_vec_lab.Mag());

                // ttbar system
                TLorentzVector ttbar_p4 = topPlus_p4 + topMinus_p4;

                // Boost to ttbar CM frame
                TVector3 beta_ttbar = ttbar_p4.BoostVector();

                topPlus_p4.Boost(-beta_ttbar);
                topMinus_p4.Boost(-beta_ttbar);
                lepPlus_p4.Boost(-beta_ttbar);
                lepMinus_p4.Boost(-beta_ttbar);

                // Boost leptons to their parent top rest frames
                TVector3 beta_topPlus = topPlus_p4.BoostVector();
                TVector3 beta_topMinus = topMinus_p4.BoostVector();

                TLorentzVector lepPlus_topCM = lepPlus_p4;
                TLorentzVector lepMinus_topCM = lepMinus_p4;

                lepPlus_topCM.Boost(-beta_topPlus);
                lepMinus_topCM.Boost(-beta_topMinus);

                // -------------------------------------------------------------
                // Define the full spin quantization basis (k, n, r)
                // Exactly matching CMS/ATLAS definitions (arXiv:1907.03729)
                // -------------------------------------------------------------
                TVector3 p_p(0, 0, 1); // Direction of one of the initial protons (Z-axis)
                
                TVector3 lepPlus_vec = lepPlus_topCM.Vect().Unit();
                TVector3 lepMinus_vec = lepMinus_topCM.Vect().Unit();

                // === TOP QUARK AXES ===
                TVector3 k_plus = topPlus_p4.Vect().Unit(); 
                
                // Calculate production angle theta_t for top
                double cos_theta_t_plus = p_p.Dot(k_plus);
                double sin_theta_t_plus = sqrt(1.0 - cos_theta_t_plus * cos_theta_t_plus);
                
                // R-axis definition: (p_p - cos(theta_t)*k) / sin(theta_t)
                TVector3 r_plus = (p_p - cos_theta_t_plus * k_plus);
                if (sin_theta_t_plus > 1e-6) {
                    r_plus = r_plus * (1.0 / sin_theta_t_plus);
                }
                r_plus = r_plus.Unit(); // Ensure normalization

                // N-axis definition: k x r
                TVector3 n_plus = k_plus.Cross(r_plus).Unit();

                // Top lepton projections (cosines)
                double cos_theta_k_plus = lepPlus_vec.Dot(k_plus);
                double cos_theta_n_plus = lepPlus_vec.Dot(n_plus);
                double cos_theta_r_plus = lepPlus_vec.Dot(r_plus);

                // === ANTI-TOP QUARK AXES ===
                TVector3 k_minus = topMinus_p4.Vect().Unit(); 
                
                // Calculate production angle theta_t for anti-top
                double cos_theta_t_minus = p_p.Dot(k_minus);
                double sin_theta_t_minus = sqrt(1.0 - cos_theta_t_minus * cos_theta_t_minus);
                
                // R-axis definition: (p_p - cos(theta_t)*k) / sin(theta_t)
                TVector3 r_minus = (p_p - cos_theta_t_minus * k_minus);
                if (sin_theta_t_minus > 1e-6) {
                    r_minus = r_minus * (1.0 / sin_theta_t_minus);
                }
                r_minus = r_minus.Unit();

                // N-axis definition: k x r
                TVector3 n_minus = k_minus.Cross(r_minus).Unit();

                // Anti-top lepton projections (cosines)
                double cos_theta_k_minus = lepMinus_vec.Dot(k_minus);
                double cos_theta_n_minus = lepMinus_vec.Dot(n_minus);
                double cos_theta_r_minus = lepMinus_vec.Dot(r_minus);

                // === SPIN CORRELATIONS (Products) ===
                double c_kk = cos_theta_k_plus * cos_theta_k_minus;
                double c_nn = cos_theta_n_plus * cos_theta_n_minus;
                double c_rr = cos_theta_r_plus * cos_theta_r_minus;

                // Opening angle in top rest frames (in ttbar CM frame)
                double cos_opening_angle = lepPlus_vec.Dot(lepMinus_vec);

                // Fill histograms
                fillWithUFOF(h_dphi, dphi_parton);
                fillWithUFOF(h_cos_opening_angle_lab, cos_opening_angle_lab);
                fillWithUFOF(h_cos_opening_angle, cos_opening_angle);
                
                fillWithUFOF(h_c_kk, c_kk);
                fillWithUFOF(h_c_nn, c_nn);
                fillWithUFOF(h_c_rr, c_rr);
                
                fillWithUFOF(h_lep_costheta_plus, cos_theta_k_plus);
                fillWithUFOF(h_lep_costheta_minus, cos_theta_k_minus);

                fillWithUFOF(h_lep_costheta_n_plus, cos_theta_n_plus);
                fillWithUFOF(h_lep_costheta_n_minus, cos_theta_n_minus);

                fillWithUFOF(h_lep_costheta_r_plus, cos_theta_r_plus);
                fillWithUFOF(h_lep_costheta_r_minus, cos_theta_r_minus);

                nEventsPartonLevel++;
            }
        }

    } // End event loop

    cout << endl;
    cout << "Event processing complete." << endl;
    cout << endl;

    // Print debugging info
    cout << "=== Event Selection Statistics ===" << endl;
    cout << "Parton level events passed: " << nEventsPartonLevel << endl;
    cout << "Events rejected (photon parent): " << nRejectedPhotonParent << endl;
    cout << "Events rejected (charge inconsistent): " << nRejectedChargeInconsistent << endl;
    cout << "Events rejected (multiple leptons): " << nRejectedMultipleLeptons << endl;
    cout << endl;

    // ----------------------------------------------------------------------
    // Print the Spin Correlation Coefficients (Mean) and their uncertainties
    // ----------------------------------------------------------------------
    cout << "==========================================================" << endl;
    cout << "      SPIN CORRELATION COEFFICIENTS & UNCERTAINTIES       " << endl;
    cout << "==========================================================" << endl;
    cout << fixed << setprecision(5);
    cout << "C_kk Mean (helicity axis)   : " << -9*h_c_kk->GetMean() << " +/- " << -9*h_c_kk->GetMeanError() << endl;
    cout << "C_nn Mean (transverse axis) : " << -9*h_c_nn->GetMean() << " +/- " << -9*h_c_nn->GetMeanError() << endl;
    cout << "C_rr Mean (orthogonal axis) : " << -9*h_c_rr->GetMean() << " +/- " << -9*h_c_rr->GetMeanError() << endl;
    cout << "----------------------------------------------------------" << endl;
    cout << "cos(theta_k_plus) Mean      : " << 3*h_lep_costheta_plus->GetMean() << " +/- " << 3*h_lep_costheta_plus->GetMeanError() << endl;
    cout << "cos(theta_k_minus) Mean     : " << 3*h_lep_costheta_minus->GetMean() << " +/- " << 3*h_lep_costheta_minus->GetMeanError() << endl;
    cout << "cos(theta_n_plus) Mean      : " << 3*h_lep_costheta_n_plus->GetMean() << " +/- " << 3*h_lep_costheta_n_plus->GetMeanError() << endl;
    cout << "cos(theta_n_minus) Mean     : " << 3*h_lep_costheta_n_minus->GetMean() << " +/- " << 3*h_lep_costheta_n_minus->GetMeanError() << endl;
    cout << "cos(theta_r_plus) Mean      : " << 3*h_lep_costheta_r_plus->GetMean() << " +/- " << 3*h_lep_costheta_r_plus->GetMeanError() << endl;
    cout << "cos(theta_r_minus) Mean     : " << 3*h_lep_costheta_r_minus->GetMean() << " +/- " << 3*h_lep_costheta_r_minus->GetMeanError() << endl;
    cout << "==========================================================" << endl;
    cout << endl;
 
    // Normalize histograms to unit area
    if (h_dphi->Integral() > 0) h_dphi->Scale(1.0 / h_dphi->Integral(), "width");
    if (h_cos_opening_angle_lab->Integral() > 0) h_cos_opening_angle_lab->Scale(1.0 / h_cos_opening_angle_lab->Integral(), "width");
    if (h_cos_opening_angle->Integral() > 0) h_cos_opening_angle->Scale(1.0 / h_cos_opening_angle->Integral(), "width");
    if (h_c_kk->Integral() > 0) h_c_kk->Scale(1.0 / h_c_kk->Integral(), "width");
    if (h_c_nn->Integral() > 0) h_c_nn->Scale(1.0 / h_c_nn->Integral(), "width");
    if (h_c_rr->Integral() > 0) h_c_rr->Scale(1.0 / h_c_rr->Integral(), "width");
    if (h_lep_costheta_plus->Integral() > 0) h_lep_costheta_plus->Scale(1.0 / h_lep_costheta_plus->Integral(), "width");
    if (h_lep_costheta_minus->Integral() > 0) h_lep_costheta_minus->Scale(1.0 / h_lep_costheta_minus->Integral(), "width");
    if (h_lep_costheta_n_plus->Integral() > 0) h_lep_costheta_n_plus->Scale(1.0 / h_lep_costheta_n_plus->Integral(), "width");
    if (h_lep_costheta_n_minus->Integral() > 0) h_lep_costheta_n_minus->Scale(1.0 / h_lep_costheta_n_minus->Integral(), "width");
    if (h_lep_costheta_r_plus->Integral() > 0) h_lep_costheta_r_plus->Scale(1.0 / h_lep_costheta_r_plus->Integral(), "width");
    if (h_lep_costheta_r_minus->Integral() > 0) h_lep_costheta_r_minus->Scale(1.0 / h_lep_costheta_r_minus->Integral(), "width");

    // Save output
    TFile* outFile = new TFile(outputFile, "RECREATE");
    h_dphi->Write();
    h_cos_opening_angle_lab->Write();
    h_cos_opening_angle->Write();
    h_c_kk->Write();
    h_c_nn->Write();
    h_c_rr->Write();
    h_lep_costheta_plus->Write();
    h_lep_costheta_minus->Write();
    h_lep_costheta_n_plus->Write();
    h_lep_costheta_n_minus->Write();
    h_lep_costheta_r_plus->Write();
    h_lep_costheta_r_minus->Write();

    outFile->Close();

    cout << "========================================" << endl;
    cout << "Analysis complete!" << endl;
    cout << "Output saved to: " << outputFile << endl;
    cout << "========================================" << endl;
}