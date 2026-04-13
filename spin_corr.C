/*
 * TTbar Spin Correlation Analysis for Delphes
 * Converted from Rivet CMS_2016_I1413748
 *
 * Modified version:
 * - Removed h_dphidressedleptons (particle-level)
 * - Removed h_lep_costheta_CPV
 * - Added h_cos_opening_angle_lab (opening angle in lab frame before boost)
 *
 * Usage:
 *   root -l examples/spin_corr.C'("ttbar.root","output.root")'
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

#include "TFile.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include "TSystem.h"

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
    cout << "TTbar Spin Correlation Analysis" << endl;
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

    TH1D* h_c1c2 = new TH1D("h_c1c2",
        "c_{1}c_{2};c_{1}c_{2};1/#sigma d#sigma/d(c_{1}c_{2})",
        nbins_c1c2, bins_c1c2);

    TH1D* h_lep_costheta_plus = new TH1D("h_lep_costheta_plus",
        "cos(#theta_{l});cos(#theta_{l});1/#sigma d#sigma/dcos(#theta_{l})",
        nbins_costheta, bins_costheta);
    TH1D* h_lep_costheta_minus = new TH1D("h_lep_costheta_minus",
        "cos(#theta_{l});cos(#theta_{l});1/#sigma d#sigma/dcos(#theta_{l})",
        nbins_costheta, bins_costheta);

    // Set sum of weights squared for proper normalization
    h_dphi->Sumw2();
    h_cos_opening_angle_lab->Sumw2();
    h_cos_opening_angle->Sumw2();
    h_c1c2->Sumw2();
    h_lep_costheta_plus->Sumw2();
    h_lep_costheta_minus->Sumw2();


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
                    // For t → W+ → l+ (both positive)
                    // For tbar → W- → l- (both negative)
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

                // Calculate spin observables in top rest frames
                TVector3 lepPlus_vec = lepPlus_topCM.Vect();
                TVector3 lepMinus_vec = lepMinus_topCM.Vect();
                TVector3 topPlus_vec = topPlus_p4.Vect();
                TVector3 topMinus_vec = topMinus_p4.Vect();

                double lepPlus_costheta = lepPlus_vec.Dot(topPlus_vec) /
                                         (lepPlus_vec.Mag() * topPlus_vec.Mag());

                double lepMinus_costheta = lepMinus_vec.Dot(topMinus_vec) /
                                          (lepMinus_vec.Mag() * topMinus_vec.Mag());

                double c1c2 = lepPlus_costheta * lepMinus_costheta;

                // Opening angle in top rest frames (in ttbar CM frame)
                double cos_opening_angle = lepPlus_vec.Dot(lepMinus_vec) /
                                          (lepPlus_vec.Mag() * lepMinus_vec.Mag());

                // Fill histograms
                fillWithUFOF(h_dphi, dphi_parton);
                fillWithUFOF(h_cos_opening_angle_lab, cos_opening_angle_lab);
                fillWithUFOF(h_cos_opening_angle, cos_opening_angle);
                fillWithUFOF(h_c1c2, c1c2);
                fillWithUFOF(h_lep_costheta_plus, lepPlus_costheta);
                fillWithUFOF(h_lep_costheta_minus, lepMinus_costheta);

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

    // Normalize histograms to unit area
    if (h_dphi->Integral() > 0) {
        h_dphi->Scale(1.0 / h_dphi->Integral(), "width");
    }
    if (h_cos_opening_angle_lab->Integral() > 0) {
        h_cos_opening_angle_lab->Scale(1.0 / h_cos_opening_angle_lab->Integral(), "width");
    }
    if (h_cos_opening_angle->Integral() > 0) {
        h_cos_opening_angle->Scale(1.0 / h_cos_opening_angle->Integral(), "width");
    }
    if (h_c1c2->Integral() > 0) {
        h_c1c2->Scale(1.0 / h_c1c2->Integral(), "width");
    }
    if (h_lep_costheta_plus->Integral() > 0) {
        h_lep_costheta_plus->Scale(1.0 / h_lep_costheta_plus->Integral(), "width");
    }
    if (h_lep_costheta_minus->Integral() > 0) {
        h_lep_costheta_minus->Scale(1.0 / h_lep_costheta_minus->Integral(), "width");
    }

    // Save output
    TFile* outFile = new TFile(outputFile, "RECREATE");
    h_dphi->Write();
    h_cos_opening_angle_lab->Write();
    h_cos_opening_angle->Write();
    h_c1c2->Write();
    h_lep_costheta_plus->Write();
    h_lep_costheta_minus->Write();


    outFile->Close();

    cout << "========================================" << endl;
    cout << "Analysis complete!" << endl;
    cout << "Output saved to: " << outputFile << endl;
    cout << "========================================" << endl;
    cout << endl;

    // Print final statistics
    cout << "=== Final Histogram Statistics ===" << endl;
    cout << "h_dphi entries: " << h_dphi->GetEntries() << endl;
    cout << "h_cos_opening_angle_lab entries: " << h_cos_opening_angle_lab->GetEntries() << endl;
    cout << "h_cos_opening_angle entries: " << h_cos_opening_angle->GetEntries() << endl;
    cout << "h_c1c2 entries: " << h_c1c2->GetEntries() << endl;
    cout << "h_lep_costheta_plus entries: " << h_lep_costheta_plus->GetEntries() << endl;
    cout << "h_lep_costheta_minus entries: " << h_lep_costheta_minus->GetEntries() << endl;

    cout << endl;
}
