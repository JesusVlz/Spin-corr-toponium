// Based on main372.cc (Sjostrand & Preuss, Pythia8 TopThreshold examples).
// Generates parton-level LHE samples (no ISR/FSR/MPI/hadronization) for
// two configurations of the tt-bar threshold model, following the exchange
// with Torbjorn Sjostrand (July 2026):
//
//   model 1: "simple Coulomb enhancement" (no top-width correction needed,
//             the double-counting issue does not apply here).
//   model 2: "Green's function over the whole threshold region", with the
//             top width and the Green's function width set separately TO
//             AVOID DOUBLE-COUNTING:
//                 gammat      = 1.24   (= 1.34 - gammatGreen)
//                 gammatGreen = 0.10
//             (values suggested explicitly by Torbjorn, 22-Jul-2026;
//              NOT derived here, just applied as given).
//
// IMPORTANT (read before using):
//    because with PartonLevel:all = off both coincide, and &pythia.process
//    is the standard Pythia8 convention for dumping the hard process to
//    LHEF.
//  - Check your installed Pythia8 version (>= 8.317 required for the
//    gammat/gammatGreen split). If you have 8.318, ask Torbjorn whether
//    this manual workaround is still needed or whether model=3 already
//    handles it.
//
// Usage:
//   ./main_topThreshold_lhe <model>   with <model> = 1 or 2
//
// Compile (example, adjust paths to your Pythia8 installation):
//   g++ main_topThreshold_lhe.cc -o main_topThreshold_lhe \
//       -I$PYTHIA8/include -O2 -std=c++17 \
//       -L$PYTHIA8/lib -Wl,-rpath,$PYTHIA8/lib -lpythia8 -ldl

#include "Pythia8/Pythia.h"
using namespace Pythia8;

int main(int argc, char* argv[]) {

  // ------------------------------------------------------------------
  // Command-line argument: which model to run.
  // ------------------------------------------------------------------
  if (argc != 2) {
    cout << "Usage: " << argv[0] << " <model>   (model = 1 or 2)" << endl;
    return 1;
  }
  int topModel = atoi(argv[1]);
  if (topModel != 1 && topModel != 2) {
    cout << "Error: model must be 1 (Coulomb) or 2 (Green's function)."
         << endl;
    return 1;
  }

  // ------------------------------------------------------------------
  // Common configuration.
  // ------------------------------------------------------------------
  int    nEvent          = 500000;   // TUNE: see statistics note below
  double mt               = 172.5;
  double eCM              = 13000.;
  bool   thresholdOnly    = true;    // TUNE: see mHat window note below
  int    alphasOrder      = 2;
  double alphasValue      = 0.118;
  double ggSingletFrac    = 2./7.;
  double qqSingletFrac    = 0.;
  double thresholdRegion  = 10.;

  // Top width and Green's function width, depending on the model.
  // tWidthGreen is set unconditionally in both cases (as in main372.cc);
  // for model 1 it has no practical effect since gammatGreen == gammat.
  double gammat, gammatGreen;
  if (topModel == 1) {
    // Simple Coulomb: no double-counting issue, nominal width.
    gammat      = 1.34;
    gammatGreen = 1.34;
  } else {
    // Full Green's function, width split as suggested by Torbjorn
    // (22-Jul-2026) to avoid double-counting of the top width.
    gammat      = 1.24;   // = 1.34 - gammatGreen
    gammatGreen = 0.10;
  }

  string outFile = "events_topThreshold_model" + std::to_string(topModel)
                  + ".lhe";

  // ------------------------------------------------------------------
  // Generator.
  // ------------------------------------------------------------------
  Pythia pythia;

  // Process: gg -> ttbar, qqbar -> ttbar; t -> b W -> b (e/mu) nu.
  pythia.readString("Top:gg2ttbar = on");
  pythia.readString("Top:qqbar2ttbar = on");
  pythia.readString("6:onMode = off");
  pythia.readString("6:onIfAll = 24 5");
  pythia.readString("24:onMode = off");
  pythia.readString("24:onIfAny = 11 13");

  // Kinematics.
  pythia.settings.parm("Beams:eCM", eCM);

  // Restrict to threshold region.
  // thresholdOnly = true  : narrow window around 2*mt (300-400 GeV), matches
  //                         the setup used for the fig372 comparison plots.
  // thresholdOnly = false : open window from 200 GeV upwards (no explicit
  //                         mHatMax, i.e. up to eCM), to keep the full
  //                         above-threshold spectrum for a sanity check
  //                         against the narrow-window sample. TUNE both
  //                         edges to whatever region you actually need.
  if (thresholdOnly) {
    pythia.readString("PhaseSpace:mHatMin = 300.");
    pythia.readString("PhaseSpace:mHatMax = 400.");
  }
  else pythia.readString("PhaseSpace:mHatMin = 200.");

  // Threshold model.
  pythia.settings.mode("TopThreshold:model", topModel);
  pythia.particleData.m0(6, mt);
  pythia.readString("6:doForceWidth = true");
  pythia.particleData.mWidth(6, gammat);
  pythia.settings.parm("TopThreshold:tWidthGreen", gammatGreen);
  pythia.settings.parm("TopThreshold:thrRegion", thresholdRegion);
  pythia.settings.mode("TopThreshold:alphasOrder", alphasOrder);
  pythia.settings.parm("TopThreshold:alphasValue", alphasValue);
  pythia.settings.parm("TopThreshold:ggSingletFrac", ggSingletFrac);
  pythia.settings.parm("TopThreshold:qqSingletFrac", qqSingletFrac);

  // Pure parton level: no ISR/FSR/MPI/hadronization.
  pythia.readString("PartonLevel:ISR = off");
  pythia.readString("PartonLevel:FSR = off");
  pythia.readString("PartonLevel:MPI = off");
  pythia.readString("HadronLevel:all = off");
  pythia.readString("Next:numberCount = 20000");

  if (!pythia.init()) {
    cout << "Error: Pythia failed to initialize (model=" << topModel << ")."
         << endl;
    return 1;
  }

  // ------------------------------------------------------------------
  // LHEF writing from Pythia's hard-process record.
  // ------------------------------------------------------------------
  LHAupFromPYTHIA8 myLHA(&pythia.process, &pythia.info);
  myLHA.openLHEF(outFile.c_str());
  myLHA.setInit();
  myLHA.initLHEF();

  int nAccepted = 0;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) continue;
    myLHA.setEvent();
    myLHA.eventLHEF();
    ++nAccepted;
  }

  pythia.stat();
  myLHA.updateSigma();
  myLHA.closeLHEF(true);

  cout << "\nmodel = " << topModel
       << "  gammat = " << gammat
       << "  gammatGreen = " << gammatGreen
       << "  accepted events = " << nAccepted << " / " << nEvent
       << "\nLHE file: " << outFile << endl;

  return 0;
}