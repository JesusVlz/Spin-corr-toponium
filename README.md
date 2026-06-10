# Spin Correlations in tt̄ & Toponium

Analysis to study top–antitop spin-correlation observables near threshold and the effect of **toponium** formation, by comparing Monte Carlo samples against public **CMS** data.

The flow is: generated events (LHE)→ Delphes → `spin_corr.C` (fills histograms) → two ROOT files → `plotting.ipynb` (overlays MC vs CMS data with ratio panels).

## Repository structure

```
.
├── analysis/spin_corr.C        # ROOT/Delphes macro: event loop → spin-correlation histograms
├── plotting.py        # PyROOT: MC-vs-data comparison plots
├── figures/             # output PDFs
├── madgraph/             # sample card generation
└── refs/              # reference papers (see below)
```

## Components

### `spin_corr.C`
ROOT macro adapted from the Rivet routine `CMS_2016_I1413748`. Loops over the Delphes `Particle` (truth) branch, selects dileptonic tt̄ events (two opposite-sign prompt leptons, photon-conversion and multi-lepton vetoes), boosts leptons into the tt̄ CM and parent-top rest frames, and fills six unit-normalized histograms:

- `h_dphi` — Δφ(ℓ⁺, ℓ⁻)
- `h_cos_opening_angle_lab` — cos of the lepton opening angle in the lab frame
- `h_cos_opening_angle` — same, in the top rest frames
- `h_c1c2` — cosθ₁·cosθ₂
- `h_lep_costheta_plus` / `h_lep_costheta_minus` — cosθℓ per lepton

Run:
```bash
root -l 'spin_corr.C("ttbar.root","spin_corr_ttbar.root")'
root -l 'spin_corr.C("toponium.root","spin_corr_topponium.root")'
```

### `plotting.ipynb`
PyROOT notebook that overlays the two MC samples (toponium signal in blue, tt̄ in red) on CMS data points and draws a Data/MC ratio panel. Data are read from HEPData ROOT files (record `ins1742786`). Edit the paths and the `histogram_configs` list at the top, then run all cells; PDFs land in `plots/`.

## Requirements
- ROOT (with the Delphes shared library `libDelphes` reachable)
- Delphes (to produce the input ROOT files)
- Python 3 + PyROOT, Jupyter
- HEPData ROOT export for the CMS measurement (`HEPData-ins1742786-*.root`)

## Reference papers
- **1907.03729** — CMS, *Top polarization and tt̄ spin correlations in dilepton final states at 13 TeV* (PRD 100, 072002). Source of the data and the measured observables.
- **2407.20330** — Aguilar-Saavedra, *Toponium Hunter's Guide*. Strategy + observables for a near-threshold excess.
- **2602.23426** — Antozzi et al., *Extracting a Toponium Signal with Spin and Quantum Information Tools*.
