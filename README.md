# Animal Coordination Model Simulations

**Companion code for:**
> Cox, C., Ravignani, A., Sørensen, S., Zamm, A., Jensen, F., Bryant, G., Cristia, A., & Fusaroli, R.
> *Breaking the Taxonomic Barrier: A Cross-Species Systematic Review of Formal Models of Coordination*
> Preprint, PsyArxiv.

---

## What is this?

This repository contains all simulation and visualisation code for a systematic review of how animals coordinate their vocalisations in time — from insect choruses to primate duets. The central question is: **what computational mechanism underlies temporal coordination, and can we tell them apart from data?**

We formalise three mechanistically distinct model families, simulate their behaviour, and show that they leave distinguishable fingerprints in inter-call interval distributions and phase dynamics. We then use Approximate Bayesian Computation with Random Forests (ABC-RF) to perform likelihood-free model inference and parameter estimation.

---

## The Three Models

| Model | Core Idea | Signature |
|---|---|---|
| **Event-Triggered** | Discrete, stimulus-driven resets to an endogenous oscillator | Bimodal ICI distribution, oscillatory phase jitter |
| **Tempo-Tracking** | Continuous mutual adjustment of calling frequencies (Kuramoto-style) | Unimodal ICI, smooth phase convergence toward π |
| **State-Dependent** | Probabilistic gating across internal states (Active / Listening / Refractory) | Exponential ICI, structureless phase scatter |

All three models can produce superficially similar alternation patterns. The point of this work is to show they are formally and empirically distinguishable.

---

## Repository Structure

```
AnimalCoordinationModelSimulations/
│
├── 01_SimulationVisualisationCode.Rmd   # Main analysis file (start here)
│
├── SimulationFunctions/
│   ├── 01_eventtriggered.R              # Model 1: Inhibitory reset oscillator
│   ├── 02_tempotracking.R              # Model 2: Kuramoto coupled oscillator
│   └── 03_statedependent.R             # Model 3: Hidden Markov / state-dependent
│
├── CustomFunctions/
│   └── functions.R                      # Plotting, summary statistics, ABC helpers
│
├── ref_tables/                          # CSV checkpoints from ABC reference table
│   ├── ref_event.csv
│   ├── ref_tempo.csv
│   └── ref_state.csv
│
└── tax_full.csv                         # Taxonomy heatmap data
```

---

## Getting Started

### Requirements

All code is written in R. You will need the following packages:

```r
install.packages(c(
  "ggplot2", "dplyr", "cowplot", "patchwork",
  "abcrf", "diptest", "latex2exp",
  "here", "tidyr", "pacman"
))
```

### Running the code

Open `01_SimulationVisualisationCode.Rmd` in RStudio and run the chunks in order. The file is structured as a reproducible R Markdown document covering:

1. Simulating the three coordination models
2. Generating all diagnostic figures (timelines, ICI distributions, phase dynamics)
3. Building the ABC reference table
4. Training the random forest classifier
5. Performing model selection and parameter recovery
6. Generating all manuscript figures

> **Note on the reference table:** Building the reference table (500 simulations per model × 1000 seconds each) takes some time. Results are automatically checkpointed to `ref_tables/` as CSV files so the loop can be safely interrupted and resumed.

---

## Key Parameters

### Event-Triggered (`simulate_event_triggered`)
| Parameter | Description | Default |
|---|---|---|
| `T1`, `T2` | Natural period of each individual (seconds) | 1.0 |
| `rho` | Rebound ratio: how much the period shortens after reset | 0.75 |
| `d_call` | Call duration / inhibition delay (seconds) | 0.2 |
| `noise_sd` | Gaussian noise on each period | 0.2 |

### Tempo-Tracking (`simulate_coupled_oscillator`)
| Parameter | Description | Default |
|---|---|---|
| `freq1_init`, `freq2_init` | Initial calling frequencies (Hz) | 0.8, 1.2 |
| `K` | Phase coupling strength (negative = anti-phase) | −0.1 |
| `epsilon` | Tempo adaptation rate | 0.05 |
| `phase_noise_sd` | Noise in phase dynamics per timestep | 0.08 |

### State-Dependent (`simulate_state_dependent`)
| Parameter | Description | Default |
|---|---|---|
| `lambda_active` | Calling rate when Active (Hz) | 1.0 |
| `p_active_to_listening_partner` | Prob. of suppression when partner calls | 0.3 |
| `refrac_dur` | Minimum refractory period (seconds) | 0.3 |
| `partner_window` | Duration of partner influence (seconds) | 0.5 |

---

## The ABC Pipeline

The model inference approach uses **Approximate Bayesian Computation with Random Forests** ([Pudlo et al. 2016](https://doi.org/10.1093/bioinformatics/btv684)):

1. **Simulate** a reference table of ~1500 datasets per model by sampling random parameters from broad uniform priors
2. **Compute** 8 summary statistics per dataset capturing ICI shape and phase dynamics
3. **Train** a random forest to classify which model generated each dataset
4. **Apply** to observed data to get posterior model probabilities
5. **Rejection ABC** — find the closest reference simulations to observed data; their parameters form the approximate posterior

The 8 summary statistics are:

| Statistic | What it captures |
|---|---|
| `mean_ici` | Average inter-call interval |
| `sd_ici` | Variability in inter-call intervals |
| `bimod` | Bimodality (skewness-kurtosis index) |
| `prop_short` | Proportion of short intervals |
| `dip` | Hartigan's dip statistic for multimodality |
| `phase_dist` | Mean phase distance from anti-phase (π) |
| `phase_var` | Variance in phase differences |
| `phase_acf1` | Lag-1 autocorrelation of phase differences |

---

## Figures

The main analysis file reproduces all figures in the manuscript:

- **Figure 1** — Vocalization timelines, ICI distributions, and phase dynamics for each model
- **Figure 2** — Cross-taxonomic heatmap of coordination models across species
- **Figure 3** — ABC model inference: mechanistic cards, summary statistic scatterplots, variable importance, and cross-model posterior predictive fit

---

## Citation

If you use this code, please cite:

```
Cox, C., Ravignani, A., Sørensen, S., Zamm, A., Jensen, F., Cristia, A., & Fusaroli, R. (in press).
Breaking the Taxonomic Barrier: A Cross-Species Systematic Review of Formal Models of Coordination.
PsyArxiv.
```

---

## Contact

Christopher Cox — [chris.mm.cox@gmail.com](mailto:chris.mm.cox@gmail.com)
Department of Linguistics & Cognitive Science / Interacting Minds Centre, Aarhus University