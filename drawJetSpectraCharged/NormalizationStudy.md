# Jet Cross-Section Normalization Factor Study

## Goal

Measure the inclusive charged jet cross-section:

```
d²σ           N_jet^{INEL}(pT)
────── = ──────────────────────
dpT dη      L_INEL × ΔpT × Δη
```

We measure `N_jet^{sel}` (unfolded, at selected-event level). We need to correct it to INEL level:

```
N_jet^{INEL}(pT) = N_jet^{sel}(pT) / ε_sel^{jet}(pT)
```

**Key principle**: corrections are applied to **jet counts**, not to luminosity. This avoids multiplicity bias — `ε_sel^{jet}(pT)` is measured at the jet level in MC, naturally including correlations between jet production and event selection (e.g., high-pT jets almost always trigger TVX, low-pT jets don't).

---

## Two Independent Approaches

### Approach 1: Data-Driven (BC-level luminosity)

Luminosity is computed from BC-level counters. Corrections are embedded in L itself.

```
σ(pT) = N_jet^{sel}(pT) / L_correct / ε_zvtx / Δη / ΔpT
```

Here L_correct already accounts for TVX, sel8 BC cuts, and RCT. Only ε_zvtx is applied explicitly. **PV reconstruction efficiency is NOT corrected** — this is a known missing piece (~5%).

### Approach 2: MC-Driven (jet-level efficiency correction)

Correct jet counts directly using MC-measured efficiency, pT bin by bin:

```
σ(pT) = N_jet^{sel}(pT) / ε_sel^{jet}(pT) × σ_INEL / N_INEL / Δη / ΔpT
```

`ε_sel^{jet}(pT)` captures ALL selection effects at the jet level, including PV reco, TVX trigger, quality cuts, z-vertex — and their correlations with jet pT.

---

## Approach 1: Data-Driven — Detail

### Flow

```
                          BC-level (data counters)
                          ─────────────────────────

    All BCs in run
         │
         ▼  TVX trigger
    hCounterTVX ──────────────────────── hLumiTVX
         │
         ▼  BC quality (TFBorder + ITSROFBorder)
    hCounterTVXafterBCcuts ──────────── hLumiTVXafterBCcuts
         │                                    │
         │            ε_sel8^{BC} = 0.873     │
         │                                    │
         ▼  RCT quality (run/TF flag)         │
    (no counter)  ───────────────────── hLumiTVXafterBCcutsRCT
                                        (TH2D, Y = CBT_hadronPID)
                    ε_RCT^{BC} = 0.375
                          │
                          ▼  σ_vis rescaling
                    L_correct = L_RCT × (σ_vis^{CCDB} / σ_vis^{correct})
                          │
                          │
                  Collision-level (data)
                  ─────────────────────
                          │
                          ▼  z-vertex acceptance
                    ε_zvtx = 0.949
                    (Gaussian fit on h_collisions_zvertex)
                          │
                          ▼
              ┌─────────────────────────────────┐
              │  σ = N_jet / L_correct / ε_zvtx │
              │      / Δη / ΔpT                 │
              └─────────────────────────────────┘

    External inputs: σ_vis = 53.4 mb (van der Meer scan, 2023)
```

### What's included in L_correct
- TVX trigger efficiency → via σ_vis (σ_vis = σ_INEL × ε_TVX)
- BC-level sel8 (TFBorder, ITSROFBorder) → via hLumiTVXafterBCcuts
- RCT quality selection → via hLumiTVXafterBCcutsRCT
- Pile-up correction → per-BC from CCDB (built into hLumi histograms)

### What's NOT included
- **PV reconstruction efficiency** (~5% jet loss): L is BC-level, N_jet is collision-level
- Correlation between event selection and jet production (assumed negligible at BC level)

### Values (LHC23, R=0.4, data file 608785)
| Variable | Value | Source |
|----------|-------|--------|
| L_TVX | 120440 μb⁻¹ | hLumiTVX |
| L_sel8 | 105137 μb⁻¹ | hLumiTVXafterBCcuts |
| L_RCT | 39412 μb⁻¹ | hLumiTVXafterBCcutsRCT (Y=5) |
| σ_vis_CCDB_est | 58.97 mb | hCounterTVX / hLumiTVX |
| σ_vis_correct | 53.4 mb | van der Meer scan |
| rescale | 1.1043 | σ_vis_CCDB / σ_vis_correct |
| L_correct | 43523 μb⁻¹ | L_RCT × rescale |
| ε_zvtx | 0.9494 | Gaussian fit |
| **normFactor** | **2.420 × 10⁻⁸** | 1 / L_correct[mb] / ε_zvtx |

---

## Approach 2: MC-Driven — Detail

### Flow

```
                   MC simulation (event-level)
                   ────────────────────────────

    INEL collisions (all generated)
    ┌─────────────────────────────────────────────┐
    │ h_mccollisions_eventselection bin 1: N_INEL  │
    │ h2_jet_pt_part_eventselection bin 1: jets    │
    └─────────────────────────────────────────────┘
         │
         │              ε^{evt}     ε^{jet}(20-30 GeV)
         │              ───────     ──────────────────
         ▼  Has reco collision (PV)
         │               0.748          0.951
         ▼  TVX trigger on associated BC
         │               0.898          0.987
         ▼  TF border cut
         │               0.975          0.973
         ▼  ITS ROF border cut
         │               0.845          0.844
         ▼  z-vertex |z| < 10 cm
         │               0.806          1.000 (*)
         │
    ┌─────────────────────────────────────────────┐
    │ ε_total^{evt} = 0.443                        │
    │ ε_total^{jet} = 0.757                        │
    │                                              │
    │ (*) ε_zvtx^{jet} = 1 because jet finder     │
    │     already applies z-vtx cut                │
    └─────────────────────────────────────────────┘
         │
         ▼  Correct jet counts
    ┌──────────────────────────────────────────────────┐
    │                                                  │
    │  N_jet^{INEL}(pT) = N_jet^{sel}(pT)             │
    │                      ─────────────────           │
    │                       ε_sel^{jet}(pT)            │
    │                                                  │
    │  σ(pT) = N_jet^{INEL}(pT) × σ_INEL              │
    │           ──────────────────────────              │
    │            N_INEL × Δη × ΔpT                     │
    │                                                  │
    └──────────────────────────────────────────────────┘

    External inputs: σ_INEL = 78.58 mb (PYTHIA8 Monash)
```

### Why jet-level correction, not event-level?

**Event-level correction introduces multiplicity bias:**
- Events failing TVX have low multiplicity → fewer jets
- Assuming "missing events had same jet rate" overestimates N_jet^{INEL}
- Especially severe at low pT where soft events contribute

**Jet-level correction avoids this:**
- ε_sel^{jet}(pT) directly measures: "of all INEL jets at this pT, what fraction survives selection?"
- Naturally includes correlation: high-pT jets almost always trigger TVX (ε ≈ 1), low-pT jets don't (ε < 1)
- pT-dependent: no single global factor, each bin gets its own correction

### ε_sel^{jet}(pT) — pT dependence

From `h2_jet_pt_part_eventselection` (605134, sel8):

| pT range | ε_sel^{jet} | Comment |
|----------|-------------|---------|
| 5-10 GeV | ~0.72 | Soft events often miss PV reco + TVX |
| 10-20 GeV | ~0.75 | |
| 20-30 GeV | 0.757 | Reference value |
| 30-50 GeV | ~0.76 | |
| 50-100 GeV | ~0.76-0.78 | Almost all trigger TVX |
| >100 GeV | ~0.78 | Near-plateau |

### Values (LHC23k4h, 605134 sel8)
| Variable | Value | Source |
|----------|-------|--------|
| N_INEL | 5.011 × 10⁸ | h_mccollisions_eventselection bin 1 |
| ε_sel^{evt} | 0.4432 | bin 7 / bin 1 |
| ε_sel^{jet}([20,30]) | 0.7570 | h2 bin 7 / bin 1 |
| N_sel_data | 2.078 × 10⁹ | data h_collisions bin 3.5 |
| σ_INEL | 78.58 mb | PYTHIA8 Monash |
| **normFactor** | **2.214 × 10⁻⁸** | σ_INEL × ε_sel^{evt} / (N_sel × ε_sel^{jet}) |

---

## Comparison

```
┌──────────────────────────────────────────────────────────┐
│  normFactor_data / normFactor_MC = 1.093 (9.3%)          │
│                                                          │
│  Breakdown:                                              │
│  ┌────────────────────────────────────────────────────┐  │
│  │ PV reconstruction efficiency        ~5.0%          │  │
│  │   Data-driven: not corrected (BC-level L)          │  │
│  │   MC-driven: included in ε_sel^{jet}               │  │
│  │                                                    │  │
│  │ Pile-up correction model            ~2.5%          │  │
│  │   Data-driven: per-BC CCDB values                  │  │
│  │   MC-driven: no pile-up in MC                      │  │
│  │                                                    │  │
│  │ MC period mismatch                  ~1-2%          │  │
│  │   MC: LHC23k4h (605134)                            │  │
│  │   Data: LHC23 (608785)                             │  │
│  └────────────────────────────────────────────────────┘  │
│                                                          │
│  Action items:                                           │
│  1. Apply PV reco correction to data-driven method       │
│     σ_corrected = σ_measured / ε_PVreco^{jet} ≈ ×1.051  │
│  2. After correction: expect ~4% remaining difference    │
│  3. Use MC-driven as systematic cross-check              │
└──────────────────────────────────────────────────────────┘
```

---

## Validation: MC Truth Cross-Section

Three independent PYTHIA8 Monash calculations at [20,30] GeV/c, R=0.4:

| Method | σ [mb/GeV] | vs HY |
|--------|-----------|-------|
| MC truth INEL (605134 + ε_zvtx correction) | 3.976 × 10⁻³ | 0.992 |
| HY standalone (601994, 4.4B MB events) | 4.008 × 10⁻³ | 1.000 |
| KIAF MB standalone (10M events) | 3.924 × 10⁻³ | 0.979 |

**All agree within 2%** — confirms the MC jet-level efficiency approach is correct.

Rebinned comparison across full pT range (MCtruth/HY ratio):
- pT 20-70 GeV: **0.99-1.01** (< 1% agreement)
- pT > 70 GeV: statistics-limited

---

## Summary Table

| Aspect | Data-Driven | MC-Driven |
|--------|-------------|-----------|
| **Concept** | L from BC counters | ε from MC step-by-step |
| **Correction target** | Luminosity | Jet counts |
| **PV reco efficiency** | ❌ Missing | ✅ Included |
| **Pile-up** | Per-BC CCDB | N/A (no pile-up in MC) |
| **pT dependence** | None (global L) | ✅ ε_sel^{jet}(pT) |
| **External input** | σ_vis (53.4 mb) | σ_INEL (78.58 mb) |
| **Multiplicity bias** | Possible (BC vs coll.) | ✅ Avoided (jet-level) |
| **normFactor** | 2.420 × 10⁻⁸ | 2.214 × 10⁻⁸ |
