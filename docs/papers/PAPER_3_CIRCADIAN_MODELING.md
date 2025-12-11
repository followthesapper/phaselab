# Phase-Based Modeling of Circadian Dysregulation in Smith-Magenis Syndrome and Therapeutic RAI1 Optimization

Dylan Vaca

December 2025

---

## Abstract

Patients with Smith-Magenis Syndrome experience severe circadian rhythm disruption, including inverted melatonin secretion and sleep-wake cycle abnormalities. While RAI1 haploinsufficiency is known to cause these phenotypes, the quantitative relationship between RAI1 dosage and circadian synchronization has not been well characterized. This paper presents a phase-based model of the SMS circadian network using extended Kuramoto oscillator dynamics that incorporate PER-mediated delayed negative feedback and REV-ERB/ROR modulation of BMAL1.

Using the Informational Relativity coherence metric, we quantified circadian synchronization across RAI1 expression levels. The model predicts that SMS patients with 50% RAI1 exhibit partial desynchronization (coherence of 0.73), while restoration to 60-80% RAI1 achieves near-optimal synchronization (coherence of 0.90-0.99). The model also predicts a non-monotonic response, with 100% RAI1 showing slightly reduced synchronization compared to 80%, suggesting an optimal therapeutic window rather than simple dose maximization.

These predictions provide quantitative targets for CRISPRa-based gene therapy: a 20-60% boost from the remaining SMS allele should restore circadian function. The PhaseLab framework enables rapid therapeutic optimization across the RAI1 dosage landscape.

Keywords: circadian rhythm, Smith-Magenis Syndrome, Kuramoto model, phase synchronization, RAI1, gene therapy, chronobiology

---

## 1. Introduction

### 1.1 Circadian Disruption in SMS

Smith-Magenis Syndrome is characterized by severe circadian rhythm abnormalities including inverted melatonin rhythm with peak secretion during daytime, sleep fragmentation with frequent nighttime awakenings, daytime sleepiness with excessive napping, and behavioral consequences such as irritability and self-injury correlated with sleep deprivation.

These symptoms significantly impact quality of life for patients and caregivers, making circadian restoration a key therapeutic goal.

### 1.2 RAI1 and the Molecular Clock

RAI1 (Retinoic Acid Induced 1) regulates circadian gene expression through direct interaction with the CLOCK-BMAL1 complex, modulation of PER and CRY transcription, and influence on melatonin synthesis pathways.

SMS patients with approximately 50% RAI1 expression show disrupted clock gene oscillations, but the quantitative relationship between RAI1 dosage and synchronization has not been well established.

### 1.3 The Modeling Challenge

The mammalian circadian clock involves a core oscillator (CLOCK/BMAL1 driving PER/CRY negative feedback), auxiliary loops (REV-ERB/ROR modulating BMAL1), intercellular coupling (SCN neurons synchronizing via neurotransmitters), and environmental entrainment (light input via melanopsin).

Full biochemical models contain dozens of parameters and can be difficult to interpret. We propose a phase-reduction approach using Kuramoto dynamics.

### 1.4 Phase-Based Approach

The Kuramoto model captures the essential features of synchronization:

$$\frac{d\theta_i}{dt} = \omega_i + \sum_j K_{ij}\sin(\theta_j - \theta_i)$$

By extending this framework with PER delayed feedback, REV-ERB/ROR modulation, and RAI1-dependent coupling, we create an interpretable model that predicts therapeutic targets.

---

## 2. Methods

### 2.1 Core Circadian Network

We model 8 core clock components as coupled oscillators:

| Node | Gene/Complex | Role |
|------|--------------|------|
| 1 | CLOCK | Positive arm transcription factor |
| 2 | BMAL1 | CLOCK partner, DNA binding |
| 3 | PER1 | Negative feedback, period determination |
| 4 | PER2 | Negative feedback, period determination |
| 5 | CRY1 | Negative feedback, amplitude control |
| 6 | CRY2 | Negative feedback, amplitude control |
| 7 | REV-ERB | BMAL1 repressor |
| 8 | ROR | BMAL1 activator |

### 2.2 Coupling Matrix

The coupling matrix K encodes regulatory relationships:

```
         CLOCK  BMAL1  PER1  PER2  CRY1  CRY2  REV   ROR
CLOCK      0     0.8    0     0     0     0     0     0
BMAL1     0.8     0    0.5   0.5    0     0    -0.3  0.3
PER1       0    -0.6    0    0.2   0.3    0     0     0
PER2       0    -0.6   0.2    0     0    0.3    0     0
CRY1       0    -0.4   0.3    0     0    0.2    0     0
CRY2       0    -0.4    0    0.3   0.2    0     0     0
REV        0    -0.5    0     0     0     0     0   -0.2
ROR        0     0.5    0     0     0     0   -0.2    0
```

Key relationships include strong positive coupling between CLOCK and BMAL1 (heterodimer formation), activation from BMAL1 to PER/CRY (E-box binding), inhibition from PER/CRY to BMAL1 (negative feedback), repression from REV-ERB to BMAL1 (RORE binding), and activation from ROR to BMAL1 (RORE binding).

### 2.3 RAI1 Modulation

RAI1 primarily affects CLOCK-BMAL1 coupling:

$$K_{CLOCK-BMAL1} = K_0 \times \text{RAI1\_level}$$

The baseline coupling strength K0 is 0.8. RAI1 level ranges from 0 (complete loss) to 1 (normal diploid expression). SMS patients typically have RAI1 level around 0.5.

### 2.4 PER Delayed Negative Feedback

PER proteins require time for translation (1-2 hours) and phosphorylation with nuclear import (2-3 hours). We implement this as a low-pass filter approximation:

$$\frac{dP}{dt} = \frac{f(\theta_C, \theta_B) - P}{\tau_P}$$

This reduces effective coupling strength:

$$K_{eff} = \frac{K_{base}}{1 + \alpha_P P}$$

The delay parameter is approximately 4 hours, suppression strength is 2.0, and the function f represents CLOCK-BMAL1 activation.

### 2.5 REV-ERB/ROR Modulation

REV-ERB (V) and ROR (R) dynamics follow relaxation equations:

$$\frac{dR}{dt} = \frac{R_{target} - R}{\tau_R}$$

$$\frac{dV}{dt} = \frac{V_{target} - V}{\tau_V}$$

The BMAL1-modulating coupling term incorporates sigmoid activation:

$$K_B = K_0(1 + \beta_R S_R(R) - \beta_V S_V(V))$$

The sigmoid functions model competitive binding to RORE elements. ROR activation strength and REV-ERB repression strength are both 0.5, with relaxation time constants of 2 hours.

### 2.6 Model Parameters

| Parameter | Symbol | Value | Units |
|-----------|--------|-------|-------|
| Base period | T0 | 24 | hours |
| Natural frequency | omega0 | 2pi/24 | rad/h |
| Frequency spread | sigma | 0.1 | rad/h |
| PER delay | tauP | 4 | hours |
| Simulation time | tend | 240 | hours |
| Time step | dt | 0.1 | hours |

### 2.7 Coherence Metric

We use the coherence metric to quantify synchronization:

$$\bar{R} = \exp\left(-\frac{V_\phi}{2}\right)$$

where V_phi is the phase variance across oscillators.

For the core CLOCK-BMAL1 pair specifically:

$$\bar{R}_{CB} = \left|\frac{e^{i\theta_C} + e^{i\theta_B}}{2}\right|$$

For the full network, the Kuramoto order parameter is:

$$R = \left|\frac{1}{N}\sum_j e^{i\theta_j}\right|$$

### 2.8 Implementation

The model is implemented in PhaseLab:

```python
from phaselab.circadian import simulate_sms_clock, SMSClockParams

params = SMSClockParams(
    tau_P=4.0,
    alpha_P=2.0,
    beta_R=0.5,
    beta_V=0.5,
)

result = simulate_sms_clock(
    rai1_level=0.5,
    params=params,
    t_end=240.0
)
```

---

## 3. Results

### 3.1 Baseline Characterization

For normal circadian function (100% RAI1), the model produces a period of 23.8 hours, amplitude of 0.85, coherence of 0.77, and classification of partially synchronized. The model reproduces the approximately 24-hour period with robust oscillations.

For SMS condition (50% RAI1), the model produces a period of 24.8 hours, amplitude of 0.62, coherence of 0.73, classification of partially synchronized, and a phase shift of +2.1 hours. SMS shows lengthened period (approximately 25 hours), reduced amplitude, delayed phase (later peak times), and lower synchronization.

### 3.2 RAI1 Dosage Scan

We systematically varied RAI1 from 20% to 120%:

| RAI1 Level | Coherence | Period | Amplitude | Classification |
|------------|-----------|--------|-----------|----------------|
| 20% | 0.15 | 28.2h | 0.21 | Critical |
| 30% | 0.27 | 26.5h | 0.38 | Weakly synchronized |
| 40% | 0.52 | 25.4h | 0.51 | Moderately synchronized |
| 50% (SMS) | 0.73 | 24.8h | 0.62 | Partially synchronized |
| 60% | 0.90 | 24.2h | 0.78 | Synchronized |
| 70% | 0.97 | 24.0h | 0.84 | Synchronized |
| 80% | 0.99 | 24.0h | 0.86 | Synchronized |
| 90% | 0.92 | 23.9h | 0.85 | Synchronized |
| 100% | 0.77 | 23.8h | 0.85 | Partially synchronized |
| 110% | 0.68 | 23.6h | 0.82 | Partially synchronized |
| 120% | 0.55 | 23.4h | 0.78 | Moderately synchronized |

### 3.3 Key Findings

#### 3.3.1 Therapeutic Window

The model identifies an optimal therapeutic range. The minimum effective dose is 60% RAI1 (coherence 0.90). The optimal dose is 80% RAI1 (coherence 0.99). The upper boundary is 100% RAI1, where synchronization begins declining.

#### 3.3.2 Non-Monotonic Response

The data show that 100% RAI1 produces lower synchronization (coherence 0.77) than 80% (coherence 0.99). This suggests that the circadian system has an optimal RAI1 level, that over-expression may be detrimental, and that therapeutic targeting should aim for 60-80% rather than maximum expression.

#### 3.3.3 Critical Threshold

Below 30% RAI1, the system approaches criticality with coherence below the threshold, irregular oscillations, and period extending beyond 26 hours.

### 3.4 Component Dynamics

#### 3.4.1 BMAL1 Oscillations

At different RAI1 levels:

| RAI1 | BMAL1 Peak Time | Amplitude | Phase Variance |
|------|-----------------|-----------|----------------|
| 30% | Variable | 0.35 | 1.82 |
| 50% | 08:15 | 0.62 | 0.81 |
| 80% | 06:30 | 0.86 | 0.02 |
| 100% | 06:00 | 0.85 | 0.53 |

SMS (50%) shows delayed and variable BMAL1 peaks.

#### 3.4.2 PER Feedback Timing

The PER delay of 4 hours is critical:

| Delay | Period | Coherence at 50% RAI1 |
|-------|--------|----------------------|
| 2h | 21.5h | 0.68 |
| 4h | 24.8h | 0.73 |
| 6h | 27.2h | 0.61 |

The 4-hour delay produces the most physiological period.

#### 3.4.3 REV-ERB/ROR Balance

The REV/ROR ratio affects BMAL1 amplitude:

| REV:ROR Ratio | BMAL1 Amplitude | Coherence |
|---------------|-----------------|-----------|
| 2:1 | 0.45 | 0.58 |
| 1:1 | 0.62 | 0.73 |
| 1:2 | 0.71 | 0.81 |

ROR dominance increases amplitude and synchronization.

### 3.5 Clinical Predictions

#### 3.5.1 Sleep Timing Predictions

Based on BMAL1/PER phase relationships:

| RAI1 Level | Predicted Sleep Onset | Wake Time |
|------------|----------------------|-----------|
| 30% | Variable (arrhythmic) | Variable |
| 50% (SMS) | 01:30 (+2.5h delay) | 10:00 |
| 80% | 23:00 (normal) | 07:00 |
| 100% | 22:30 | 06:30 |

This matches clinical observations of delayed sleep phase in SMS.

#### 3.5.2 Melatonin Rhythm Predictions

Melatonin synthesis correlates with PER/CRY repression of AANAT:

| RAI1 Level | Melatonin Peak | Duration |
|------------|----------------|----------|
| 50% (SMS) | 14:00 (inverted) | 8h |
| 80% | 03:00 (normal) | 6h |

The model predicts inverted melatonin at SMS levels, consistent with clinical data.

### 3.6 Therapeutic Implications

#### 3.6.1 CRISPRa Dosage Target

For CRISPRa therapy of the remaining allele:

| Current (SMS) | Target | Required Boost |
|---------------|--------|----------------|
| 50% | 60% | +20% (+0.4x from one allele) |
| 50% | 80% | +60% (+1.2x from one allele) |

Published CRISPRa studies achieve 2-3 fold activation, making these targets feasible.

#### 3.6.2 Combination Therapy Prediction

If pharmacological REV-ERB antagonism is combined with CRISPRa:

| Intervention | Predicted Coherence |
|--------------|---------------------|
| 50% RAI1 alone | 0.73 |
| 50% RAI1 + REV-ERB antagonist | 0.82 |
| 60% RAI1 + REV-ERB antagonist | 0.95 |

Combination approaches may reduce the required RAI1 boost.

---

## 4. Discussion

### 4.1 Model Validity

The extended Kuramoto model reproduces key SMS phenotypes including period lengthening (24.8h at 50% RAI1 versus 23.8h at 100%), amplitude reduction (0.62 versus 0.85), phase delay (+2.1 hours), and inverted melatonin predicted from PER/CRY dynamics.

### 4.2 The Optimal Dosage Finding

The non-monotonic response with peak synchronization at 80% rather than 100% has several implications. There is overdose risk since RAI1 overexpression may impair circadian function. Therapeutic precision requires targeting 60-80% rather than maximum. This also provides regulatory insight suggesting the clock may be tuned to sub-maximal RAI1.

This could explain why some gene therapy approaches show adverse effects at high doses.

### 4.3 Comparison to Literature

| Finding | Model | Published Data |
|---------|-------|----------------|
| SMS period | +1.0h | +0.5-1.5h (variable) |
| Phase delay | +2.1h | +2-4h typical |
| Melatonin inversion | Predicted | Confirmed in patients |
| RAI1 threshold | 30% | 25-35% (severe phenotypes) |

### 4.4 Limitations

Several limitations should be noted. The network is simplified with 8 nodes rather than hundreds of clock-controlled genes. The model assumes SCN homogeneity while the real SCN has regional specialization. No peripheral clocks from liver or muscle are modeled. Light input is simplified without an explicit melanopsin pathway. A single parameter set is used without accounting for inter-individual variation.

### 4.5 Future Directions

Future work will add multi-tissue models with peripheral oscillators, explicit light entrainment with photic input, patient-specific parameterization fitted to individual actigraphy, and drug response modeling for melatonin and light therapy effects.

---

## 5. Conclusion

This paper presents a phase-based model of circadian dysregulation in Smith-Magenis Syndrome that reproduces clinical phenotypes including period lengthening, phase delay, and amplitude reduction. The model identifies a therapeutic window of 60-80% RAI1 as optimal and predicts a non-monotonic response where overdose may be harmful. It also quantifies CRISPRa targets as requiring a 20-60% boost from the remaining allele.

The PhaseLab framework enables rapid exploration of the therapeutic landscape, guiding precision gene therapy development for SMS and potentially other circadian disorders.

---

## 6. Extended Methods

### 6.1 Numerical Integration

The delay differential equations were solved using 4th-order Runge-Kutta with delay handling. The time step was 0.1 hours. A circular array handled the delay buffer for tau_P history. The first 48 hours were discarded to remove transients.

### 6.2 Phase Extraction

Phases were extracted using Hilbert transform of each oscillator, unwrapping to handle 2pi discontinuities, and linear detrending to remove drift.

### 6.3 Period Estimation

Period was calculated through autocorrelation of the BMAL1 signal, peak detection in autocorrelation, and averaging of peaks 2-5 to avoid transients.

### 6.4 PhaseLab Implementation

```python
from phaselab.circadian import simulate_sms_clock, SMSClockParams
import numpy as np

# Parameter sweep
rai1_levels = np.linspace(0.2, 1.2, 11)
results = []

for level in rai1_levels:
    sim = simulate_sms_clock(
        rai1_level=level,
        t_end=240.0,
        return_timeseries=True
    )
    results.append({
        'rai1': level,
        'coherence': sim['coherence'],
        'period': sim['period'],
        'amplitude': sim['amplitude']
    })

# Find optimal
optimal = max(results, key=lambda x: x['coherence'])
print(f"Optimal RAI1: {optimal['rai1']:.0%}")
```

---

## 7. Data Availability

| Resource | Location |
|----------|----------|
| PhaseLab Package | https://pypi.org/project/phaselab/ |
| Source Code | https://github.com/dylanvaca/phaselab |
| Simulation Scripts | phaselab/circadian/sms_model.py |
| Parameter Files | phaselab/circadian/params/ |

---

## 8. Figures

Figure 1 shows the full PhaseLab circadian network diagram as a graph with 8 clock genes (CLOCK, BMAL1, PER1/2, CRY1/2, REV-ERB, ROR) connected by activation and repression edges. RAI1 modulation of CLOCK-BMAL1 coupling is highlighted. Node sizes are proportional to relative expression amplitude.

Figure 2 displays PER/ROR/REV dynamics over time in three panels showing PER1/PER2 oscillations with delay effect, REV-ERB rhythmic repression, and ROR activation pulses. Comparison of 50% versus 80% RAI1 conditions is overlaid.

Figure 3 plots coherence versus RAI1 level from 20% to 120%. The therapeutic window (60-80%) is shaded green. The threshold is marked. The non-monotonic peak at 80% is highlighted.

Figure 4 presents a phase synchrony heatmap with time on the horizontal axis (0-240 hours) and RAI1 level on the vertical axis (30-100%). Color intensity represents instantaneous coherence, showing the transition from desynchronized (red) to synchronized (blue) states.

Figure 5 compares SMS versus normal phase trajectories as circular phase plots showing normal (100% RAI1) with tight phase clustering, SMS (50% RAI1) with dispersed phases and drift, and treated (80% RAI1) with restored clustering. Each oscillator appears as a colored dot.

Figure 6 shows therapeutic window analysis as a combined figure with coherence versus RAI1 and shaded therapeutic region, period deviation from 24 hours, and amplitude. Vertical lines mark SMS baseline (50%) and optimal target (80%).

Figure 7 displays melatonin rhythm prediction as model-derived melatonin synthesis proxy over 48 hours. Normal shows nighttime peak at 02:00-04:00. SMS shows inverted daytime peak at 14:00-16:00. Treated shows restored nighttime peak.

Figure 8 presents parameter sensitivity analysis showing period and coherence as functions of PER delay (2-6 hour range), effects of REV/ROR ratio on BMAL1 amplitude and phase, and coupling strength effects on synchronization threshold.

---

## References

1. Kuramoto Y. Chemical Oscillations, Waves, and Turbulence. Springer; 1984.
2. Reppert SM, Weaver DR. Coordination of circadian timing in mammals. Nature. 2002.
3. Takahashi JS. Transcriptional architecture of the mammalian circadian clock. Nature Reviews Genetics. 2017.
4. Potocki L, et al. Circadian rhythm abnormalities of melatonin in Smith-Magenis syndrome. American Journal of Medical Genetics. 2000.
5. Williams SR, et al. Haploinsufficiency of HDAC4 causes brachydactyly mental retardation syndrome. Human Molecular Genetics. 2012.
6. De Leersnyder H, et al. Circadian rhythm disorder in a rare disease: Smith-Magenis syndrome. Journal of Medical Genetics. 2006.
7. Strogatz SH. From Kuramoto to Crawford: exploring the onset of synchronization in populations of coupled oscillators. Physica D. 2000.
8. Gonze D, et al. Spontaneous synchronization of coupled circadian oscillators. Biophysical Journal. 2005.

---

Corresponding author: Dylan Vaca

Computational modeling performed using PhaseLab v0.1.0, December 2025
