Circadian Clock Modeling
========================

This tutorial covers modeling circadian rhythms for SMS gene therapy research.

Smith-Magenis Syndrome
----------------------

Smith-Magenis Syndrome (SMS) is caused by RAI1 haploinsufficiency, leading to:

- Inverted melatonin rhythm
- Shortened circadian period
- Sleep disturbances
- Reduced circadian coherence

PhaseLab models these effects for CRISPRa gene therapy optimization.

Basic Clock Simulation
----------------------

.. code-block:: python

   from phaselab.circadian import simulate_sms_clock

   # Simulate wild-type (100% RAI1)
   wt = simulate_sms_clock(rai1_level=1.0)
   print(f"Wild-type period: {wt['period']:.2f} hours")
   print(f"Wild-type coherence: {wt['coherence']:.4f}")

   # Simulate SMS (50% RAI1)
   sms = simulate_sms_clock(rai1_level=0.5)
   print(f"SMS period: {sms['period']:.2f} hours")
   print(f"SMS coherence: {sms['coherence']:.4f}")
   print(f"Phase shift: {sms['phase_shift']:.2f} hours")

Custom Clock Parameters
-----------------------

.. code-block:: python

   from phaselab.circadian import simulate_sms_clock, SMSClockParams
   import numpy as np

   # Custom parameters
   params = SMSClockParams(
       tau_P=4.0,           # PER delay time constant (hours)
       alpha_P=2.0,         # PER suppression strength
       beta_R=0.5,          # RORa effect on BMAL1
       beta_V=0.5,          # REV-ERBa effect on BMAL1
       K_light=0.1,         # Light sensitivity
       omega_0=2*np.pi/24,  # Base angular frequency
   )

   result = simulate_sms_clock(rai1_level=0.5, params=params)

Time Series Output
------------------

Get full time series for visualization:

.. code-block:: python

   from phaselab.circadian import simulate_sms_clock
   import matplotlib.pyplot as plt

   results = simulate_sms_clock(
       rai1_level=0.5,
       t_end=240.0,           # 10 days
       dt=0.1,                # Time step (hours)
       return_timeseries=True
   )

   # Extract time series
   t = results['time']
   bmal1 = results['BMAL1']
   per = results['PER']
   rev_erb = results['REV_ERB']
   ror = results['ROR']

   # Plot
   fig, axes = plt.subplots(2, 2, figsize=(12, 8))

   axes[0, 0].plot(t, bmal1)
   axes[0, 0].set_title('BMAL1')
   axes[0, 0].set_xlabel('Time (hours)')

   axes[0, 1].plot(t, per)
   axes[0, 1].set_title('PER')

   axes[1, 0].plot(t, rev_erb)
   axes[1, 0].set_title('REV-ERB')

   axes[1, 1].plot(t, ror)
   axes[1, 1].set_title('ROR')

   plt.tight_layout()
   plt.savefig('sms_clock.png')

Therapeutic Dose-Response
-------------------------

Find the RAI1 level needed for GO coherence:

.. code-block:: python

   from phaselab.circadian import simulate_sms_clock
   from phaselab.core.constants import E_MINUS_2
   import numpy as np

   # Scan RAI1 levels
   rai1_levels = np.linspace(0.3, 1.0, 15)
   results = []

   for level in rai1_levels:
       sim = simulate_sms_clock(rai1_level=level)
       results.append({
           'rai1': level,
           'period': sim['period'],
           'coherence': sim['coherence'],
           'go': sim['coherence'] > E_MINUS_2
       })

   # Find therapeutic threshold
   import pandas as pd
   df = pd.DataFrame(results)
   threshold = df[df['go']]['rai1'].min()
   print(f"Minimum RAI1 for GO: {threshold:.1%}")

   # Plot dose-response
   import matplotlib.pyplot as plt

   fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

   ax1.plot(df['rai1'] * 100, df['coherence'], 'bo-')
   ax1.axhline(E_MINUS_2, color='r', linestyle='--', label=f'GO threshold ({E_MINUS_2:.3f})')
   ax1.axvline(threshold * 100, color='g', linestyle=':', label=f'Min RAI1 ({threshold:.0%})')
   ax1.set_xlabel('RAI1 Level (%)')
   ax1.set_ylabel('Coherence R')
   ax1.legend()
   ax1.set_title('Coherence Dose-Response')

   ax2.plot(df['rai1'] * 100, df['period'], 'go-')
   ax2.axhline(24, color='gray', linestyle='--', label='24h')
   ax2.set_xlabel('RAI1 Level (%)')
   ax2.set_ylabel('Period (hours)')
   ax2.legend()
   ax2.set_title('Period Dose-Response')

   plt.tight_layout()
   plt.savefig('dose_response.png')

Kuramoto Oscillators
--------------------

Model coupled oscillator networks:

.. code-block:: python

   from phaselab.circadian import KuramotoNetwork

   # Create network of 100 oscillators
   network = KuramotoNetwork(
       n_oscillators=100,
       coupling_strength=0.5,
       natural_frequency_spread=0.1
   )

   # Simulate
   results = network.simulate(t_end=100.0)

   print(f"Order parameter R: {results['order_parameter']:.4f}")
   print(f"Synchronized: {'Yes' if results['synchronized'] else 'No'}")

Heterogeneous Networks
^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   from phaselab.circadian import KuramotoNetwork
   import numpy as np

   # Bimodal frequency distribution (morning/evening types)
   frequencies = np.concatenate([
       np.random.normal(1.0, 0.05, 50),   # Morning types
       np.random.normal(1.1, 0.05, 50)    # Evening types
   ])

   network = KuramotoNetwork(
       n_oscillators=100,
       coupling_strength=0.8,
       natural_frequencies=frequencies
   )

   # With light-dark forcing
   def light_forcing(t):
       return 0.1 * np.sin(2 * np.pi * t / 24)

   results = network.simulate(
       t_end=240.0,
       external_forcing=light_forcing
   )

Multi-Tissue Modeling
---------------------

Model inter-tissue coupling:

.. code-block:: python

   from phaselab.circadian import MultiTissueModel

   # Create multi-tissue model
   model = MultiTissueModel(
       tissues=['SCN', 'liver', 'muscle', 'adipose'],
       coupling_matrix=[
           [0.0, 0.8, 0.6, 0.5],  # SCN coupling to others
           [0.2, 0.0, 0.3, 0.2],  # Liver
           [0.3, 0.3, 0.0, 0.2],  # Muscle
           [0.2, 0.4, 0.2, 0.0],  # Adipose
       ]
   )

   # Simulate SMS condition
   results = model.simulate(
       rai1_level=0.5,
       t_end=240.0
   )

   for tissue, data in results['tissues'].items():
       print(f"{tissue}: period={data['period']:.2f}h, coherence={data['coherence']:.3f}")

Chronotherapy Optimization
--------------------------

Optimize dosing schedule:

.. code-block:: python

   from phaselab.circadian import optimize_dosing_schedule

   # Optimize CRISPRa dosing
   schedule = optimize_dosing_schedule(
       target_gene='RAI1',
       delivery_method='AAV',
       tissue_targets=['SCN', 'liver'],
       constraints={
           'max_doses_per_day': 1,
           'min_interval_hours': 24,
           'treatment_duration_days': 7
       }
   )

   print(f"Optimal dose times: {schedule['dose_times']}")
   print(f"Predicted final coherence: {schedule['predicted_coherence']:.4f}")

Integration with CRISPR Design
------------------------------

Combine circadian modeling with guide design:

.. code-block:: python

   from phaselab.crispr import design_guides
   from phaselab.circadian import simulate_sms_clock

   # Design RAI1-targeting guides
   rai1_sequence = "..."  # RAI1 promoter
   guides = design_guides(rai1_sequence, tss_index=200)

   # Filter GO candidates
   go_guides = guides[guides['go_no_go'] == 'GO']

   # For each guide, predict therapeutic outcome
   for idx, row in go_guides.head(5).iterrows():
       # Assume guide efficacy correlates with combined_score
       predicted_rai1 = 0.5 + 0.35 * row['combined_score']  # 50% + restoration

       sim = simulate_sms_clock(rai1_level=predicted_rai1)

       print(f"Guide: {row['sequence'][:10]}...")
       print(f"  Predicted RAI1: {predicted_rai1:.1%}")
       print(f"  Predicted period: {sim['period']:.2f}h")
       print(f"  Predicted coherence: {sim['coherence']:.4f} ({sim['status']})")
       print()

Next Steps
----------

- :doc:`crispr_guide_design` - Design guides for RAI1
- :doc:`quantum_validation` - Validate on hardware
