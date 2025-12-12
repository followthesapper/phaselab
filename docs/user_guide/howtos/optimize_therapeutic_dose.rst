Optimize Therapeutic Dose
=========================

This guide shows how to use circadian modeling to optimize gene therapy dosage.

Finding Therapeutic Threshold
-----------------------------

Determine the minimum RAI1 level for GO coherence:

.. code-block:: python

   from phaselab.circadian import simulate_sms_clock
   from phaselab.core.constants import E_MINUS_2
   import numpy as np

   # Scan RAI1 levels
   rai1_levels = np.linspace(0.3, 1.0, 50)

   results = []
   for level in rai1_levels:
       sim = simulate_sms_clock(rai1_level=level)
       results.append({
           'rai1': level,
           'coherence': sim['coherence'],
           'period': sim['period'],
           'is_go': sim['coherence'] > E_MINUS_2
       })

   # Find threshold
   import pandas as pd
   df = pd.DataFrame(results)

   threshold = df[df['is_go']]['rai1'].min()
   print(f"Minimum RAI1 for GO: {threshold:.1%}")

   # Safety margin (20% above threshold)
   target_rai1 = threshold * 1.2
   print(f"Recommended target: {target_rai1:.1%}")

Dose-Response Curve
-------------------

Generate a dose-response curve:

.. code-block:: python

   import matplotlib.pyplot as plt
   from phaselab.circadian import simulate_sms_clock
   from phaselab.core.constants import E_MINUS_2
   import numpy as np

   rai1_levels = np.linspace(0.3, 1.0, 30)
   coherences = []
   periods = []

   for level in rai1_levels:
       sim = simulate_sms_clock(rai1_level=level)
       coherences.append(sim['coherence'])
       periods.append(sim['period'])

   # Plot
   fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

   ax1.plot(rai1_levels * 100, coherences, 'b-', linewidth=2)
   ax1.axhline(E_MINUS_2, color='r', linestyle='--', label=f'GO threshold')
   ax1.fill_between(rai1_levels * 100, 0, coherences,
                    where=[c > E_MINUS_2 for c in coherences],
                    alpha=0.3, color='green', label='GO region')
   ax1.set_xlabel('RAI1 Level (%)')
   ax1.set_ylabel('Coherence (R)')
   ax1.legend()
   ax1.set_title('Coherence Dose-Response')

   ax2.plot(rai1_levels * 100, periods, 'g-', linewidth=2)
   ax2.axhline(24, color='gray', linestyle='--', label='24h')
   ax2.set_xlabel('RAI1 Level (%)')
   ax2.set_ylabel('Period (hours)')
   ax2.legend()
   ax2.set_title('Period Dose-Response')

   plt.tight_layout()
   plt.savefig('dose_response.png', dpi=150)

Guide Efficacy to RAI1 Model
----------------------------

Estimate RAI1 restoration from guide efficacy:

.. code-block:: python

   from phaselab.crispr import design_guides
   from phaselab.circadian import simulate_sms_clock

   def estimate_rai1_restoration(guide_score: float, baseline: float = 0.5) -> float:
       """
       Estimate RAI1 level from CRISPRa guide efficacy.

       Model: RAI1 = baseline + efficacy * (1 - baseline) * guide_score
       - baseline: SMS condition (50% RAI1)
       - efficacy: CRISPRa maximum restoration (~40%)
       - guide_score: Combined guide score (0-1)
       """
       max_restoration = 0.40  # 40% maximum restoration from CRISPRa
       return baseline + max_restoration * guide_score

   # Design guides
   guides = design_guides(rai1_promoter, tss_index)
   go_guides = guides[guides['go_no_go'] == 'GO'].head(10)

   # Predict therapeutic outcomes
   for idx, row in go_guides.iterrows():
       predicted_rai1 = estimate_rai1_restoration(row['combined_score'])
       sim = simulate_sms_clock(rai1_level=predicted_rai1)

       print(f"Guide: {row['sequence'][:15]}...")
       print(f"  Score: {row['combined_score']:.3f}")
       print(f"  Predicted RAI1: {predicted_rai1:.1%}")
       print(f"  Predicted coherence: {sim['coherence']:.4f}")
       print(f"  Status: {sim['status']}")
       print()

Multi-Tissue Optimization
-------------------------

Optimize for multiple tissues:

.. code-block:: python

   from phaselab.circadian import MultiTissueModel

   # Create multi-tissue model
   model = MultiTissueModel(
       tissues=['SCN', 'liver', 'muscle'],
       coupling_matrix=[
           [0.0, 0.8, 0.6],  # SCN
           [0.2, 0.0, 0.3],  # Liver
           [0.3, 0.3, 0.0],  # Muscle
       ]
   )

   # Find dose that achieves GO in all tissues
   for rai1 in [0.6, 0.7, 0.8, 0.9]:
       results = model.simulate(rai1_level=rai1)

       all_go = all(
           data['coherence'] > E_MINUS_2
           for data in results['tissues'].values()
       )

       print(f"RAI1 {rai1:.0%}: {'All GO' if all_go else 'Some NO-GO'}")
       for tissue, data in results['tissues'].items():
           status = 'GO' if data['coherence'] > E_MINUS_2 else 'NO-GO'
           print(f"  {tissue}: R={data['coherence']:.3f} [{status}]")

Therapeutic Window
------------------

Identify the therapeutic window:

.. code-block:: python

   from phaselab.circadian import simulate_sms_clock
   import numpy as np

   # Scan for therapeutic window
   rai1_levels = np.linspace(0.5, 1.2, 50)  # Include over-expression

   therapeutic_window = []
   for level in rai1_levels:
       sim = simulate_sms_clock(rai1_level=level)

       # Therapeutic criteria:
       # 1. Coherence > GO threshold
       # 2. Period within normal range (23-25h)
       # 3. No over-expression toxicity (< 110%)

       is_therapeutic = (
           sim['coherence'] > E_MINUS_2 and
           23 <= sim['period'] <= 25 and
           level <= 1.1
       )

       therapeutic_window.append({
           'rai1': level,
           'is_therapeutic': is_therapeutic,
           'coherence': sim['coherence'],
           'period': sim['period']
       })

   df = pd.DataFrame(therapeutic_window)
   therapeutic = df[df['is_therapeutic']]

   print(f"Therapeutic window: {therapeutic['rai1'].min():.0%} - {therapeutic['rai1'].max():.0%}")

See Also
--------

- :doc:`/user_guide/tutorials/circadian_modeling` - Full circadian tutorial
- :doc:`filter_guides_by_coherence` - Filter guides by coherence
