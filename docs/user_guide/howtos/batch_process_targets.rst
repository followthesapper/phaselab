Batch Process Targets
=====================

This guide shows how to efficiently process multiple gene targets.

Basic Batch Processing
----------------------

Process multiple genes:

.. code-block:: python

   from phaselab.crispr import design_guides
   import pandas as pd

   # Multiple target sequences
   targets = {
       'RAI1': {
           'sequence': 'ATGCGATCGATCGATCGATCG...',
           'tss': 200,
       },
       'CLOCK': {
           'sequence': 'GCTAGCTAGCTAGCTAGCTAG...',
           'tss': 150,
       },
       'BMAL1': {
           'sequence': 'TAGCTAGCTAGCTAGCTAGCT...',
           'tss': 180,
       },
   }

   # Process each target
   all_guides = []
   for gene, info in targets.items():
       guides = design_guides(info['sequence'], info['tss'])
       guides['target_gene'] = gene
       all_guides.append(guides)

   # Combine results
   combined = pd.concat(all_guides, ignore_index=True)
   print(f"Total guides designed: {len(combined)}")

   # Best guide per gene
   best_per_gene = combined.groupby('target_gene').apply(
       lambda x: x.nlargest(1, 'combined_score')
   ).reset_index(drop=True)

   print(best_per_gene[['target_gene', 'sequence', 'combined_score']])

Parallel Processing
-------------------

Use multiprocessing for large batches:

.. code-block:: python

   from concurrent.futures import ProcessPoolExecutor, as_completed
   from phaselab.crispr import design_guides
   import pandas as pd

   def process_gene(gene_info):
       """Process a single gene target."""
       gene, info = gene_info
       try:
           guides = design_guides(info['sequence'], info['tss'])
           guides['target_gene'] = gene
           return guides
       except Exception as e:
           print(f"Error processing {gene}: {e}")
           return pd.DataFrame()

   # Process in parallel
   with ProcessPoolExecutor(max_workers=4) as executor:
       futures = {
           executor.submit(process_gene, (gene, info)): gene
           for gene, info in targets.items()
       }

       results = []
       for future in as_completed(futures):
           gene = futures[future]
           try:
               result = future.result()
               results.append(result)
               print(f"Completed: {gene}")
           except Exception as e:
               print(f"Failed: {gene} - {e}")

   combined = pd.concat(results, ignore_index=True)

Batch with Coherence Validation
-------------------------------

Add quantum coherence to batch results:

.. code-block:: python

   from phaselab.crispr import design_guides, compute_coherence_batch
   import pandas as pd

   # Process genes
   all_guides = []
   for gene, info in targets.items():
       guides = design_guides(info['sequence'], info['tss'])
       guides['target_gene'] = gene

       # Filter to GO candidates
       go_guides = guides[guides['go_no_go'] == 'GO'].head(10)
       all_guides.append(go_guides)

   combined = pd.concat(all_guides, ignore_index=True)

   # Batch quantum validation (efficient)
   combined['quantum_coherence'] = compute_coherence_batch(
       combined['sequence'].tolist(),
       mode="quantum"
   )

   # Rank by quantum coherence
   combined = combined.sort_values('quantum_coherence', ascending=False)

Loading from File
-----------------

Process targets from a file:

.. code-block:: python

   import pandas as pd
   from phaselab.crispr import design_guides

   # Load target list
   targets_df = pd.read_csv('targets.csv')
   # Expected columns: gene_name, sequence, tss_position

   results = []
   for idx, row in targets_df.iterrows():
       print(f"Processing {row['gene_name']}...")

       guides = design_guides(
           row['sequence'],
           row['tss_position']
       )
       guides['target_gene'] = row['gene_name']

       # Take top 5 per gene
       top_guides = guides.nlargest(5, 'combined_score')
       results.append(top_guides)

   # Save results
   combined = pd.concat(results, ignore_index=True)
   combined.to_csv('designed_guides.csv', index=False)

Progress Tracking
-----------------

Track progress for long batches:

.. code-block:: python

   from tqdm import tqdm
   from phaselab.crispr import design_guides
   import pandas as pd

   results = []
   for gene, info in tqdm(targets.items(), desc="Designing guides"):
       guides = design_guides(info['sequence'], info['tss'])
       guides['target_gene'] = gene
       results.append(guides[guides['go_no_go'] == 'GO'].head(5))

   combined = pd.concat(results, ignore_index=True)

Error Handling
--------------

Handle failures gracefully:

.. code-block:: python

   from phaselab.crispr import design_guides
   import pandas as pd
   import logging

   logging.basicConfig(level=logging.INFO)
   logger = logging.getLogger(__name__)

   successful = []
   failed = []

   for gene, info in targets.items():
       try:
           # Validate input
           if len(info['sequence']) < 100:
               raise ValueError(f"Sequence too short for {gene}")

           guides = design_guides(info['sequence'], info['tss'])

           if len(guides) == 0:
               raise ValueError(f"No guides found for {gene}")

           guides['target_gene'] = gene
           successful.append(guides)
           logger.info(f"Success: {gene} ({len(guides)} guides)")

       except Exception as e:
           failed.append({'gene': gene, 'error': str(e)})
           logger.error(f"Failed: {gene} - {e}")

   # Report
   print(f"Successful: {len(successful)} genes")
   print(f"Failed: {len(failed)} genes")

   if failed:
       print("Failed genes:")
       for f in failed:
           print(f"  {f['gene']}: {f['error']}")

Export Results
--------------

Export batch results in various formats:

.. code-block:: python

   # CSV export
   combined.to_csv('guides.csv', index=False)

   # Excel with multiple sheets
   with pd.ExcelWriter('guides.xlsx') as writer:
       # Summary sheet
       summary = combined.groupby('target_gene').agg({
           'combined_score': ['count', 'mean', 'max'],
           'quantum_coherence': 'mean'
       })
       summary.to_excel(writer, sheet_name='Summary')

       # Per-gene sheets
       for gene in combined['target_gene'].unique():
           gene_guides = combined[combined['target_gene'] == gene]
           gene_guides.to_excel(writer, sheet_name=gene[:30], index=False)

   # JSON export
   combined.to_json('guides.json', orient='records', indent=2)

See Also
--------

- :doc:`filter_guides_by_coherence` - Filter results by coherence
- :doc:`compare_crispr_modalities` - Compare modality approaches
