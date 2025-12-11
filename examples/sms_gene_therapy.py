#!/usr/bin/env python3
"""
Example: Smith-Magenis Syndrome CRISPRa Guide Design

This example demonstrates the full PhaseLab pipeline for designing
CRISPRa guide RNAs to upregulate RAI1 expression in SMS patients.

Based on E200 experiment from Informational Relativity framework.
"""

from phaselab.crispr import design_guides, GuideDesignConfig
from phaselab.circadian import simulate_sms_clock, therapeutic_scan
from phaselab.io import export_guides_fasta, export_crispor_batch

# RAI1 promoter sequence (chr17:17680958-17681958, GRCh38)
RAI1_PROMOTER = """
TGTCTCTTCCCACCAGGATGCCAGCTCCACGACAGCAGGGACCGCATTGTCTCTGCTGATTCCTCAAGACTCAGCAAGGTGCCAGGTGCCCAGTAAATGAAGTTTACATGAAGCCAATGGAGTTTTGAAATCGAAAATGTTTCTCAAACTGCAAAGAAGTGGGCACAGGTGAGGGCTAATGGTGGGGTGAGGGTGGGGGCCTGGCCTGACGCTGGAAGCTCCTGTACGCAGGGCGCTAGTGTCCTCCCAGGCCGCACCGGGATCCCACCGCAGAGGCCCCCTCCCGCAGGCCGGCCAGCGGGCGGGGGGCGCGCAGTGGGGCGCGCGGACCGCCAGGGGGTGCTGTCGGGGCGCGCACGGCCCACCTGGCGCCCGGGGCCGCGAGCGAGCGCAGCCAGGCCGCCCGGGAACTTTTTTCCTCGCGCTCTTGCTCTCCTTCTCCCCCCGCGCGTCGCCCGCCCGCCCGTCCGCTCGTCCGCTCTTCCTGTGCGGAGGGGATCGCATTCCTCGCGGCACCCGCGGACACCAGGCAGGCCCGGCCGGGTGCGCGGCCAAGGCCCGCTCCCCGCGTCCCCGGGCGCCGCCCCCGCCCGCGGCTGGGCTCCGAGAGACGAGTGGGAGAGCGAGTGCAGCGAGGGCGAGAGGCGAGCCGAGCAGGCCGCCCTGCCCGCGGCCCTAGCGCCGGCGCGAGGAGGGGGCGCCGCGGCCCACCCTCCTTCCTGCCTGGCCGCGGGCCGGCCGGGCCGCCGCGCCCCGACCCCCATGGCCACCCAGGCCTCCGGGCCGCGAAGTCGCAGCGCCAGACCCAAGGCCCCCGAGTGAGCGCGGGCGCCGAGGTGAGCAGCGAGCGCCGGGGCGCGGGGGGCGCAGTGTATCCTGCCCCGGGGCGGGGGCGGCGCGAGGCGTCGGGATCCTAGGCCGGGTGATGAGCCGAGGCGGGTTCGGACGCCGAGCGCGCGCGCTGGGTGGAGATGTGGGGTCCGCGGGGTCCTTGCTTTGTG
""".replace("\n", "").strip()


def main():
    print("=" * 60)
    print("PhaseLab: SMS Gene Therapy CRISPRa Pipeline")
    print("=" * 60)

    # TSS is at position 500 in this 1kb promoter
    tss_index = 500

    # Step 1: Design guide RNAs
    print("\n[1] Designing CRISPRa guide RNAs...")
    config = GuideDesignConfig(
        pam="NGG",
        crispr_window=(-400, -50),
        min_gc=0.4,
        max_gc=0.7,
        compute_coherence=True,
    )

    guides = design_guides(
        sequence=RAI1_PROMOTER,
        tss_index=tss_index,
        config=config,
        verbose=True,
    )

    print(f"\nTop 5 guide candidates:")
    print("-" * 80)
    for i, row in guides.head(5).iterrows():
        print(f"{i+1}. {row['sequence']}")
        print(f"   Position: {row['position']:+d} bp | GC: {row['gc']:.0%} | "
              f"MIT: {row['mit_score']:.0f} | R̄: {row['coherence_R']:.3f} | {row['go_no_go']}")

    # Step 2: Circadian model - find therapeutic window
    print("\n[2] Simulating circadian clock response...")
    scan = therapeutic_scan(
        rai1_levels=[0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
        n_trials=3,
    )

    print("\nRAI1 Level vs Circadian Synchronization:")
    print("-" * 40)
    for level, R_bar, cls in zip(scan['levels'], scan['R_bars'], scan['classifications']):
        bar = "█" * int(R_bar * 20)
        print(f"  {level*100:3.0f}% | {R_bar:.3f} {bar} {cls}")

    print(f"\nSMS baseline (50%): R̄ = {scan['sms_baseline_R']:.3f}")
    print(f"Optimal RAI1 level: {scan['optimal_level']*100:.0f}%")
    if scan['required_boost']:
        print(f"Required boost from SMS: +{scan['required_boost']*100:.0f}%")

    # Step 3: Detailed simulation at therapeutic level
    print("\n[3] Detailed simulation at optimal RAI1 level...")
    result = simulate_sms_clock(
        rai1_level=scan['optimal_level'],
        t_end=240,  # 10 days
    )
    print(f"Final synchronization: R̄ = {result['final_R_bar']:.3f}")
    print(f"Classification: {result['classification']}")
    print(f"GO/NO-GO: {result['go_no_go']}")

    # Step 4: Export for CRISPOR validation
    print("\n[4] Exporting for CRISPOR validation...")
    export_crispor_batch(RAI1_PROMOTER, "rai1_crispor.fa", name="RAI1_promoter")
    export_guides_fasta(guides, "rai1_guides.fa", include_pam=True)
    print("  -> rai1_crispor.fa (for CRISPOR batch)")
    print("  -> rai1_guides.fa (guide sequences)")

    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"Top candidate: {guides.iloc[0]['sequence']}")
    print(f"  Position: {guides.iloc[0]['position']:+d} bp from TSS")
    print(f"  IR Coherence: R̄ = {guides.iloc[0]['coherence_R']:.3f}")
    print(f"  Status: {guides.iloc[0]['go_no_go']}")
    print(f"\nTherapeutic target: {scan['optimal_level']*100:.0f}% RAI1 expression")
    print(f"Required CRISPRa boost: +{(scan['optimal_level'] - 0.5)*100:.0f}% from baseline")
    print("\nNext steps:")
    print("  1. Submit rai1_crispor.fa to crispor.tefor.net")
    print("  2. Validate top guides with genome-wide off-target check")
    print("  3. Partner with wet lab for cell culture validation")


if __name__ == "__main__":
    main()
