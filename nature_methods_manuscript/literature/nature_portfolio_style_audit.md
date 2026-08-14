# Nature Portfolio Writing and Figure Audit

## What the Target Journal Demands

Nature Methods accepts computational and machine-learning methods for
biological data, but asks a method Article to establish performance,
reproducibility, general applicability and potential to enable new biology.
The Article format permits an unreferenced 150-word abstract, a 3,000-word main
text target and up to six figures and/or tables. These constraints favor a
single evidence ladder rather than an internal report with many equivalent
plots.

## Verified Comparative Corpus

| Paper | Narrative movement | Figure logic relevant to this project | Specific lesson |
| --- | --- | --- | --- |
| Avsec et al., *Effective gene expression prediction from sequence by integrating long-range interactions*, Nat. Methods 18, 1196-1203 (2021), doi:10.1038/s41592-021-01252-x | Regulatory gap -> minimal architectural move -> held-out prediction -> enhancer and variant validation. | Fig. 1 combines architecture with global held-out performance; later figures move from mechanism to eQTL and direct mutagenesis. | Do not isolate a network cartoon from its first decisive validation. |
| Linder et al., *Predicting RNA-seq coverage from DNA sequence as a unifying model of gene regulation*, Nat. Genet. 57, 949-961 (2025), doi:10.1038/s41588-024-02053-6 | RNA coverage is both difficult and useful -> model -> held-out coverage fidelity -> transcript biology -> variant endpoints. | The opening figure pairs architecture, a real held-out locus and global track/gene statistics; later figures separate distinct regulatory endpoints. | Use both aggregate metrics and observed-versus-predicted coverage at selected held-out loci. |
| Avsec et al., *Advancing regulatory variant effect prediction with AlphaGenome*, Nature 649, 1206-1218 (2026), doi:10.1038/s41586-025-10014-0 | Unified model and training regime -> broad benchmark -> track fidelity -> variant mechanisms -> ablations. | Fig. 1 gives the architecture, training and broad external-performance overview; track, QTL, locus and ablation evidence then have distinct roles. | Keep external benchmark, biological utility and component attribution in separate figures. |
| Dalla-Torre et al., *Nucleotide Transformer: building and evaluating robust foundation models for human genomics*, Nat. Methods 22, 287-297 (2025), doi:10.1038/s41592-024-02523-z | Benchmark first -> representations -> variant use -> efficiency and scaling. | Training/task design appears before downstream claims; scale and resource requirements are visible rather than buried in Methods. | Make the matched from-scratch control, training budget and compute disclosure explicit. |
| Fudenberg et al., *Predicting 3D genome folding from DNA sequence with Akita*, Nat. Methods 17, 1111-1117 (2020), doi:10.1038/s41592-020-0958-x | Model -> biological correspondence -> sequence perturbation -> interpretation -> cross-system evidence. | The paper moves beyond correlations to perturbation-linked interpretation. | Attribution is useful only when it reaches an independently testable regulatory endpoint. |
| Zhou and Troyanskaya, *Predicting effects of noncoding variants with deep learning-based sequence model*, Nat. Methods 12, 931-934 (2015), doi:10.1038/nmeth.3547 | Task -> single-nucleotide prediction -> eQTL/GWAS utility. | Three compact figures carry method concept, predictive accuracy and variant utility. | Give each display item one conclusion; avoid redundant performance panels. |
| Kempynck et al., *CREsted: modeling genomic and synthetic cell-type-specific enhancers across tissues and species*, Nat. Methods 23, 946-959 (2026), doi:10.1038/s41592-026-03057-2 | Workflow -> biological systems -> interpretation -> public baseline -> in vivo validation. | Cross-context generalization and experimental validation are separate evidence tiers. | An internal illustration cannot replace independent validation. |
| Spiro et al., *A scalable approach to investigating sequence-to-function predictions from personal genomes*, Nat. Methods 23, 1308-1312 (2026), doi:10.1038/s41592-026-03124-8 | Define the correct personal-genome split -> external cohort validation -> success and non-generalizing boundary. | The split itself is part of the first method figure. | State exactly what is unseen before reporting performance. |
| Huang et al., *Personal transcriptome variation is poorly explained by current genomic deep learning models*, Nat. Genet. 55, 2056-2059 (2023), doi:10.1038/s41588-023-01574-w | Separate across-gene prediction from across-person prediction and allele effects. | Distinct endpoints are not collapsed into one correlation. | For variant use, report direction, ranking and expression association separately. |
| Park et al., *Automated neuron tracking inside moving and deforming C. elegans using deep learning and targeted augmentation*, Nat. Methods 21, 142-149 (2024), doi:10.1038/s41592-023-02096-3 | Worm-specific technical barrier -> method -> benchmark -> biological output -> usable interface. | *C. elegans* is introduced as a rigorous experimental system with practical payoff. | Frame the organism as a tractable testbed for a reusable biological method. |

## Story and Tone Rules Derived from the Corpus

1. Start with a specific predictive limitation, not a claim that the field needs
   another model.
2. State the technical move once, then spend the Results on evidence rather than
   repeated architecture description.
3. Put the independent split and fair-comparison contract in the first two
   figures, before ablation or interpretation.
4. Use one direct held-out locus illustration alongside quantitative summaries;
   it gives the reader an observable object behind a correlation.
5. Treat components, external generalization and biological use as different
   questions with different evidence and different display items.
6. Let external biology close the Results. It should test a prespecified use of
   the model, not merely decorate it with a familiar gene or pathway.
7. Write confident, concrete sentences whose scope equals the frozen evidence.
   Avoid generic promotional phrases, defensive scene-setting and redundant
   transitions.

## Visual Grammar for This Manuscript

- Use a schematic-led Fig. 1, then quantitative grids with one dominant hero
  panel per figure.
- Directly label full Model B and major comparators. Keep background controls
  grey and use the focal model color consistently.
- Show paired fold effects where the comparison is matched. Show per-track
  distributions as descriptive heterogeneity, not as inflated replication.
- Use genome-browser-style coverage panels only for preselected held-out loci;
  annotate the context, assay and prediction/observation relationship clearly.
- Reserve strong color for the proposed model and the biologically validated
  signal. Keep audit and control panels visually quieter.

## Sources

- Nature Methods scope and Article requirements: https://www.nature.com/nmeth/aims and https://www.nature.com/nmeth/content
- Nature Methods formatting guidance: https://www.nature.com/nmeth/submission-guidelines/aip-and-formatting
- Individual official paper URLs are given through their DOI links above.
