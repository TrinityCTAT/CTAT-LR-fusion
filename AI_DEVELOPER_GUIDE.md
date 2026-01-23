# CTAT-LR-Fusion: AI Developer Onboarding Guide

## Project Overview

**CTAT-LR-Fusion** is a fusion transcript detection pipeline for long-read RNA-seq data (PacBio Iso-seq and Oxford Nanopore) that's part of the Trinity Cancer Transcriptome Analysis Toolkit (CTAT). It can optionally incorporate Illumina short reads for additional validation.

**Current Version:** v1.2.1 (as of Sept 2025)

**Primary Language:** Perl (main driver), with Python utilities

**Main Driver Script:** [ctat-LR-fusion](ctat-LR-fusion) - Perl pipeline orchestrator

---

## Architecture: Two-Phase Detection Pipeline

### Phase 1: Initial Candidate Detection (Lines 420-570 in ctat-LR-fusion)

**Goal:** Identify preliminary fusion candidates from chimeric alignments

**Steps:**
1. **Align long reads** to reference genome using `ctat-minimap2` (custom minimap2)
2. **Extract chimeric alignments** (reads with SA supplementary alignment tags)
3. **Convert to GFF3 format** for structured parsing
4. **Generate `chims_described` report** - maps each chimeric read to gene annotations
5. **Identify preliminary candidates** - aggregate by fusion pair and apply initial filters
6. **Annotate with FusionAnnotator** - add known fusion annotations
7. **Filter by annotation rules** - remove NEIGHBORS_OVERLAP, etc.
8. **Limit candidates** - cap at `max_phase1_candidates` (default: 10,000)
9. **Extract candidate reads** - prepare for phase 2 (unless `--max_rigor` mode)

### Phase 2: Precise Breakpoint Resolution (Lines 580-700 in ctat-LR-fusion)

**Goal:** Create fusion contigs and realign for precise breakpoint determination

**Steps:**
1. **Build fusion contigs** - create mini-genomes for each fusion candidate
2. **Shrink long introns** (optional, default: shrink to ≤1000bp for faster alignment)
3. **Create minimap2 index** for fusion contigs
4. **Add genome decoy** (optional, if `--max_rigor` or target list provided)
5. **Realign reads** against fusion contigs (±genome)
6. **Extract fusion alignments** - parse alignments to fusion contigs
7. **Filter by sequence similarity** - exclude evidence from seq-similar regions
8. **Run FusionInspector** on short reads (if provided)
9. **Merge evidence** from long and short reads
10. **Annotate and filter** - FusionAnnotator, FusionFilter, abundance filters
11. **Generate outputs** - TSV predictions and optional IGV visualization

---

## Key Components

### Core Executables

- **[ctat-LR-fusion](ctat-LR-fusion)** - Main driver (Perl, 1180 lines)
- **[ctat-minimap2/](ctat-minimap2/)** - Custom minimap2 with `--only_chimeric` flag
- **[FusionInspector/](FusionInspector/)** - Submodule for short-read validation
- **[FusionAnnotator/](FusionAnnotator/)** - Submodule for fusion annotation
- **[FusionFilter/](FusionFilter/)** - Submodule for blast/promiscuity filtering

### Critical Utilities (util/)

**Phase 1 - Candidate Identification:**
- `SAM_to_gxf.pl` - Convert BAM/SAM to GFF3 alignment format
- `genome_gff3_to_chim_summary.pl` - **[KEY]** Map chimeric alignments to genes (682 lines)
- `identify_prelim_fusion_transcript_candidates.pl` - **[KEY]** Aggregate and filter candidates (325 lines)
- `retrieve_reads_for_fusion_transcript_candidates.pl` - Extract reads for phase 2
- `revise_fusion_reads_fasta.pl` - Update read set based on candidates

**Phase 2 - Contig Alignment:**
- `LR-FI_fusion_align_extractor.pl` - Extract fusion evidence from contig alignments
- `incorporate_LR_FFPM.pl` - Add FFPM (fusion fragments per million) calculations
- `merge_mm2fusion_FI.py` - Merge long-read and short-read evidence

**Filtering:**
- `filter_LR_fusions_by_evidence_abundance.py` - Apply read count and FFPM thresholds
- `filter_low_pct_dom_iso.py` - Filter by dominant isoform fraction
- `filter_max_candidate_fusions.pl` - Cap number of candidates

**Visualization:**
- `create_ctat-LR-fusion_inspector_igvjs.py` - Generate IGV-reports JSON
- `LR_sam_fusion_read_extractor.pl` - Extract reads for IGV visualization

---

## Critical Algorithm: Chimeric Alignment to Fusion Candidates

### Step 1: Generate `chims_described` ([util/genome_gff3_to_chim_summary.pl](util/genome_gff3_to_chim_summary.pl))

**Input:** GFF3 alignments from minimap2, reference GTF

**Process:**
1. Build interval tree of gene coordinates for fast overlap detection
2. For each read with ≥2 alignment segments:
   - Order segments by transcript coordinates
   - For each adjacent segment pair (left, right):
     - Map to overlapping genes using interval tree
     - Calculate **delta** = distance from alignment breakpoint to nearest exon boundary
     - Consider both sense and antisense orientations
     - Report mapping(s) with minimum combined delta (deltaA + deltaB)

**Key Variables:**
- `deltaA` - Distance from left breakpoint to nearest exon boundary in geneA
- `deltaB` - Distance from right breakpoint to nearest exon boundary in geneB
- `trans_brkptA` - Breakpoint coordinate on transcript for geneA
- `trans_brkptB` - Breakpoint coordinate on transcript for geneB

**Output Format:**
```
#transcript  num_alignments  align_descr(s)  geneA;deltaA;trans_brkptA;chrA:coordA;geneB;deltaB;trans_brkptB;chrB:coordB;geneA--geneB
```

### Step 2: Aggregate and Compute Statistics ([util/identify_prelim_fusion_transcript_candidates.pl](util/identify_prelim_fusion_transcript_candidates.pl))

**Aggregation (lines 126-172):**
- Group all reads by fusion pair (geneA--geneB)
- Collect arrays of: deltaA values, deltaB values, read names, transcript breakpoint distances

**Statistics Computed (lines 175-186):**
- `median_deltaA`, `median_deltaB` - Median distances
- `min_deltaA`, `min_deltaB` - **Minimum distances (used for filtering)**
- `median_trans_brkpt_delta` - Median distance between breakpoints
- `min_trans_brkpt_delta` - Minimum breakpoint separation
- `num_reads` - Total read count

### Step 3: Filter Preliminary Candidates (lines 234-264)

**Filtering Logic (ALL conditions must be met):**

```perl
FFPM = (num_reads / num_total_reads) × 1,000,000

Pass if:
  FFPM >= min_FFPM  AND
  num_reads >= min_num_LR  AND
  (
    (min_deltaA <= MAX_EXON_DELTA  AND  min_deltaB <= MAX_EXON_DELTA)
      OR
    (One delta <= MAX_EXON_DELTA  AND  other <= 1000bp  AND  num_reads > 1)
  )
```

**Rationale for Dual-Threshold Logic:**
- Allows fusions where one breakpoint is very precise (≤50bp from exon)
- But the other is less precise (≤1000bp) if supported by multiple reads
- Captures real fusions with ambiguous breakpoints while maintaining specificity

---

## Key Parameters and Their Effects

### Long Read Evidence Thresholds
- `--min_num_LR` (default: 1) - Minimum long reads with canonical splicing
- `--min_LR_novel_junction_support` (default: 2) - Min LR for non-canonical junctions
- `--min_per_id` (default: 70) - Minimum percent identity for alignments
- `--min_trans_overlap_length` (default: 100) - Min read overlap per gene

### FFPM and Phase Control
- `--min_FFPM` (default: 0.1) - Minimum fusion fragments per million
- `--frac_FFPM_phase1` (default: 0.6) - Fraction of min_FFPM for phase 1 (soft threshold)
- `--max_phase1_candidates` (default: 10,000) - Max candidates entering phase 2

### Breakpoint Precision
- `--max_exon_delta` (default: 50) - Max distance from exon boundary for initial search
- `--snap_dist` (default: 3) - Snap breakpoints within this distance to splice sites

### Short Read Thresholds (if Illumina data provided)
- `--min_J` (default: 1) - Minimum junction reads
- `--min_sumJS` (default: 1) - Minimum sum of junction + spanning reads
- `--min_novel_junction_support` (default: 1) - Min reads for novel junctions

### Performance and Rigor
- `--max_rigor` - Use entire read set + genome decoy in phase 2 (slower, more sensitive)
- `--max_intron_length` (default: 100,000) - Max intron in initial alignment
- `--shrink_intron_max_length` (default: 1,000) - Shrink introns to this length in phase 2
- `--no_shrink_introns` - Disable intron shrinking (slower but preserves intron structure)

### Filtering Options
- `--no_annot_filter` - Skip annotation-based filtering
- `--no_abundance_filter` - Skip final abundance-based filtering
- `--min_frac_dom_iso_expr` (default: 0.05) - Min expression of dominant isoform

---

## Data Flow and File Naming

### Phase 1 Intermediate Files (in fusion_intermediates_dir/)

```
long_reads.fq.mm2.prelim.bam                    # Initial minimap2 alignment
long_reads.fq.mm2.bam                           # Chimeric alignments only
long_reads.fq.mm2.gff3                          # GFF3 format alignments
long_reads.fq.mm2.chims_described               # Mapped to gene annotations
chimeric_read_candidates.preliminary_candidates_info_from_chims_described
                                                # Pre-filter summary statistics
chimeric_read_candidates.preliminary_candidates_info_from_chims_described.read_support_filtered
                                                # Post-filter candidates (FI_listing)
chimeric_read_candidates.preliminary_candidates_info_from_chims_described.read_support_filtered.wAnnot
                                                # With FusionAnnotator annotations
chimeric_read_candidates.preliminary_candidates_info_from_chims_described.read_support_filtered.wAnnot.annot_filter.pass
                                                # After annotation filtering
chimeric_read_candidates.preliminary_candidates_info_from_chims_described.read_support_filtered.wAnnot.annot_filter.pass.max_phase1_candidates
                                                # Final phase 1 candidate list
chimeric_read_candidates.transcripts.fa         # Reads supporting candidates
```

### Phase 2 Intermediate Files

```
LR-FI_targets.fa                                # Fusion contig sequences
LR-FI_targets.gtf                               # Fusion contig annotations
LR-FI_targets.fa.mm2                            # Minimap2 index
LR-FI_targets.gtf.mm2.splice.bed               # Splice junction BED
LR-FI_targets.seqsimilar_regions.gff3          # Sequence-similar regions
LR-FI.mm2.bam                                   # Realignments to fusion contigs
LR-FI.mm2.gff3                                  # GFF3 format
LR-FI.mm2.fusion_transcripts                    # Extracted fusion evidence
LR-FI.mm2.fusion_transcripts.breakpoint_info.tsv # Breakpoint details
LR-FI.mm2.fusion_transcripts.breakpoint_info.tsv.w_LR_FFPM # With FFPM
mm2_and_FI_fusions_merged.tsv                   # Merged LR+SR evidence
```

### Final Outputs

```
ctat-LR-fusion.fusion_predictions.preliminary.tsv        # Pre-filter predictions
ctat-LR-fusion.fusion_predictions.preliminary.abridged.tsv # Without read names
ctat-LR-fusion.fusion_predictions.tsv                     # Final predictions
ctat-LR-fusion.fusion_predictions.abridged.tsv           # Final without read names
ctat-LR-fusion.fusion_inspector_web.html                 # IGV visualization (if --vis)
```

---

## Output Format: Final Predictions TSV

**Core Columns:**
- `#FusionName` - GeneA--GeneB
- `num_LR` - Number of supporting long reads
- `LR_FFPM` - Fusion fragments per million (long reads)
- `LeftGene`, `RightGene` - Gene symbols
- `LeftBreakpoint`, `RightBreakpoint` - Genomic coordinates
- `SpliceType` - ONLY_REF_SPLICE, INCL_NON_REF_SPLICE
- `LR_accessions` - Comma-separated read names (omitted in .abridged.tsv)

**If Short Reads Included:**
- `JunctionReadCount`, `SpanningFragCount` - Read counts
- `JunctionReads`, `SpanningFrags` - Read names (omitted in abridged)

**Annotation Columns (from FusionAnnotator):**
- `annots` - Known fusion annotations (CHIMERDB, MITELMAN, etc.)

**Coding Effect (if --examine_coding_effect):**
- Coding region predictions and frame information

---

## Common Troubleshooting Scenarios

### No Candidates Found After Phase 1
**Check:**
1. Total read count - verify reads were counted correctly
2. FFPM threshold - may be too high for low-depth samples
3. Alignment percent identity - adjust `--min_per_id` if reads have lower quality
4. `chims_described` file - verify chimeric alignments were detected

### Too Many Phase 1 Candidates (>10,000)
**Solutions:**
1. Increase `--min_FFPM` or `--frac_FFPM_phase1`
2. Increase `--min_num_LR`
3. Decrease `--max_exon_delta` for stricter breakpoint requirements
4. Increase `--max_phase1_candidates` if you have compute resources

### Phase 2 Takes Too Long
**Solutions:**
1. Use `--no_shrink_introns` flag is OFF (default shrinks introns)
2. Avoid `--max_rigor` unless specifically needed
3. Reduce candidates in phase 1 (see above)
4. Check if `--LR_bam` option can be used to skip phase 1 on reruns

### Memory Issues
**Check:**
1. Not loading entire alignment file - processes one target at a time (fixed in v1.0.2)
2. ctat-minimap2 memory leak - update to latest version
3. Consider reducing `--max_phase1_candidates`

---

## Special Modes

### Pre-aligned BAM Input (`--LR_bam`)
- Skip initial minimap2 alignment
- Requires properly formatted BAM with SA tags for chimeric alignments
- Useful for reanalysis or when using externally aligned data

### Max Rigor Mode (`--max_rigor`, or implied by target lists)
- Uses entire read set in phase 2 (not just candidate reads)
- Adds genome as decoy to fusion contigs
- More sensitive but much slower
- Automatically enabled with `--incl_fusion_targets` or `--only_fusion_targets`

### Candidate-Only Mode (`--chim_candidates_only`)
- Stops after phase 1 candidate identification
- Useful for quick exploration or troubleshooting
- Output: `chimeric_read_candidates.FI_listing`

### Target-Specific Analysis
- `--incl_fusion_targets <file>` - Include specific fusions in addition to discovered ones
- `--only_fusion_targets <file>` - Only analyze specified fusions
- Format: Tab-separated fusion pairs (GeneA--GeneB)

---

## Dependencies and Submodules

### Git Submodules
- `FusionInspector/` - Short read validation (branch: devel)
- `FusionAnnotator/` - Fusion annotation database
- `FusionFilter/` - Blast and promiscuity filtering
- `ctat-minimap2/` - Fork of lh3/minimap2 (branch: master)

### External Tools Required
- `samtools` - BAM/SAM manipulation
- `minimap2` (or use bundled ctat-minimap2)
- `igv-reports` (if using `--vis` option)
- STAR (if using short reads via FusionInspector)

### Perl Modules (PerlLib/)
- `Pipeliner.pm` - Pipeline execution framework with checkpointing
- `Process_cmd.pm` - Command execution utilities
- `Fasta_reader.pm`, `Fastq_reader.pm` - Sequence file parsing
- `DelimParser.pm` - Tab-delimited file parsing
- Various SAM/BAM and GFF/GTF parsing utilities

---

## Testing

**Test Directory:** [testing/](testing/)

**Quick Test:**
```bash
cd testing
# See testing README for specific test datasets and commands
```

---

## Docker and Singularity

**Build Scripts:** [Docker/](Docker/)
- `build_docker.sh` - Build Docker image
- `push_docker.sh` - Push to Docker Hub
- `make_simg.sh` - Create Singularity image
- `VERSION.txt` - Version for tagging images

**Image Location:** `trinityctat/ctat_lr_fusion:latest`

---

## Version History Highlights

**v1.2.1 (Sept 2025):**
- Phase 1 candidate limit enforcement
- Faster contig builder with preloaded genome
- Gzipped intermediate files for space savings

**v1.1.0 (Feb 2025):**
- Added `--LR_bam` for pre-aligned input

**v1.0.2 (Jan 2025):**
- Memory optimization (process alignments one target at a time)
- ctat-mm2 memory leak fix

**v1.0.0 (Oct 2024):**
- Stable release for ONT and PacBio
- Sequence-similarity filtering
- Excludes neighbor-overlap gene pairs

---

## Quick Reference: Key Code Sections

### Main Driver ([ctat-LR-fusion](ctat-LR-fusion))
- Lines 1-90: Variable declarations and defaults
- Lines 92-223: Usage documentation
- Lines 225-360: Argument parsing and validation
- Lines 420-570: Phase 1 (candidate detection)
- Lines 580-700: Phase 2 (fusion contig alignment)
- Lines 700-830: Filtering and final output
- Lines 830-880: IGV report generation
- Lines 900-1180: Helper functions

### genome_gff3_to_chim_summary.pl
- Lines 60-160: Parse GTF and build interval trees
- Lines 163-295: Main alignment evaluation loop
- Lines 413-580: Map alignments to exon junctions (key algorithm)

### identify_prelim_fusion_transcript_candidates.pl
- Lines 99-123: Main workflow
- Lines 126-186: Parse and aggregate fusion candidates
- Lines 234-264: Filter candidates by FFPM and delta thresholds
- Lines 267-325: Write output summaries

---

## Design Principles

1. **Two-phase approach:** Broad discovery (phase 1) → precise validation (phase 2)
2. **Checkpointing:** All steps use Pipeliner with checkpoints for resumability
3. **Minimum delta strategy:** Use best-aligned breakpoints for candidate selection
4. **Flexible breakpoint tolerance:** Accommodate real fusions with ambiguous junctions
5. **Resource management:** Limit phase 1 candidates to prevent combinatorial explosion
6. **Evidence integration:** Combine long and short read evidence when available
7. **Sequence-aware filtering:** Exclude evidence from sequence-similar regions

---

## Additional Resources

- **Wiki:** https://github.com/TrinityCTAT/CTAT-LR-fusion/wiki
- **CTAT Genome Lib:** https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/
- **CHANGELOG:** [CHANGELOG.txt](CHANGELOG.txt)
- **Issues:** Use GitHub issues for bug reports and feature requests

---

*This guide was created to facilitate AI-assisted development and troubleshooting. Keep it updated as the codebase evolves.*
