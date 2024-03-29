# PAW_phospho

Software for processing TMT-labeled phospho peptide enriched proteomics data processed with Proteome Discoverer 3.1 (specifically, the PSM export files).

## Background

Reversible phosphorylation is one of the most important regulatory biological processes. Mass spectrometry is an important technique for characterizing post-translational modifications like phosphorylation. PTMs can be identified and (sometimes) localized to specific residues in peptides from enzymatically digested proteins. Considerable effort over the past couple of decades has gone into techniques to enrich phospho peptides and configure the mass spectrometer to sequence phosphorylated peptides.

Enrichment is necessary because most PTMs, including phosphorylation, are substoichiometric (small numbers of modified proteins out of total protein). Phospho peptides are generally harder to detect than their unmodified counterparts. The phospho group (particularly for S and T residues) is labile and makes peptide fragmentation trickier. Collision energy can be dissipated by loss of the phospho group leaving insufficient energy for peptide bond cleavage.

Enrichment methods have matured and there are many kits and choices available. This is not to say that efficient enrichment is trivial. Much has been learned over the years about better modes of fragmentation for phosphorylated peptides and most mass specs have standard options now. Studying phosphorylation is a key area of proteomics, and techniques for sample processing and mass spectrometry have advanced tremendously.

## The Problem

Unfortunately, data analysis methods for phospho peptide datasets have not kept pace. One can argue that they never really got there in the first place. Find any published proteomics paper studying phosphorylation and you will see the phospho peptides are identified with a tool/pipeline designed to report proteins. If you are interested in peptides, does it matter that the tool/pipeline is protein-centric? Yes, Virginia, there is a Santa Claus and you do not want to use protein-centric tools for peptide-centric studies.

Protein inference and parsimony logic are critical steps in working backwards from the peptide digest (as determined by the sequenced peptides from the mass spec) to the (likely) proteins in the sample that could have give rise to those peptides. None of that applies to sets of enriched phospho peptides. Phosphorylation is (somewhat) specific (only at a few residues in a protein) and enrichment methods are very powerful (capable of pulling out very low abundance phospho peptides). Many phospho peptides will map to phospho proteins without any other peptide evidence. There is litle information to drive protein inference algorithms.

The typical tools/pipelines first figure out the final protein list. This may involve a protein ranking function/heuristic and some attempt at protein false discovery rate control using target/decoy protein counting. Once the proteins are decided, tables of peptides and PSMs associated with the proteins in the final list are generated. Some tools/pipelines may also make conditional PTM reports. The important point is that any PSM, peptide, or PTM reports are dependent on the final protein list. They will not include all PSMs that passed any PSM-level filtering. The protein inference, parsimony logic, protein identification criteria, and any protein FDR control will reject some PSMs. Some of the phospho peptides will be lost.

Phospho peptide experiments need peptide-centric results reporting. This might include some peptide grouping of quantitative measurements (less noise and fewer multiple tests to correct for), modification reports, site reports, site localizations, etc. No tools/pipelines do this - even after more than two decades of phosphorylation studies.   

## Poor data analysis practices are common

The fact that protein inference is irrelevant (and detrimental) in phospho peptide studies and all phospho peptide studies are processed with tools/pipelines that are 100% dependent on protein inference should have been a **giant red flag**. How can it be that so many very smart folks for decades have not realized that their data analyses of phospho peptide data were fundamentally flawed? The short answer is that tool/pipeline failures are not catastrophic and common bad practices in tool/pipeline design give the appearance that they work for peptide-centric data.

All tools/pipelines are primarily designed to report parsimonious lists of inferred proteins. Protein identification criteria (minimum number of distinct peptides per protein being the most important here) and any protein FDR thresholds affect the final list of inferred proteins. The peptide and PTM reports commonly used in downstream statistical analyses of phospho peptides are conditioned on the final protein lists. Only the peptides and PTMs associated the final list of proteins are reported. The default protein ID criteria in most tools is what makes them (sort of) work for phospho peptide data.

### Single-peptide-per-protein identifications

What protein inference/criteria do tools/pipelines adopt that lets most phospho peptides end up in the peptide and PTM tables? It is the adoption of single-peptide-per-protein identifications and the abandonment of the "two peptide rule". The two peptide rule was the de facto standard for many years until [this paper by Gupta and Pevzner](https://pubs.acs.org/doi/abs/10.1021/pr9004794) argued against the two-peptide rule. Yes, some "one-hit-wonders" might be correct, but the difficulty in distinguishing them from incorrect proteins is why they were discarded in the place. Quantification is what is import these day, not identification. Single peptide per protein IDs will be the lowest abundance proteins and will seldom be quantifiable. Think about it, 100 proteins with a single PSM each is a total of 100 crappy PSMs in datasets that are typically 10s to 100s of thousands of PSMs.

> Even earlier, [ICAT (isotope coded affinity tags)](https://pubmed.ncbi.nlm.nih.gov/10504701/), a peptide-centric cysteine labeling and enrichment strategy produced data where single peptides were observed from proteins. The two-peptide rule was not good for those experiments. That is why the [ProteinProphet paper](http://tools.proteomecenter.org/publications/Nesvizhskii.AnalChem.03.pdf) was designed to report proteins with single peptides.

Rigorous mathematical arguments and practical heuristics are often a poor mix. Students with the maths skills frequently lack the experience to understand/appreciate practical heuristics. An insufficient understanding of the starting data for algorithms and equally poor understanding of how algorithm outputs are used by biologists also contribute to the poor state of proteomics software.

Peptide sequencing errors do not have a simple relationship to protein identification errors. It is not a simple scoring function problem like PSM FDR and not a simple target/decoy counting exercise. It is a probability problem:

![Protein Errors](images/protein_errors.png)

There is a one-to-one relationship between PSM errors and protein errors when single peptide per protein results are allowed. A dataset with one hundred thousand identified MS2 scans at a 1% PSM error rate produces 1,000 incorrect PSMs that result in roughly 1,000 incorrect proteins. If we assume the data is from human and we have about 20,000 proteins in the canonical FASTA file, then 1,000 incorrect PSMs would give rise to about 25 incorrect proteins if two peptides per protein was required. For fun, if we required three peptides per protein, the estimated number of incorrect proteins would be less than one (0.4). The number of peptides per protein is a very powerful protein noise filter (reducing incorrect proteins from 1,000 to 25 to 0.4, for this example).

This poor understanding of the relationship between PSM errors and protein errors, has led to protein error control processes analogous to PSM FDR control. PSM error control relies on some scoring function to distinguish incorrect matches from correct matches. Note that PSM score distributions involve all PSMs in an experiment so that the target incorrect distribution can be compared to the decoy distribution and the suitability of the decoy FASTA sequences verified. By analogy, an *ad hoc* protein "scoring" function must be constructed to distinguish confident protein IDs from less confident protein IDs. Common functions are sums of search engine scores, sums of some PSM probability, etc., where the PSMs have already undergone FDR filtering (this is a biased subset of protein IDs). This is just glorified spectral counting. The target and decoy proteins get ranked by these scores and protein FDR is computed analogously to PSM FDR via list traversals. The protein ranking mostly puts proteins with single peptides on the bottom, then proteins with two peptides, then proteins with three peptides, etc. The protein FDR trimming mostly removes single peptide proteins that had lower scoring PSMs. This process masquerades as a statistical analysis but is little more than a weighted "peptides per protein" rule.

That said, single peptide per protein identifications followed by some protein FDR analysis retains most PSMs and only kicks out some of the lowest scoring matches. The PSM loss is not great and largely goes un-noticed. However, all tools/pipelines probably have a minimum peptides per protein cutoff (usually defaulting to one) that can be set to higher values. Any setting other than one peptide per protein will have serious negative effects on phospho peptides. Tools/pipelines appearing to work for phospho peptides mostly do so completely by accident.

### Site numbers depend on the choice of protein sequences

When was the last time (maybe the fist time?) you read the Methods section of a phospho proteomics paper and there was any mention of how the protein sequences were selected so that the modified peptide site numbers reported by the search engine would facilitate comparison to literature and downstream databases (Uniprot.org, PhosphoSitePlus.org, PathwayCommons.org, etc.). Phosphorylation residues and site numbers have an implicit dependence on a single specific protein sequence. Do you know the full, exact protein sequence that was used to define a residue/position for a peptide/PTM in a publication or database resource?

It seems mind boggling to me that there is no discussion of protein sequence choice in phospho peptides studies that I have ever seen. Amino acid residue and residue position counted from the protein N-terminus is how phospho sites are reported and tabulated. For example, UniProt records have one protein sequence listed and modification sites are given in that protein sequence reference frame. Any alternatively spliced protein forms are described in Swiss-Prot records but no records list sites in any alternative protein reference frames. Alternative protein forms (insertions, deletions, truncations, protein processing, etc.) have different sequences and produce different site numbers (for some of the sites). Canonical sequences may or may not have initial Met residues that can cause off-by-one site number errors. Secreted proteins have signal peptides that are not present in the secreted protein forms but those proteins will probably have the signal peptide at the beginning of the FASTA sequence. There are also many mature proteins processed from longer protein sequences. What is the site numbering reference frame for those mature proteins?

An important part of solving problems is thinking them all the way through. If we do not ever talk about the protein reference frames for phospho site numbering, just how well has this problem been thought through? This is another red flag to me.

### Site localization is often impossible

Site localization algorithms are often (mis)used in phospho peptide studies. This is one of those "looks great on paper" but "comes up short in practice" situations. Like too many spectral processing algorithms, site localization can work for abundant peptides with high signal-to-noise and good peptide fragmentation. A lot of MS2 spectra in real world studies are not high quality.

Situations were site localization might work:

- peptides with one PTMs
- PTMs in the middle of the peptides
- possible sites in peptide are not at adjacent residues
- peptides do not have co-eluting positional isomers

What is common in phospho peptide data?

- multiple phospho groups per peptide
- phospho sites close to the peptide N- or C-terminus
- many S, T, or Y residues per peptide
- many adjacent S, T, or Y residues
- (likely) co-eluting positional isomers

PTM work in proteomics is really hard. There are many assumptions in site localization algorithms and those assumptions are only met for a subset of the peptides. Treating site localization probabilities as rigorous computed values and invoking hard site localization probability cutoffs is not a smart thing to do. How well site localization works in "the wild" has not been proven well enough for me to trust using site localization as an upstream data filter. Another red flag commonly seen in publications is using site localization probabilities to filter phospho peptide data.

## Quantitative testing of peptide-centric data

The above issues are related to identification of phospho peptides and reporting phospho peptide sites. Biological experiments involve comparing sample groups to see what phospho peptides are differentially abundant. Phospho peptide enrichment proteomics will run headlong into bottom-up quantitative proteomics. We need to think about how to summarize phospho peptide quantitative measurements for safer statistical testing.

Most quantitative proteomics data is noisy at the measurement level (individual scans) and benefits from averaging measurements (provided you don't average away biological effects). Peptides can have multiple charge states, more than one MS2 scan can be acquired from the same peptide, Met residues can be oxidized during sample handling, peptides can be present in more than one fraction (in fractionated samples), peptides can have phospho groups attached to different residues, and peptides can have different numbers of phospho groups.

Combining multiple PSMs makes sense for most of these cases. I view site localization as not very robust, so combining PSMs when the phospho site is on different residues in the same peptide sequence also makes sense (to me). The total number of phospho groups in a protein region might have biological significance, so combining peptides that have different numbers of phospho groups might not be as okay. Combining as many PSMs as possible before statistical testing will reduce measurement noise and reduce the number of statistical tests (reducing effects of multiple testing corrections).    

## Summary of main concepts

*Slide 1*

![Slide 1](images/Slide1.png)

There are many hidden assumptions in protein-centric proteomics data analyses. Their goal is not to retain all the peptides or provide the best reporting of peptides. Their focus is on the proteins/genes.

---

*Slide 2*

![Slide 2](images/Slide2.png)

The enrichment process occurs after protein digestion and pulls out phosphorylated **peptides**. The mass spectrometer sequences **peptides**. The search engine assigns **peptide** sequences to fragment ion specta. The taget/decoy method determines correctness of **peptide** assignments. Proper summarization and reporting of peptides is what should be done for phospho peptide studies.

---

*Slide 3*

![Slide 3](images/Slide3.png)

Almost every tool/pipeline commonly used in proteomics data analyses is designed to report **proteins**. A large fraction of the data processing in these tools/pipelines is devoted to the protein part. The important question is whether or not this protein part distorts the summarized and reported peptides.

---

*Slide 4*

![Slide 4](images/Slide4.png)

Whether or not these common tools/pipelines (Proteome Discoverer, MaxQuant, the Trans-Proteomic Pipeline, FragPipe, etc.) can report acceptable results for phospho peptide studies depends on how they are configured. The settings need to pass through as many phospho peptides as possible into the final results files. Some loss of usable phospho peptides will always occur, but it can be kept to low levels with careful configuration of the tools/pipelines.

---

*Slide 5*

![Slide 5](images/Slide5.png)

That leaves us yet again with no proper tool for a very common type of proteomics experiment. I did an original version of this Python script for Proteome Discoverer 1.4 back in 2015 when I could not find any appropriate tools for processing TMT-labeled phospho peptide data. Nine years later I am updating the PD 1.4 script for PD 3.1 because I still do not see any appropriate tools for the job. It would sure be nice if the decades of bioinformatic/computational work in proteomics solved some of these basic problems...

---

*Slide 6*

![Slide 6](images/Slide6.png)

Multiple plex TMT experiments can be done for **protein-level** quantitation easily using the internal reference scaling experimental design. There are many aspects of protein-level summarization that help this method work. Most of those aspects are lost at the peptide level. An IRS-like approach to peptide-centric TMT data will not work as well as it does for protein-centric studies. It is best to restrict phospho peptide studies to single plex experiments. Luckily, we now have 18 channels to use.  

---

## What `PD3.1_TMT_phospho_processer` script does

Proteome Discoverer (PD) v3.1 can be configured to provide phospho peptide PSM lists with TMT intensities and site localization information from ptmRS. Those PSM lists can be filtered by Percolator q-value to provide an appropriately small number of incorrect PSMs. PSMs can be combined to improve quantitation with some loss of site localization information. The large PD export tables can be simplified for easier biological interpretation.

PD SEQUEST setting should be adjusted to provide better PSM ID rates (see below). Variable PTMs are set to the minimum (phosphorylation at S, T, or Y and oxidation of Met). Some consensus workflow settings also need adjustment: protein ID criteria is one peptide per protein and the peptide and protein FDR settings in the consensus workflow are set to 0.99 to make sure all identified PSMs make it into the PSM results table.

> PD default parameters should not be assumed to be optimal in real world samples. For example, in SPS-MS3 TMT data where the ion trap is used for peptide IDs, the default SEQUEST settings mass tolerance setting of 10 PPM precursor and 0.6 Da fragment ion are not as good as 20 PPM and 1.0005 Da. A minimum peptide length of 7 instead of 6 should always be used. The default reporter ion quantitative measure defaults to "automatic" (signal-to-noise ratio for Orbitraps) and should be set to "intensities" (reporter ion peak heights).

The script reads in a tab-delimited text version of the PSM results table and the protein FASTA file used in the SEQUEST search (for protein descriptions and to map peptide site numbers into protein site numbers). Other results tables from PD are irrelevant.

Peptide IDs are filtered to a user set q-value cutoff and delta mass window (variables set at the top of the Python script). If using a 20 PPM precursor tolerance in the SEQUEST search, a +/- 10 PPM cutoff would be a good choice (depending on the instrument's mass calibration). The q-value cutoff choice depends on the dataset size and how many PSMs are identified. For example, 50,000 PSMs at a q-value cutoff of 0.01 would result in 500 incorrect phospho peptides. More strict q-value cutoffs may need to be set. The script has some console output that will list the final number of phosphorylated PSMs and an estimate of the number of incorrect matches. Q-value cutoff of 0.005 or 0.001 may be needed to get the final number of incorrect PSM down to an acceptable level. [I have no idea what an acceptable maximum level of incorrect phospho sites should be, though.]

In addition to the delta mass and q-value filters, non-phosphorylated peptides are (optionally) excluded (they can be kept if desired). Enrichment methods may not be as specific for phosphorylated peptides as advertised and the number of non-phosphorylated peptides after enrichment may be large. The filtered PSMs are written to a new PSM results table. PSMs and their reporter ion signals are combined. The combining algorithm uses the unmodified peptide sequence string and number of phospho groups as a combination key for PSM grouping. Multiple copies of peptides, peptides in multiple charge states, peptides with or without oxidized Met, and peptides with phospho groups at different locations within the peptide will be combined. Information (lists) about what peptides (PSMs) were combined is also retained. A new PSM results table for combined peptides is produced. A simplified (only essential columns) combined PSM table is written for use in statistical testing notebooks.

The script lets the user browse to the PSM export file and browse to the FASTA sequence file. Any PSM filters are hard coded at the top of the script. The PSM results file parsing is not adaptive. It expects all of the default columns typically present in the PSM results table from PD 3.1. The tab-delimited text output files written by the script use the base file name of the PD PSM export table and append some suffixes.

If you are interested in using this script, you might want to contact me by email for assistance. This tool is in its early stages. It may need tweaks to run without errors on other datasets. The PD results files depend on how PD is configured. They may have fewer or more columns that the data I tested with.

---

Phil Wilmarth <br> PSR Core, OHSU <br> February 29, 2024
