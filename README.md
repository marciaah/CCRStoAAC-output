# Output from CCRStoAAC
=======================
This repository contains the description for the **CCRs_output.tsv** file, the output from mapping the Constrained Coding Regions (CCRs) in the Human genome, based on [gnomAD](https://gnomad.broadinstitute.org/) and its version 3.0 (76.156  Human genomes, GRCh38), to the amino acids in canonical proteins of [UniProtKB](https://www.uniprot.org/) (version 10-2020), plus extra information regarding gnomAD allele count and allele frequency, quality of regions, protein features and inter-species conservation. For further details, see the [CCRtoAAC R package repository](https://github.com/marciaah/CCRStoAAC.git) and our **manuscript**. We keep the output in this separate repository, so you can also have access to this information without having to install the CCRtoAAC R package.


### Overview
================
Constrained Coding Regions (CCRs) are focal regions in the Human coding genome depleted of protein changing variations (i.e. missense, stop gain/loss, frameshift indels), see [Havrilla et al., 2019, Nature Genetics](https://doi.org/10.1038%2Fs41588-018-0294-6) ( [GitHub](https://github.com/quinlan-lab/ccr) ) for further details. 

Briefly, in the CCRs model, each of the variant-depleted (constrained) regions in the Human coding genome is weighted based on 
 - its length in base pairs 
 - the fraction of individuals (above 50% of total individuals) having at least a 10x sequencing coverage at each bp of the region. 
Then, a linear regression is calculated comparing the weights and the CpG density of the regions, as an indicator of the region mutability upon spontaneous deamination of methylated cytosines. Regions with a greater weighted distance between protein-changing variants than expected based upon their CpG density (residual from the linear regression), are predicted to be under the greatest constraint. The residuals of the regression are ranked in CCRs percentiles (CCRpct) from 0 to 100, where: 

 - CCRpct = 0 : unconstrained regions (i.e. regions *with* gnomAD variants)
 - >0 CCRpct ≤ 100 : constrained regions (i.e. regions *without* gnomAD variants). 
 
 The longer a constrained region and the bigger its CpG content, in general, the higher its CCRpct will be.


These regions were originally identified using the whole exome sequencing data from large cohorts of healthy control populations aggregated in [gnomAD](https://gnomad.broadinstitute.org/) (The Genome Aggregation Database) version 2.0, with GRCh37/hg19 reference genome, including 125.748 human exomes. 

Here, we extended this by re-calculating the CCRs using [gnomAD3.0](https://gnomad.broadinstitute.org/news/2019-10-gnomad-v3-0/) (GRCh38, 71.702 samples of Human genomes) and mapping these regions to the amino acids in Human protein sequences of UniProtKB. We believe that CCRs have the potential to highlight key functional amino acids in both ordered and intrinsically disordered proteins, lying in protein regions which are constrained.


Citation
--------
If you find this information useful for your work, please mention us: 
**Manuscript in preparation**
*Mapping the Constrained Coding Regions in the human genome to their corresponding proteins*,
Marcia A. Hasenahuer, Alba Sanchis-Juan, Roman A. Laskowski, James A. Baker, James D. Stephenson, Christine A. Orengo,F Lucy Raymond, Janet M. Thornton


### Sources of information
==========================
The information presented in this file was produced or obtained using the following methods and databases:

| Method/Database | Description | Version/date of accession |
| --- | --- | --- |
| [CCRs](https://doi.org/10.1038%2Fs41588-018-0294-6) | The [CCRs model pipeline](https://github.com/quinlan-lab/ccr) for [new Ensembl builds](https://quinlan-lab.github.io/ccr/examples/updates/) | Based on: [gnomAD v3.0](https://gnomad.broadinstitute.org/), [Human Genome GRCh38](ftp://ftp.ensembl.org/pub/release-98/gtf/homo_sapiens/Homo_sapiens.GRCh38.98.gtf.gz), annotated with [VEP v101](https://www.ensembl.org/info/docs/tools/vep/index.html) for [Ensembl v101](http://aug2020.archive.ensembl.org/info/data/index.html) transcripts in the [GENCODE basic v35](http://www.ensembl.org/info/genome/genebuild/transcript_quality_tags.html#basic) set |
| [Ensembl](https://www.ensembl.org/index.html) | [GTF files](http://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/) | v101 |
| [gnomAD](https://gnomad.broadinstitute.org/) | VCF files for v3.0 are not longer available to download, refer to [v3.1.2](https://gnomad.broadinstitute.org/downloads) or request me a copy | v3.0 |
| [UniProtKB](https://www.uniprot.org/) | Protein features obtained via [REST API](https://www.ebi.ac.uk/proteins/api/doc/index.html#/features) | March 2021 |
| [VarSite](https://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/VarSite/GetPage.pl?home=TRUE) | data set provided by Roman A. Laskowski | March 2021 |
| [M-CSA](https://www.ebi.ac.uk/thornton-srv/m-csa/) | Mechanism and Catalytic Site Athlas database, [list of catalytic residues](https://www.ebi.ac.uk/thornton-srv/m-csa/api/residues/?format=json) | June 2021 |
| [BioLip](https://zhanggroup.org/BioLiP/) | Database for biologically relevant ligand-protein binding interactions, [list of artifact ligands](https://zhanggroup.org/BioLiP/ligand_list) | April 2021 |
| [MobiDB](https://mobidb.bio.unipd.it/) | Disorder and mobility annotations | April 2021 |
| [ELM](http://elm.eu.org/searchdb.html) | [all instances](http://elm.eu.org/instances.html?q=*) filtered for Human | February 2021 | 

Further updates will include newer versions of the data.

### The columns in the CCRStoAAC_output.tsv file
================================================
In this file, each line represents an amino acid position (uniprotAcc + sequence position) from a UniProtKB/SwissProt canonical protein, with different types of annotations organized 47 in columns as follows:
 
- **uniprotAcc**: UniProtKB/SwissProt canonical protein identifier. Includes the extension “-#” when there are more isoforms available for this protein

- **aac_pos**: amino acid position

- **aac**: amino acid in single letter nomenclature

- **ensembl_gene_name**: Ensembl gene name

- **aac_weighted_pct**: the CCRs weighted percentile (CCRpct) transferred to the amino acid position. For those cases where a position was in the boundary between two regions, i.e. having a **weighted_pct**=0 (unconstrained) for one/two bases of its codon and a  **weighted_pct**>0 (constrained), for the other one/two bases, such amino acids were assigned as unconstrained with **aac_weighted_pct**=0.

- **aac_weighted_pct_cat**: attempts to generate different categories, “discretizing” the “aac_weighted_pct” (see our **Manuscript**), to facilitate further analysis. Different numbers indicate that the residue position is part of:

  - 6 : CCRpct in [99,100], the top most highly constrained regions
  - 5 : CCRpct in [95,99), a highly constrained region
  - 4 : CCRpct in [90,95), a moderately constrained region
  - 3 : CCRpct in [60,90), a lowly constrained region
  - 2 : CCRpct in [30,60), a lowly constrained region
  - 1 : CCRpct in (0,30), a lowly constrained region
  - 0 : CCRpct = 0 , an unconstrained region, i.e. there is at least one protein changing variant affecting any of the genomic bases of its codon
  - -1 : CCRpct not available, because the amino acid is encoded by a genomic base that is part of a genomic region that corresponds to a segmental tandem duplication or that is highly similar to other region in the genome (>=90% identity)  (we identified 827.952 amino acid positions with this situation)
  - -2 : CCRpct not available, because the amino acid position is encoded by a base with low sequencing coverage in gnomAD3.0 (<50% individuals with at least 10x depth of sequencing coverage) (5.056 amino acid positions with this situation)
  - -3 : CCRpct not available, because of none of the previous reasons (44.177 amino acid positions with this situation). This situation is still under investigation, but most of the cases (so far, manually but not exhaustively checked) coincide with having inframe/frameshift indels in gnomAD3.0 affecting such positions. Therefore, it seems that the CCRs pipeline is failing sometimes to deal with these variants. Consider the possibility of marking these amino acid positions as unconstrained
   

- **region_start_end**: captures the boundaries of the CCRs regions, in amino acid start and end positions
- **region_length**: captures the length of the CCRs regions, in amount of consecutive amino acid positions
- **region_flag**: captures if it was possible to map the CCRs percentile for all the amino acids of the region, with the following code:

  - 1 : all the amino acids in the region were correctly assigned a CCRpct
  - 0 : the region is fragmented and some amino acids have *CCRpct not available*
  - -1 : captures the specific amino acid positions in the region that have *CCRpct not available*

Positions with *CCRpct not available* also have -1 in the columns **region_start_end**, **region_length** and **region_flag**

- **warning_prot**: some proteins carry a **warning_prot**=-1 flag. Please, be cautious with these 45 proteins listed in **warning_proteins.tsv** when using their CCRpct. They have either more than one transcript or more than one gene that can encode for them, hence they might have some amino acids with morethan one CCRpct. As an example, the protein [P62805](https://www.uniprot.org/uniprotkb/P62805/entry) (H4 clustered histone) is the product of 14 different genes and some amino acids you will see 14 rows with information.
- **miss_aac**: comma separated list of gnomAD missense variants for the amino acid site
- **missAC**: sum of gnomAD allele counts for the missense variants changing the amino acid site
- **missAF**: comma separated list of gnomAD allele frequencies for the missense variants affecting the amino acid site
- **CONSERV_score**: the [ScoreCons](https://pubmed.ncbi.nlm.nih.gov/12112692/) inter-species conservation score, as obtained from [VarSite](https://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/VarSite/GetPage.pl?uniprot_acc=NONE&template=home.html)

The last 33 columns assign the presence (>=1) or absence (-1) of a protein functional feature, informing if an amino acid site is part of or has:
- **DOMAIN**: a domain annotated in Pfam and/or CATH, information obtained from VarSite
- **REPEAT**: a repeated sequence motif or domain, as annotated in UniProtKB
- **TRANSMEM**: a transmembrane region, as annotated in UniProtKB
- **COILED**: a coiled-coil region, as annotated in UniProtKB
- **LOWCOMPLEXITY**: a low-complexity region in terms of amino acid composition, as annotated in UniProtKB and/or MobiDB (prediction_low_complexity_merge)
- **CATRES**: a catalytic site, as annotated in M-CSA and/or UniProtKB 
- **METAL**: an interaction with a metal ion, as annotated in UniProtKB and/or VarSite
- **LIGAND**: an interaction with a small molecule, as annotated in VarSite and filtered with BioLip, in order to remove biologically non-relevant molecules 
- **PROTEIN**: an interaction with another protein, as annotated in VarSite
- **CROSSLNK**: a protein cross-linking interaction, as annotated in UniProtKB
- **DISULPHIDE**: a disulfide bond, as annotated in UniProtKB and/or VarSite
- **DNA_RNA**: an interaction with DNA or RNA, as annotated in UniProtKB and/or VarSite
- **MOTIF**: a linear motif, as annotated in UniProtKB and/or ELM
- **LIP**: a linear interacting peptite, as annotated in [MobiDB](https://mobidb.bio.unipd.it/)
- **Dis**: an instrinsically disordered or highly mobile region, as assigned in [different features of MobiDB](https://mobidb.bio.unipd.it/about/mobidb). The values in this column represent:

  - -1: not disordered (i.e ordered)
  - 1: disordered, assigned based on the MobiDB features:
    - prediction_disorder_mobidb_lite
    - curated_disorder_merge
    - homology_disorder_merge
    - prediction_disorder_th_50
  - 2: missing coordinates in experimentally determined structures, based on:
    - derived_missing_residues_context_dependent_th_90
    - derived_missing_residues_th_90
  - 3: highly mobile in NMR structures, based on:
    - derived_mobile_context_dependent_th_90
    - derived_mobile_th_90

- **DtoD**: a region that remains disordered upon binding with a protein partner, as annotated in MobiDB with any of these features:
  
  - curated_binding_mode_disorder_to_disorder_merge
  - derived_binding_mode_disorder_to_disorder_priority

- **DtoO**: a region that structurally transitions from disordered to ordered upon binding with a protein partner, obtained from MobiDB:

  - derived_binding_mode_disorder_to_order_priority

- **context_dep**: a region that has alternative binding configurations, which change with the cellular conditions and different partners, obtaiend from MobiDB:

  - derived_binding_mode_context_dependent_priority

- **LLPS**: a region that drives liquid-liquid phase separation (LLPS), as derived from MobiDB and PhasePro, using the features:

  - curated_phase_separation_merge
  - homology_phase_separation_merge

- **SIGNAL**: a signal peptide, as annotated in UniProtKB
- **PROPEP**: a propeptide (cleaved during protein "maturation"), as annotated in UniProtKB
- **TRANSIT**: a region neccesary for the transport of a protein encoded by a nuclear gene to a particular organelle  (e.g. peroxisome, mitochondrion), as annotated in UniProtKB
- **PHOSPHO**: a site that is post-traslationally modified with phosphorylation, as annotated in UniProtKB
- **LIPID**: a site that adquires covalently attached lipid group(s)
- **CARBOHYD**: a glycosylated site (covalently attached glycan group (mono-, di-, or polysaccharide)) as annotated in UniProtKB
- **PTM**: other post-translational modifications, obtained from UniProtKB 
- **REGION**: a region of interest that cannot be described in other subsections of UniProtKB
- **SITE**: a site of interest that cannot be described in other subsections of UniProtKB
- **CDSjunction**: a CDS (CoDing Sequence) junction, involving the boundaries of exons in a transcript, obtained from Ensembl
- **NO_feature**: having none of the previous 29 features
- **Pathogenic_miss**: at least one Pathogenic/Likely_Pathogenic missense variant in ClinVar
- **Benign_miss**: at least one Benign/Likely_Benign missense variant in ClinVar
- **VUS_conflict_miss**: at least one VUS (variant of uncertain significance) or conflictive interpretations missense variant in ClinVar

### How much of the Human UniProt/SP proteins were effectively mapped with CCRpct?
==================================================================================
We were able to assign CCRs percentiles (column **aac_weighted_pct** >= 0) to  at least one amino acid position in **17,366 UniProt/SP canonical proteins**, corresponding to **17,372 genes in human chromosomes 1-22 and X**. This corresponds to **9.825.893 amino acid positions** (unique combinations of **uniprotAcc**+**aac_pos**)

There were in total  **877,185 amino acid positions**, which have **aac_weighted_pct** < 0, i.e. *CCRpct not available* 

### Why was it not possible to obtain CCRs for all amino acids in all UniProtSP canonical sequences?
====================================================================================================
Some amino acid positions are not annotated with CCRpct (**aac_weighted_pct**), this can be explaine by any of these three reasons:

- Mitochondrial and chromosome Y genes are not considered in the CCRs model, because they lack good coverage in gnomAD3.0, therefore the corresponding proteins are not included.

- If there were mismatches between Ensembl protein and UniProt sequences, i.e. <100% coverage and <100% identity, such proteins were not included in the mapping from genomic coordinates to UniProt protein sequences

- Genomic segmental duplications (tandem repeated regions), highly similar regions across the genome (>=90% identity) and also positions with low sequencing coverage (<50% samples with >=10x depth) are not considered by the CCRs pipeline. This is because such regions are of low quality and variants observed there are not reliable. Such unreliable regions generated absence of CCRs percentile (**aac_weighted_pct** < 0) for regions spanning from single amino acids up to complete protein sequences.


