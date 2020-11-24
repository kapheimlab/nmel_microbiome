# Reanalysis of _Nomia melanderi_ microbiome data with Qiime2
## Karen M. Kapheim
## July 3, 2020

NOTE: This is a clean version of the nmel microbiome qiime2 pipeline.
I ended up deleting all output files and rerunning the whole pipeline
after feeling like I had made a mess of things trying some alternative methods


Working through tutorials on https://docs.qiime2.org/2018.8/interfaces/q2cli/.

Have updated to the newest QIIME version since starting the reanalysis last November, so re-doing now.

Working in interactive mode.

## Set-up

Working through tutorials on https://docs.qiime2.org/2018.8/interfaces/q2cli/.

```
module load anaconda3/2019.03 qiime2/2019.4
cd ../kapheim-group1/nmel_microbiome/
```

## Import sequence data

Working from https://docs.qiime2.org/2018.8/tutorials/importing/.

We have data in the `PairedEndFastqManifestPhred33` format, according to report from UIUC.

> In this variant of the fastq manifest format, there must be forward and reverse \
read fastq.gz fastq files for each sample id. As a result, each sample id is \
represented twice in this file: once for its forward reads, and once for its \
reverse reads. This format assumes that the PHRED offset used for the positional \
quality scores in all of the fastq.gz fastq files is 33.

Make a sample manifest for all V4 seqs in `seqs_V4`.

headers:
`sample-id,absolute-filepath,direction`

```
cd qiime2
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path ../nmel_microbiome_manifest.csv \
  --output-path ./import/v4-pe-demux.qza \
  --input-format PairedEndFastqManifestPhred33
```
This step took a long time to complete.


## Trim adapters

Following
https://forum.qiime2.org/t/demultiplexing-and-trimming-adapters-from-reads-with-q2-cutadapt/2313  https://docs.qiime2.org/2018.8/plugins/available/cutadapt/trim-paired/

When running through qiime1, I had determined that additional adapter trimming was \
not necessary. However, these reads had been joined first. So it seems maybe this \
is necessary now? Will see what the effect is...

Just trimming from the 5' end, not the 3' end. Spot-checking, I couldn't find many
adapters/primers in the 3' ends. I used both the primer sequence and the CS primers \
that UIUC recommended trimming, because the function `--p-front0f` trims the \
sequence and everything preceding it. So this would cover catching seqs where \
one of the other sequence was missing for some reason.

```
cd ./import
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences v4-pe-demux.qza \
  --p-cores 3 \
  --p-front-f GTGYCAGCMGCCGCGGTAA \
  --p-front-f AGACCAAGTCTCTGC \
  --p-front-r GGACTACNVGGGTWTCTAAT \
  --p-front-r TGTAGAACCATGTC \
  --o-trimmed-sequences ../trimmed/v4-pe-trimmed.qza \
  --verbose
cd ../trimmed
qiime demux summarize --i-data v4-pe-trimmed.qza --o-visualization v4-pe-trimmed.qzv
qiime tools view v4-pe-trimmed.qzv
```

Saved the read counts as csv file
`per-sample-fastq-counts(1).csv` in home directory.

## Denoising and dereplicating

Using dada2

Following:
https://docs.qiime2.org/2019.7/tutorials/overview/
https://docs.qiime2.org/2019.7/tutorials/moving-pictures/
https://docs.qiime2.org/2019.7/tutorials/atacama-soils/#atacama-demux


> Dada2 is a pipeline for detecting and correcting (where possible) Illumina \
amplicon sequence data. As implemented in the q2-dada2 plugin, this quality \
control process will additionally filter any phiX reads (commonly present in \
marker gene Illumina sequence data) that are identified in the sequencing data, \
and will filter chimeric sequences.

#### Denoising
> To put it simply, these methods filter out noisy sequences, correct errors in \
marginal sequences (in the case of DADA2), remove chimeric sequences, remove \
singletons, join denoised paired-end reads (in the case of DADA2), and then \
dereplicate those sequences.

Dada2 will also join paired-ends.

Used `v4-pe-trimmed.qzv` visualization to choose length to truncate sequences to based on where median quality score
drops below 30.

```
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs v4-pe-trimmed.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 213 \
  --p-trunc-len-r 191 \
  --o-table /uufs/chpc.utah.edu/common/home/kapheim-group1/nmel_microbiome/qiime2/dada2_out/nmel_v4_q2dada2_table.qza \
  --o-representative-sequences /uufs/chpc.utah.edu/common/home/kapheim-group1/nmel_microbiome/qiime2/dada2_out/nmel_v4_q2dada2_repseqs.qza \
  --o-denoising-stats /uufs/chpc.utah.edu/common/home/kapheim-group1/nmel_microbiome/qiime2/dada2_out/nmel_v4_q2dada2_denoising-stats.qza
```

*Summarize the denoising*

First get metadata
Copied sample data from qiime1 statistical analysis
`uufs/chpc.utah.edu/common/home/kapheim-group1/nmel_microbiome/statistical_analysis/sample_data.txt`
But needed to modify a few things
1. rename sample to id and sampleID to sample
2. replace '.' in id with '_' to match sequence files
3. add in missing samples
  a. blanks and samples I had previously eliminated due to low read OTU number
    -- (samples M36, M45, M40, P1.03)
  b. did this manually from 'microbiome_samples_for_sending_11apr2017_kmk.xls'

```
cd ../dada2_out/
qiime feature-table summarize \
  --i-table nmel_v4_q2dada2_table.qza \
  --o-visualization nmel_v4_q2dada2_visualization_table.qzv \
  --m-sample-metadata-file /uufs/chpc.utah.edu/common/home/kapheim-group1/nmel_microbiome/sample_data.txt
qiime feature-table tabulate-seqs \
  --i-data nmel_v4_q2dada2_repseqs.qza \
  --o-visualization nmel_v4_q2dada2_repseqs.qzv
```

Visualize the denoising stats

```
qiime metadata tabulate \
  --m-input-file nmel_v4_q2dada2_denoising-stats.qza \
  --o-visualization nmel_v4_q2dada2_denoising-stats.qzv
```

To view these artifacts:

```
qiime tools view nmel_v4_q2dada2_denoising-stats.qzv
qiime tools view nmel_v4_q2dada2_repseqs.qzv
qiime tools view mel_v4_q2dada2_visualization_table.qzv
```

## Taxonomic Classification

#### Import classifiers for training datasets

Following https://docs.qiime2.org/2019.7/tutorials/feature-classifier/

Linked to data from https://docs.qiime2.org/2019.7/data-resources/

```
mkdir training-feature-classifiers
cd training-feature-classifiers
wget https://www.arb-silva.de/fileadmin/silva_databases/qiime/Silva_132_release.zip
unzip Silva_132_release.zip
rm -R __MACOSX/
wget ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_8_otus.tar.gz
```

#### Training on Silva 132 database

*Import sequences and taxonomy*

Sequences

```
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path /uufs/chpc.utah.edu/common/home/kapheim-group1/nmel_microbiome/qiime2/training-feature-classifiers/SILVA_132_QIIME_release/rep_set/rep_set_16S_only/99/silva_132_99_16S.fna \
  --output-path /uufs/chpc.utah.edu/common/home/kapheim-group1/nmel_microbiome/qiime2/training-feature-classifiers/SILVA_132_QIIME_release/silva_132_99_16S.qza
```

Now taxonomy

Using 7 level taxonomy

```
qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path /uufs/chpc.utah.edu/common/home/kapheim-group1/nmel_microbiome/qiime2/training-feature-classifiers/SILVA_132_QIIME_release/taxonomy/16S_only/99/taxonomy_7_levels.txt\
  --output-path /uufs/chpc.utah.edu/common/home/kapheim-group1/nmel_microbiome/qiime2/training-feature-classifiers/SILVA_132_QIIME_release/silva_132_99_16S_7level_tax.qza
```

*Extract reference reads*

>It has been shown that taxonomic classification accuracy of 16S rRNA gene sequences\
 improves when a Naive Bayes classifier is trained on only the region of the target \
 sequences that was sequenced (Werner et al., 2012). This may not necessarily \
 generalize to other marker genes (see note on fungal ITS classification below). \
 We know from the Moving Pictures tutorial that the sequence reads that we’re trying \
 to classify are 120-base single-end reads that were amplified with the 515F/806R \
 primer pair for 16S rRNA gene sequences. We optimize for that here by extracting \
 reads from the reference database based on matches to this primer pair, and then \
 slicing the result to 120 bases.

```
qiime feature-classifier extract-reads \
  --i-sequences /uufs/chpc.utah.edu/common/home/kapheim-group1/nmel_microbiome/qiime2/training-feature-classifiers/SILVA_132_QIIME_release/silva_132_99_16S.qza \
  --p-f-primer GTGYCAGCMGCCGCGGTAA \
  --p-r-primer GGACTACNVGGGTWTCTAAT \
  --p-min-length 100 \
  --p-max-length 400 \
  --o-reads /uufs/chpc.utah.edu/common/home/kapheim-group1/nmel_microbiome/qiime2/training-feature-classifiers/SILVA_132_QIIME_release/silva_132_99_16S_refseqs.qza
```
This step took a long time to complete.

*Train the classifier*



7 level taxonomy

```
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads /uufs/chpc.utah.edu/common/home/kapheim-group1/nmel_microbiome/qiime2/training-feature-classifiers/SILVA_132_QIIME_release/silva_132_99_16S_refseqs.qza \
  --i-reference-taxonomy /uufs/chpc.utah.edu/common/home/kapheim-group1/nmel_microbiome/qiime2/training-feature-classifiers/SILVA_132_QIIME_release/silva_132_99_16S_7level_tax.qza \
  --o-classifier /uufs/chpc.utah.edu/common/home/kapheim-group1/nmel_microbiome/qiime2/training-feature-classifiers/SILVA_132_QIIME_release/silva_132_99_16S_7tax_classifier.qza
```
This took a couple of hours

#### Classify rep sequences

*Test the classifier*

> Verify that the classifier works by classifying the representative sequences \
> and visualizing the resulting taxonomic assignments.

7 level taxonomy

```
cd /uufs/chpc.utah.edu/common/home/kapheim-group1/nmel_microbiome/qiime2/classified
qiime feature-classifier classify-sklearn \
  --i-classifier /uufs/chpc.utah.edu/common/home/kapheim-group1/nmel_microbiome/qiime2/training-feature-classifiers/SILVA_132_QIIME_release/silva_132_99_16S_7tax_classifier.qza \
  --i-reads /uufs/chpc.utah.edu/common/home/kapheim-group1/nmel_microbiome/qiime2/dada2_out/nmel_v4_q2dada2_repseqs.qza \
  --o-classification nmel_v4_SILVA_132_99_16S_taxonomy.qza
qiime metadata tabulate \
  --m-input-file nmel_v4_SILVA_132_99_16S_taxonomy.qza \
  --m-input-file /uufs/chpc.utah.edu/common/home/kapheim-group1/nmel_microbiome/qiime2/dada2_out/nmel_v4_q2dada2_repseqs.qza \
  --o-visualization nmel_v4_SILVA_132_99_16S_taxonomy.qzv
```

The classify-sklearn step took a long time.

Received error about white space for visualization step

>There was an issue with viewing the artifact
nmel_v4_SILVA_132_99_16S_taxonomy.qza as QIIME 2 Metadata:
CategoricalMetadataColumn does not support values with leading or trailing
whitespace characters. Column 'Taxon' has the following value:
'D_0__Bacteria;D_1__Rokubacteria;D_2__NC10;D_3__Rokubacteriales;D_4__uncultured bacterium '


## Generate a phylogenetic tree de novo


Use existing pipelines to generate a de novo tree that includes all sequences.

https://docs.qiime2.org/2020.2/tutorials/phylogeny/

[Scroll down to 'Pipelines' near the bottom]

> This pipeline will start by creating a sequence alignment using MAFFT, after \
which any alignment columns that are phylogenetically uninformative or ambiguously \
aligned will be removed (masked). The resulting masked alignment will be used to \
infer a phylogenetic tree and then subsequently rooted at its midpoint. Output files\
 from each step of the pipeline will be saved. This includes both the unmasked and \
masked MAFFT alignment from q2-alignment methods, and both the rooted and unrooted \
phylogenies from q2-phylogeny methods.

```
cd /uufs/chpc.utah.edu/common/home/kapheim-group1/nmel_microbiome/qiime2/phylogeny
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences /uufs/chpc.utah.edu/common/home/kapheim-group1/nmel_microbiome/qiime2/dada2_out/nmel_v4_q2dada2_repseqs.qza \
  --output-dir nmel_v4_mafft-fasttree-output_7tax
```

## Remove mitochondria and chloroplast

Decided to do this before looking at rarefaction curves, in case these taxa are disproportioned

Following: https://usda-ars-gbru.github.io/Microbiome-workshop/tutorials/qiime2/

#### Make a list of taxa identified as mitochondria or chloroplast

```
cd /uufs/chpc.utah.edu/common/home/kapheim-group1/nmel_microbiome/qiime2
# Export taxonomy data to tabular format
qiime tools export \
  --output-path classified \
  --input-path ./classified/nmel_v4_SILVA_132_99_16S_taxonomy.qza
# search for matching lines with grep then select the id column
cd ./classified
grep -v -i "mitochondria" taxonomy.tsv > taxonomy_no-mito.tsv
grep -v -i "chloroplast" taxonomy_no-mito.tsv > taxonomy_no-mito-chloro.tsv
cut -f 1 taxonomy_no-mito-chloro.tsv > taxIDs_no-mito-chloro.tsv
awk 'NR != 1' taxIDs_no-mito-chloro.tsv > taxIDs_no-mito-chloro_nohead.tsv
```

#### Convert feature table to biom file

```
cd /uufs/chpc.utah.edu/common/home/kapheim-group1/nmel_microbiome/qiime2/
# export data to biom
qiime tools export \
  --input-path ./feature_tables/nmel_v4_q2dada2_table.qza \
  --output-path ./feature_tables/
cd ./feature_tables/
mv feature-table.biom nmel_v4_q2dada2_table.biom
# filter biom file
cd ../feature_tables
biom subset-table \
  -i nmel_v4_q2dada2_table.biom \
  -a observation \
  --ids ../classified/taxIDs_no-mito-chloro_nohead.tsv  \
  -o nmel_v4_q2dada2_table_no-mito-chloro.biom
# create new qza with filtered biom file
qiime tools import \
  --input-path nmel_v4_q2dada2_table_no-mito-chloro.biom \
  --output-path nmel_v4_q2dada2_table_no-mito-chloro.qza \
  --type FeatureTable[Frequency]
```

## Rarefaction

For a complete list of filtering methods:

https://docs.qiime2.org/2019.10/tutorials/filtering/

```
# Without removing mito and chloro
cd /uufs/chpc.utah.edu/common/home/kapheim-group1/nmel_microbiome/qiime2/rarefaction
qiime diversity alpha-rarefaction \
    --i-table ../feature_tables/nmel_v4_q2dada2_table.qza \
    --i-phylogeny ../phylogeny/nmel_v4_mafft-fasttree-output_7tax/rooted_tree.qza \
    --p-max-depth 5000 \
    --m-metadata-file /uufs/chpc.utah.edu/common/home/kapheim-group1/nmel_microbiome/sample_data.txt \
    --o-visualization nmel_v4-alpha-rarefaction.qzv
# After removing mito and chloro
qiime diversity alpha-rarefaction \
    --i-table ../feature_tables/nmel_v4_q2dada2_table_no-mito-chloro.qza \
    --i-phylogeny ../phylogeny/nmel_v4_mafft-fasttree-output_7tax/rooted_tree.qza \
    --p-max-depth 5000 \
    --m-metadata-file /uufs/chpc.utah.edu/common/home/kapheim-group1/nmel_microbiome/sample_data.txt \
    --o-visualization nmel_v4-alpha-rarefaction_no-mito-chloro.qzv
```

They look pretty much the sample_data

## Statistical analysis

Export to R and analyze with RStudio on laptop.

file `nmel_v4_jul2020.Rmd`
