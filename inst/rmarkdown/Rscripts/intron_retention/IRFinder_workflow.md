To run any of these commands, need to activate the bioconda IRFinder environment prior to running script.

1. First script creates reference build required for IRFinder

    ```bash
    #SBATCH -t 24:00:00                    # Runtime in minutes
    #SBATCH -n 4
    #SBATCH -p medium                # Partition (queue) to submit to
    #SBATCH --mem=128G        # 128 GB memory needed (memory PER CORE)
    #SBATCH -o %j.out               # Standard out goes to this file
    #SBATCH -e %j.err               # Standard err goes to this file
    #SBATCH --mail-type=END         # Mail when the job ends 

    IRFinder -m BuildRefProcess -r reference_data/
    ```

      >**NOTE:** The files in the `reference_data` folder are sym links to the bcbio ref files and need to be named specifically `genome.fa` and `transcripts.gtf`:
      >
      >`genome.fa -> /n/app/bcbio/biodata/genomes/Hsapiens/hg19/seq/hg19.fa`
      >
      >`transcripts.gtf -> /n/app/bcbio/biodata/genomes/Hsapiens/hg19/rnaseq/ref-transcripts.gtf`

2. Second script (.sh) runs IRFinder and STAR on input file

      ```bash
      #!/bin/bash

      module load star/2.5.4a

      IRFinder -r /path/to/irfinder/reference_data \
      -t 4 -d results \
      $1
      ```

3. Third script (.sh) runs a batch job for each input file in directory

      ```bash
      #!/bin/bash

      for fq in /path/to/*fastq
      do

      sbatch -p medium -t 0-48:00 -n 4 --job-name irfinder --mem=128G -o %j.out -e %j.err --wrap="sh /path/to/irfinder/irfinder_input_file.sh $fq"
      sleep 1 # wait 1 second between each job submission

      done
      ```

4. Fourth script takes output (IRFinder-IR-dir.txt) and uses the replicates to determine differential expression using the Audic and Claverie test (# replicates < 4). analysisWithLowReplicates.pl script comes with the IRFinder github repo clone, so I cloned the repo at https://github.com/williamritchie/IRFinder/. Notes on the Audic and Claverie test can be found at: https://github.com/williamritchie/IRFinder/wiki/Small-Amounts-of-Replicates-via-Audic-and-Claverie-Test.

      ```bash
      #!/bin/bash

      #SBATCH -t 24:00:00                    # Runtime in minutes
      #SBATCH -n 4
      #SBATCH -p medium                # Partition (queue) to submit to
      #SBATCH --mem=128G        # 8 GB memory needed (memory PER CORE)
      #SBATCH -o %j.out               # Standard out goes to this file
      #SBATCH -e %j.err               # Standard err goes to this file
      #SBATCH --mail-type=END         # Mail when the job ends

      analysisWithLowReplicates.pl \
        -A A_ctrl/Pooled/IRFinder-IR-dir.txt A_ctrl/AJ_1/IRFinder-IR-dir.txt A_ctrl/AJ_2/IRFinder-IR-dir.txt A_ctrl/AJ_3/IRFinder-IR-dir.txt \
        -B B_nrde2/Pooled/IRFinder-IR-dir.txt B_nrde2/AJ_4/IRFinder-IR-dir.txt B_nrde2/AJ_5/IRFinder-IR-dir.txt B_nrde2/AJ_6/IRFinder-IR-dir.txt \
        > KD_ctrl-v-nrde2.tab
      ```

5. Output `KD_ctrl-v-nrde2.tab` file can be read directly into R for filtering and results exploration.

6. Rmarkdown workflow:

```r
## Methods

The intron retention analysis was performed using IRFinder version 1.2.3 based on the gene annotations for the hg19 reference genome.

- All potential introns were analyzed for intron retention ratios in each condition, and introns were defined as any region between two exon features in any transcript. The exon coordinates were derived from reference gene annotation files.

- Intronic regions with poor mappability based on pre-determined mapping of synthetic reads to the genome were excluded from the measurable intron area.

- Intron retention was assessed for directional RNA-Seq data (reverse strand).

- The Audic and Claverie test for small numbers of replicates (Max = 3 replicates) was applied to the quantified intron retention events to determine differential expression between Control and NRDE2-KD conditions as described in [this post](http://mimirna.centenary.org.au/irfinder/example1.html).

For more detailed descriptions, see the [IRFinder wiki](https://github.com/williamritchie/IRFinder/wiki).

```{r setup}
library(tidyverse)
library(annotables)

# Small Amounts of Replicates via Audic and Claverie Test
irfinder_AC_test <- read.table("~/bcbio/PIs/frank_slack/slack_intron_retention/irfinder/IRFinder/KD_ctrl-v-nrde2.tab", header = T)

parsed_rownames <- str_split(irfinder_AC_test$Intron.GeneName.GeneID, "/", simplify = TRUE)

irfinder_AC_test <- data.frame(irfinder_AC_test, parsed_rownames)

irfinder_AC_test <- irfinder_AC_test %>% rename(gene=X1,
           ensembl_id=X2,
           splicing=X3)
grch37_description <- grch37[, c("ensgene", "biotype", "description")]

grch37_description <- grch37_description[which(!(duplicated(grch37_description$ensgene))), ]

irfinder_ACtest_merged <- merge(irfinder_AC_test, grch37_description, by.x="ensembl_id", by.y="ensgene")

irfinder_ACtest_merged <- irfinder_ACtest_merged[, -5]

# irfinder_ACtest_merged <- irfinder_ACtest_merged[order(irfinder_ACtest_merged$p.diff), ]
# 
# irfinder_ACtest_merged$rank <- 1:2172
# 
# irfinder_ACtest_merged$factor <- 2172 / irfinder_ACtest_merged$rank
# 
# irfinder_ACtest_merged$padj <- irfinder_ACtest_merged$p.diff * irfinder_ACtest_merged$factor

irfinder_ACtest_merged$padj <- p.adjust(irfinder_ACtest_merged$p.diff, "BH")

irfinder_ACtest_merged <- irfinder_ACtest_merged[order(irfinder_ACtest_merged$padj), ]

write.csv(irfinder_ACtest_merged, "results/slack_irfinder_ACtest_all_results_padj.csv")
sig_irfinder_ACtest <- irfinder_ACtest_merged[which(irfinder_ACtest_merged$padj < 0.05), ]

sig_irfinder_ACtest <- sig_irfinder_ACtest[order(sig_irfinder_ACtest$p.diff),]

```

## Results

There were ### significantly retained introns. The results output for each intron includes the following information, also described in the [IRFinder wiki](https://github.com/williamritchie/IRFinder/wiki) and an [analysis example](http://mimirna.centenary.org.au/irfinder/example1.html). 

- **ensembl_id:** Ensembl ID
- **Chr**: chromosome
- **Start:** start coordinates of intron
- **End:** end coordinates of intron
- **Direction:** strand
- **ExcludedBases:** number of bases within the intronic region that have been excluded from the calculation of intronic coverage because of overlapping features or mapping issues. The 5bp flanking any exon are also excluded, hence all introns exclude at least 10 bases.
- **p.diff:** p-value output from the Audic and Claverie test regarding intron retention differences between conditions. Indicates to what extent intron-retention has significantly changed between Condition A (control) and Condition B (Nrde2-KD)
- **p.increased:** p-value for whether higher intron retention in Nrde2-KD relative to control
- **p.decreased:** p-value for whether lower intron retention in Nrde2-KD relative to control
- **A.IRratio:** intron retention ratio for control. Calculated by: IntronDepth / (max(number of reads that map the 3' flanking exon surrounding the intron and to another exon within the same gene, number of reads that map the 5' flanking exon surrounding the intron and to another exon within the same gene) + IntronDepth)
- **A.IRok:** Warnings about potential biases to the IR calculation. Low coverage, non-uniform coverage of the intron, etc. "NonUniformCover" indicates that the multiple places IRFinder measures the depth of the intron are not consistent. A visual check of the RNA-Seq trace is advised to rule out alternate-TSS or similar. In some cases this warning may simply be triggered by the uneven cover often seen in short-read RNA-Seq experiments.
- **A.IntronCover:** Ratio of bases with mapped reads for control 
- **A.IntronDepth:** number of reads that map over a given bp for control. IntronDepth is the median depth of the intronic region without the excluded regions. It is used to calculate the IRratio. Excluded regions comprise ExclBases and bases with the top and bottom 30% of intronic depth.
- **A.SplicesMax:** number of reads that map to any exon for control? (no documentation)
- **A.SplicesExact:** number of reads that map across the 3' and 5' flanking exons for control 
- **B.IRratio:** intron retention ratio for Nrde-KD. 
- **B.IRok:** Warnings about potential biases to the IR calculation for Nrde-KD.
- **B.IntronCover:** Ratio of bases with mapped for Nrde-KD. 
- **B.IntronDepth:** number of reads that map over a given bp for Nrde-KD.
- **B.SplicesMax:** number of reads that map to any exon for Nrde-KD? (no documentation).
- **B.SplicesExact:** number of reads that map across the 3' and 5' flanking exons for Nrde-KD.
- **replicates:** whether or not there were replicates (no documentation)
- **A1.IRratio:** IR ratio for first replicate for control condition
- **A2.IRratio:** IR ratio for second replicate for control condition
- **A3.IRratio:** IR ratio for third replicate for control condition
- **B1.IRratio:** IR ratio for first replicate for Nrde2-KD condition
- **B2.IRratio:** IR ratio for second replicate for Nrde2-KD condition
- **B3.IRratio:** IR ratio for third replicate for Nrde2-KD condition
- **gene:** associated gene name
- **splicing:** whether there are known transcripts that include this intron. `known-exon` included if there are known transcripts.
- **biotype:** type of transcript (protein_coding, mtRNA, etc.)
- **description:** short description of gene
- **padj:** p-value corrected for multiple testing using the Benjamini-Hochberg/FDR method

Similar to the criteria used in the [IRFinder paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1184-4), the significant intron retention events were defined as:

- `p.adj` < 0.05
- `IRratio` for either control or treatment > 0.1 (at least 10% of transcripts exhibit retention of the intron for at least one of the conditions)
- Retained introns exhibit a coverage (`IntronDepth`) of at least three reads across the entire intron after excluding non-measurable intronic regions

All results can be accessed using the links below the table. The top 20 most significant retained introns are displayed below:

```{r filtering_results}
# Only returning those introns represented in more than 10% of transcripts in A or B
sig_irfinder_filtered <- sig_irfinder_ACtest[which(sig_irfinder_ACtest$A.IRratio > 0.1 | sig_irfinder_ACtest$B.IRratio > 0.1), ]

# Only returning those introns with a coverage of more than three reads across the entire intron after excluding non-measurable intronic regions
sig_irfinder_filtered2 <- sig_irfinder_filtered[sig_irfinder_filtered$A.IntronDepth > 3 | sig_irfinder_filtered$B.IntronDepth > 3, ]

sig_irfinder_filtered2 <- sig_irfinder_filtered2[order(sig_irfinder_filtered2$padj), ]
sig_irfinder_filtered2$p.diff <- formatC(sig_irfinder_filtered2$p.diff, format = "e", digits = 2)
sig_irfinder_filtered2$p.increased <- formatC(sig_irfinder_filtered2$p.increased, format = "e", digits = 2)
sig_irfinder_filtered2$p.decreased <- formatC(sig_irfinder_filtered2$p.decreased, format = "e", digits = 2)
knitr::kable(sig_irfinder_filtered2[1:20, ])
write.csv(sig_irfinder_filtered2, "results/significant_irfinder_ACtest_results_padj.csv")
```

[Download all results](./results/irfinder_ACtest_all_results_padj.csv)

[Download significant results](./results/significant_irfinder_ACtest_results_padj.csv)

The significantly retained introns were explored for several genes of interest. 

```{r specific_genes}
interesting_genes <- sig_irfinder_filtered2[sig_irfinder_filtered2$gene %in% c(), ]
knitr::kable(interesting_genes)

```

```{r sessioninfo}
sessionInfo()
```
