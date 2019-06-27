# hbcABC

Compilation of rmarkdown templates for bioinformatics analysis

These templates are not ready for fully automatization,
it gives a hint about the main figures/tables we produce
for specific analysis like chipseq, small-rna, complex rnaseq ... etc

Probably all of them will need some work in some chunks, so don't 
press knitr button directly without testing the sections.

## General tips:

Please go to [general_start](http://htmlpreview.github.com/?blob/master/docs/articles/general_start.html) to get regular tips for all templates.

## Templates

* legacy fastrna differential expression: this works with fastrna pipeline using the legacy RNAseq template from Rory's code.

* Complex differential expression: This template support multiple complex or
time serie data. And uses `degPattern` to cluster genes into expression
profiles. It is required to use `bcbioRNASeq` to create the object.

* Transcript differential expression with sleuth

* Singel cell templates: Quick QC, Seurat Clustering analysis with clustree plots.
