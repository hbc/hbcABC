# hbcABC
Compilation of rmarkdown templates for bioinformatics analysis

These templates are not ready for fully automatization,
it gives a hint about the main figures/tables we produce
for specific analysis like chipseq, small-rna, complex rnaseq ... etc

Probably all of them will need some work in some chunks, so don't 
press knitr button directly without testing the sections.


## Templates

* Complex differential expression: This template support multiple complex or
time serie data. And uses `degPattern` to cluster genes into expression
profiles. It is required to use `bcbioRNASeq` to create the object.
