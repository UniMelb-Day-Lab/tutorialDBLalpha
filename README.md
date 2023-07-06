# tutorialDBLalpha

# DBLalpha TUTORIAL: Clean, Cluster, Classify, and Combine

version v1.0.0

Bioinformatic pipelines developed by the Day Lab at The University of Melbourne for analysing targeted amplicon sequencing for a defined region of the DBL&alpha; domains of the antigen-encoding *P. falciparum* *var* genes.

This tutorial provides a step by step example how to use these bioinformatic pipelines to: 
1. **Clean** raw Illumina sequence data to generate DBL&alpha; tags per isolate ;
2. **Cluster** the cleaned DBL&alpha; tags to define unique DBL&alpha; types;
3. **Classify** the unique DBL&alpha; types to their most probable DBL&alpha; domain class;
4. **Combine** the cluster and classify output files for the DBL&alpha; type analyses

> NOTE: For each step in this tutorial the input and output files are provided in the data folder. The expected run time for the tutorial is ~30 - 45 mins. 

<br />


> ## **PCR Amplification and Amplicon Library Preparation for Illumina Sequencing**
> The DBL&alpha; domain of *P. falciparum* *var* genes were amplified from genomic DNA using fusion primers for multiplexed sequencing. 
> 
> We coupled template-specific universal degenerate primer sequences (i.e.,DBL&alpha;AF and DBL&alpha;BR) ([Taylor et al. (2000)](https://www.sciencedirect.com/science/article/pii/S0166685199001590?via%3Dihub) and [Bull et al. (2005)](https://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.0010026)) with a 454 Titanium primer sequence (i,e., adaptor) and a unique 10 bp multiplex identifier (MID) tag. For additional details see [Rask et al. (2016)](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1032-7). For Illumina sequencing, the amplicon library pools were prepared by combining equimolar amounts of the PCR amplicons each with a unique MID tags. This approach was used so that reads could be demultiplexed into individual files for each isolate and then paired based on valid combinations of MID tags (i.e., exact match) in the forward and reverse reads. 
> 
> Example forward primer: **DBL&alpha;AF-MID-1**:     5'-(adaptor)(MID)(DBL&alpha; AF forward primer)-3'
>```
> 5'-(CGTATCGCCTCCCTCGCGCCATCAG)(ACGAGTGCGT)(GCACGMAGTTTYGC)-3'
>```
> Example reverse primer: **DBL&alpha;BR-MID-1**:     5'-(adaptor)(MID)(DBL&alpha; BR reverse primer)-3'
>```
> 5'-(CTATGCGCCTTGCCAGCCCGCTCAG)(ACGAGTGCGT)(GCCCATTCSTCGAACCA)-3'
>```
> The use of 454 adaptor sequences here is a remnant of previous work and is being phased out. They are not necessary for this pipeline, 
> as [DBLaCleaner](https://github.com/UniMelb-Day-Lab/DBLaCleaner) only requires that the sequences contain the forward or reverse primer sequence (i.e.,DBL&alpha;AF and DBL&alpha;BR) ([Taylor et al. (2000)](https://www.sciencedirect.com/science/article/pii/S0166685199001590?via%3Dihub) and [Bull et al. (2005)](https://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.0010026), along with the 
> unique MID tags for each isolate.

<br />

<img src="https://github.com/UniMelb-Day-Lab/tutorialDBLalpha/blob/main/figures/clean%20cluster%20classify%20figure.png" width="400">

<br />

## Step 1: Clean 
The raw Illumina sequence data are first cleaned and processed using [DBLaCleaner](https://github.com/UniMelb-Day-Lab/DBLaCleaner). 
Please visit the [DBLaCleaner](https://github.com/UniMelb-Day-Lab/DBLaCleaner) repository for detailed instructions on installation and usage and [He et al.(2018)](https://www.nature.com/articles/s41467-018-04219-3) for an indepth description of the method.

<br />

Example Input Files Step 1: 
1. First read file (-r READ1): *pool_1_R1.fastq*
2. Paired read file (-R READ2): *pool_1_R2.fastq*
3. Mapping file (-d DESC): *pool_1_mapping_file.desc*

> Note: There are examples of raw sequence data from three sequencing pools in the data folder for Step 1 (pool_1, pool_2, pool_3). All three sequencing pools need
> to be processed independently through **DBLaCleaner**.

Command line: Example Pool 1
```
python2.7 cleanDBLalpha.py -o ./ -r pool_1_R1.fastq -R pool_1_R2.fastq -d pool_1_mapping_file.desc --perID 0.96 --cpu 2 --verbose
```

Command line: Example Pool 2
```
python2.7 cleanDBLalpha.py -o ./ -r pool_2_R1.fastq -R pool_2_R2.fastq -d pool_2_mapping_file.desc --perID 0.96 --cpu 2 --verbose
```

Command line: Example Pool 3
```
python2.7 cleanDBLalpha.py -o ./ -r pool_3_R1.fastq -R pool_3_R2.fastq -d pool_3_mapping_file.desc --perID 0.96 --cpu 2 --verbose
```


Example Output Files Step 1:
This step generates a fasta file with all the cleaned DBL&alpha; tags. 
1. Cleaned DBL&alpha; tags: *pool_1_R_DBLa_cleaned.fasta*
2. Summary statistics: *pool_1_summaryStatistics.log*
3. Contaminant file: *pool_1_R_NOT_dblalpha.fasta*

> These cleaned DBLa tags are then used in Step 2. 

<br />
<br />
<br />

## Step 2: Cluster
Next, the cleaned DBL&alpha; tags are clustered using [clusterDBLalpha](https://github.com/UniMelb-Day-Lab/clusterDBLalpha) to define the unique DBL&alpha; types. 
Please visit the [clusterDBLalpha](https://github.com/UniMelb-Day-Lab/clusterDBLalpha) repository for detailed instructions on installation and usage and [[He et al.(2018)](https://www.nature.com/articles/s41467-018-04219-3) for an indepth description of the method.


### **Step 2A:** Before the DBL&alpha; tags can be clustered, these need to be combined into one fasta file. 

Example Input Files Step 2A:
1. Cleaned DBL&alpha; tags: *pool_1_R_DBLa_cleaned.fasta*
2. Cleaned DBL&alpha; tags: *pool_2_R_DBLa_cleaned.fasta*
3. Cleaned DBL&alpha; tags: *pool_3_R_DBLa_cleaned.fasta*


Command line 
```
cat pool_1_R_DBLa_cleaned.fasta pool_2_R_DBLa_cleaned.fasta pool_3_R_DBLa_cleaned.fasta > pool_combine123_DBLa_cleaned.fasta
```

Example Output Files Step 2A:
1. Combined DBL&alpha; tags: *pool_combine123_DBLa_cleaned.fasta*

<br />

### **Step 2B:** The combined DBL&alpha; tags are then clustered. 

Example Input Files Step 2B:
1. Fasta file containing combined DBL&alpha; tags (-r READ): *pool_combine123_DBLa_cleaned.fasta*


Command line 
```
python2.7 clusterDBLa.py -o ./ -r ./pool_combine123_DBLa_cleaned.fasta --perID 0.96 --cpu 2 --verbose
```

Example Output Files Step 2B:
This step generates several files. Important files include: 
1.	Unique DBL&alpha; types (binary OTU Table): *pool_combine123_DBLa_cleaned_renamed_otuTable_binary.txt* 
> Note: This file is a binary operational taxonomic unit (OTU) table where each row represents a **unique DBL&alpha; type** and each column represents a ***P. falciparum* isolate**.
2.	Unique DBL&alpha; type sequences: *pool_combine123_DBLa_cleaned_renamed_centroids.fasta*
> Note: These unique DBL&alpha; type sequences are then further classified in Step 3. 
> 
<br />
<br />
<br />

## Step 3: Classify
Following Step 2, the unique DBL&alpha; types are further classified using [classifyDBLalpha](https://github.com/UniMelb-Day-Lab/classifyDBLalpha). In this step the unique DBL&alpha; types are assigned to their most probable DBL&alpha; domain class (i.e., DBL&alpha;0, DBL&alpha;1, or DBL&alpha;2) based on E-values, and classified as either upsA (i.e., DBL&alpha;1), non-upsA (i.e., DBL&alpha;0 or DBL&alpha;2), or other (i.e., putatively non-DBLa like).  

Please visit the [classifyDBLalpha](https://github.com/UniMelb-Day-Lab/classifyDBLalpha) repository for detailed instructions on installation and usage and [Ruybal-Pesántez et al. (2017)](https://www.nature.com/articles/s41598-017-11814-9) for an indepth description of the method. 

<br />

Example Input Files Step 3:
1. Fasta file containing unique DBL&alpha; type sequences (-r READ): *pool_combine123_DBLa_cleaned_renamed_centroids.fasta*

Command line
```
python2.7 allocate_reads.py -r ./pool_combine123_DBLa_cleaned_renamed_centroids.fasta -E 1e-8 -o ./ --noUproc --splitIsolates
```
Command line
```
python2.7 reads_to_domains.py --hmm pool_combine123_DBLa_cleaned_renamed_centroids_nhmmOut.txt --out pool_combine123_DBLa_reads_to_domains.csv 
```

Example Output Files Step 3:
1.  Text file containng assigned DBL&alpha; domain classes: *pool_combine123_DBLa_reads_to_domains.csv* 

<br />
<br />
<br />

## Step 4: Combine
Finally, the output files from:
1. Step 2 (*pool_combine123_DBLa_cleaned_renamed_otuTable_binary.txt*)
2. Step 3 (*pool_combine123_DBLa_reads_to_domains.csv*) 

are combined to generate a final DBL&alpha; type dataset (*pool_combine123_DBLa_binary_ups.csv*) for analysis as required. Please see *tutorialDBLa_combine.Rmd* or below for more details.  

```
## Load data
binary_example <-fread("~/Desktop/pool_combine123_DBLa_cleaned_renamed_otuTable_binary.txt", data.table = FALSE)
ups_example <-fread("~/Desktop/pool_combine123_DBLa_reads_to_domains.csv", data.table = FALSE)

## Binary Example: Rename OTU to DBLa type
binary_example<-binary_example%>% rename_at("#OTU ID", ~"DBLa_type")

## Ups Example: Categorize the DBLa domains as upsA, non-upsA, or other
ups_example$read <- gsub(";.*", "", ups_example$read)
ups_example  <- ups_example  %>% rename_at("read", ~"DBLa_type")
ups_example$domain <- as.factor(ups_example$domain)

ups_example$Ups <- ifelse(grepl("^DBLa1", ups_example$domain, ignore.case = T), "upsA",
                    ifelse(grepl("^DBLa0", ups_example$domain, ignore.case = T), "non-upsA", 
                    ifelse(grepl("^DBLa2", ups_example$domain, ignore.case = T), "non-upsA", "Other")))
ups_example$Ups <- as.factor(ups_example$Ups)


## FINAL Combined: Combine ups grouping to the DBLa types in the isolate_example. This file used for the combined DBLa type analyses.
isolate_ups_example<-binary_example %>% inner_join(dplyr::select(ups_example, DBLa_type, Ups), by ="DBLa_type")
isolate_ups_example <- isolate_ups_example %>% relocate(Ups, .after = DBLa_type)

write.table(isolate_ups_example, file="pool_combine123_DBLa_binary_ups.csv", sep=",")
```
<br />
<br />
<br />

## References
Please use the following citations, as indicated, when referencing these bioinformatic pipelines:

1. He, Q., Pilosof, S., Tiedje, K. E., Ruybal-Pesántez, S., Artzy-Randrup, Y., Baskerville, E. B., … Pascual, M. (2018). Networks of genetic similarity reveal non-neutral processes shape strain structure in Plasmodium falciparum. Nature Communications, 9(1), 1817.(https://www.nature.com/articles/s41467-018-04219-3) (Reference for: **DBLaCleaner** and **clusterDBLalpha**)

2. Ruybal-Pesántez, S., Tiedje, K. E., Tonkin-Hill, G., Rask, T. S., Kamya, M. R., Greenhouse, B., … Day, K. P. (2017). Population genomics of virulence genes of Plasmodium falciparum in clinical isolates from Uganda. Scientific Reports, 7(1), 11810.(https://www.nature.com/articles/s41598-017-11814-9) (Reference for: **classifyDBLalpha**)
