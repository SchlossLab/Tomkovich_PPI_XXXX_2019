Download the [latest release](https://github.com/SchlossLab/new_project/releases/latest) to the directory and decompress


## The proton pump inhibitor omeprazole does not promote *Clostridium difficile* colonization in a murine model

Proton pump inhibitor (PPI) use has been associated with microbiota alterations and susceptibility to *Clostridium difficile* infections (CDIs) in humans. We assessed how PPI treatment alters the fecal microbiota and whether treatment promotes CDIs in a mouse model. Mice receiving a PPI treatment were gavaged with 40 mg/kg of omeprazole during a 7-day pretreatment phase, the day of *C. difficile* challenge, and the following 9 days. We found that mice treated with omeprazole were not colonized by *C. difficile*. When omeprazole treatment was combined with a single clindamycin treatment, one cage of mice remained resistant to *C. difficile* colonization, while the other cage was colonized. Treating mice with only clindamycin followed by challenge resulted in *C. difficile* colonization. 16S rRNA gene sequencing analysis revealed that omeprazole had minimal impact on the structure of the murine microbiota throughout the 16 days of omeprazole exposure. These results suggest omeprazole treatment alone is not sufficient to disrupt microbiota resistance to *C. difficile* infection in mice that are normally resistant in the absence of antibiotic treatment. 


### Overview

	project
	|- README          # the top level description of content (this doc)
	|- CONTRIBUTING    # instructions for how to contribute to your project
	|- LICENSE         # the license for this project
	|
	|- submission/
	| |- study.Rmd    # executable Rmarkdown for this study, if applicable
	| |- study.md     # Markdown (GitHub) version of the *.Rmd file
	| |- study.tex    # TeX version of *.Rmd file
	| |- study.pdf    # PDF version of *.Rmd file
	| |- header.tex   # LaTeX header file to format pdf version of manuscript
	| |- references.bib # BibTeX formatted references
	| |- XXXX.csl     # csl file to format references for journal XXX
	|
	|- data           # raw and primary data, are not changed once created
	| |- references/  # reference files to be used in analysis
	| |- raw/         # raw data, will not be altered
	| |- mothur/      # mothur processed data
	| +- process/     # cleaned data, will not be altered once created;
	|                 # will be committed to repo
	|
	|- code/          # any programmatic code
	|
	|- results        # all output from workflows and analyses
	| |- tables/      # text version of tables to be rendered with kable in R
	| |- figures/     # graphs, likely designated for manuscript figures
	| +- pictures/    # diagrams, images, and other non-graph graphics
	|
	|- exploratory/   # exploratory data analysis for study
	| |- notebook/    # preliminary analyses
	| +- scratch/     # temporary files that can be safely deleted or lost
	|


### How to regenerate this repository

#### Dependencies and locations
* Gnu Make (v4.2) should be located in the user's PATH
* mothur (v1.40.5) should be located in the user's PATH
* R (v. 3.5.1) should be located in the user's PATH
* R packages:
  * knitr v1.2
  * rmarkdown v1.10
  * readxl v1.1.0
  * tidyverse v1.2.1
  * cowplot v.0.9.3
  * magick v.2.0
* Analysis assumes the use of 8 processors


#### Running analysis
Download 16S rRNA sequencing dataset from the NCBI Sequence Read Archive (Accession no. PRJNA554866).
```
git clone https://github.com/SchlossLab/Tomkovich_PPI_mSphere_2019
```
Transfer 16S rRNA sequencing fastq files into Tomkovich_PPI_mSphere_2019/data/raw.
```
cd Tomkovich_PPI_mSphere_2019
```
Obtain the SILVA reference alignment from version 132 described at http://blog.mothur.org/2018/01/10/SILVA-v132-reference-files/. We will use the SEED v. 132, which contain 12,083 bacterial sequences. This also contains the reference taxonomy. We will limit the databases to only include bacterial sequences.
```
wget -N https://mothur.org/w/images/7/71/Silva.seed_v132.tgz
tar xvzf Silva.seed_v132.tgz silva.seed_v132.align silva.seed_v132.tax
mothur "#get.lineage(fasta=silva.seed_v132.align, taxonomy=silva.seed_v132.tax, taxon=Bacteria);degap.seqs(fasta=silva.seed_v132.pick.align, processors=8)"
mv silva.seed_v132.pick.align data/references/silva.seed.align
rm Silva.seed_v132.tgz silva.seed_v132.*
#Narrow to v4 region
mothur "#pcr.seqs(fasta=data/references/silva.seed.align, start=11894, end=25319, keepdots=F, processors=8)"
mv data/references/silva.seed.pcr.align data/references/silva.v4.align
```
Obtain the RDP reference taxonomy. The current version is v11.5 and we use a "special" pds version of the database files, which are described at http://blog.mothur.org/2017/03/15/RDP-v16-reference_files/.
```

wget -N https://www.mothur.org/w/images/c/c3/Trainset16_022016.pds.tgz
tar xvzf Trainset16_022016.pds.tgz trainset16_022016.pds
mv trainset16_022016.pds/* data/references/
rm -rf trainset16_022016.pds
rm Trainset16_022016.pds.tgz
```
Obtain the Zymo mock community data; note that Zymo named the 5 operon of Salmonella twice instead of the 7 operon.
```
wget -N https://s3.amazonaws.com/zymo-files/BioPool/ZymoBIOMICS.STD.refseq.v2.zip
unzip ZymoBIOMICS.STD.refseq.v2.zip
rm ZymoBIOMICS.STD.refseq.v2/ssrRNAs/*itochondria_ssrRNA.fasta #V4 primers don't come close to annealing to these
cat ZymoBIOMICS.STD.refseq.v2/ssrRNAs/*fasta > zymo_temp.fasta
sed '0,/Salmonella_enterica_16S_5/{s/Salmonella_enterica_16S_5/Salmonella_enterica_16S_7/}' zymo_temp.fasta > zymo.fasta
mothur "#align.seqs(fasta=zymo.fasta, reference=data/references/silva.v4.align, processors=12)"
mv zymo.align data/references/zymo_mock.align
rm -rf zymo* ZymoBIOMICS.STD.refseq.v2* zymo_temp.fasta
```
Run the 16S rRNA sequence data through mothur.
Generate a fasta, taxonomy and count-table file that has the chimeras and non bacterial sequences removed.
```
mothur code/get_good_seqs.batch
rm data/mothur/*.map
```
Generate a shared file and a cons.taxonomy file based on the OTU data.
```

mothur code/get_shared_otus.batch
        rm data/mothur/ppi.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.pick.count_table
        rm data/mothur/ppi.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.fasta
        rm data/mothur/ppi.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.pick.taxonomy
```
Calculate overall error rate by comparing Mock sequences to Mock reference community from Zymo.
```
mothur code/get_error.batch
```
Alpha and Beta diversity analysis:
```
mothur code/ppi_alpha_beta.batch
```

Run dist.shared and generate PCoA for 1B using just a subset of samples (Omeprazole treated mice before *C. difficile* challenge, see code/ppi_taxa.R file for list of samples).
```
mothur
set.dir(input=data/mothur, output=data/mothur/subset)
# Run dist.shared on a subet of samples, specified with groups=
dist.shared(shared=ppi.opti_mcc.shared, calc=thetayc-jclass-braycurtis, subsample=3000, groups=Opos+M1+D0-Opos+M2+D0-Opos+M3+D0-Opos+M4+D0-Opos+M5+D0-Opos+M1+D2-Opos+M2+D2-Opos+M3+D2-Opos+M4+D2-Cpos+M5+D-Opos+M1+D4-Opos+M3+D4-Opos+M4+D4-Opos+M5+D4-Opos+M11+D6-Opos+M2+D6-Opos+M3+D-Opos+M4+D6-Opos+M5+D6-Opos+M1+D7-Opos+M2+D7-Opos+M3+D7-Opos+M4+D7-Opos+M5+D7)
# Generate PCoA based on Bray-Curtis distance
set.dir(input=data/mothur/subset, output=data/mothur/subset)
pcoa(phylip=ppi.opti_mcc.braycurtis.0.03.lt.ave.dist)
```
#### Generate statisical analysis and plots for the paper
Run the following scripts in R to generate statisitical analysis and components of the figures.
* ppi_pcoa.R
* ppi_taxa.R
* ppi_cfu_over_time.R
* ppi_alpha.R
Run the figure scripts in R to generate figures for the paper.
* figure_1.R
* figure_2.R
* figure_s1.R
* figure_s2.R

#### Generate the paper
open submission/manuscript.Rmd and knit to Word or pdf file.