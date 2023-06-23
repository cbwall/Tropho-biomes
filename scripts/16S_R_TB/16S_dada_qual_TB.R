######### First in a series of scripts for R
# before you run make sure you close conda 
# in terminal: conda info --env 
# close open conda with 'conda deactivate'


#### load packages

#### load dada2 package
# install an older version of devtools from CRAN
dev<- "https://cran.r-project.org/src/contrib/Archive/devtools/devtools_2.4.3.tar.gz"
install.packages(dev, repos=NULL, type="source")
library("devtools")
library("dada2")
library("dplyr")
# (.packages()) to check if packages are loaded

start_time <- Sys.time() # track timing


#########
setwd('/projects/ps-shurinlab/users/cbwall/Trophobiomes')

# read in metadata
md<- read.csv("data/TB_metadata.csv")
  
# made a sample or control factor
md$sample_control<- md$Sample.type
md$sample_control[md$sample_control=='zooplankton' |
                    md$sample_control=='water_TB']<- "sample"

# ID the - controls
md$sample_control[md$sample_control=='negative.control'] <- "neg.controls"


#### add this SampleNames to the metadata file
S1<-sapply(strsplit(md$UCSD_metanames, "_"), `[`, 2)
S2<-sapply(strsplit(md$UCSD_metanames, "_"), `[`, 3)
S.name<-sampleNames.md<-paste(S1,S2)
sampleNames<-gsub(" ", "_", S.name) # remove space and add an underscore

md$sampleNames<-(sampleNames)
md$Miseq.ANL<-3

# rearrange
run.metaD<- md %>% 
  dplyr::select(UCSD_metanames, sampleNames, Project, Location, Experiment, Time, Tank, 
                Inoc.plank.nutr.source, Inoc.treatment, Zoop.tube_Inoc.tube.ID, Inoc.Bottle.color,
                Inoc.Bottle.number, rep, Sample.type,
                Organism, Functional.group, Number.of.individuals.or.ml, 
                Qubit.DNA..ng.ul, sample_control, Miseq.ANL)

make.fac<-c("Location", "Experiment", "Time", "Tank", 
            "Inoc.plank.nutr.source", "Inoc.treatment", "Zoop.tube_Inoc.tube.ID", "Inoc.Bottle.color",
            "Inoc.Bottle.number", "rep",  "Sample.type", "Organism", 
            "Functional.group", "sample_control", "Miseq.ANL")

run.metaD[make.fac] <- lapply(run.metaD[make.fac], factor) # make all these factors


# load in the cut-adapt samples in the "trimmed" folder
# run this in terminal
# gunzip data/ANL.sequences/*.fastq.gz
# now unzipped... proceed

# path to folder containing demultiplexed library sequencing files
# make sure unzipped
miseq_path<-"data/ANL.sequences" # CHANGE to the directory containing the fastq files after unzipping.
list.files(miseq_path)


################################# Filter and Trim
### remove low quality reads, trim to consistent length

# Sort ensures forward/reverse reads are in same order
fnFs <- sort(list.files(miseq_path, pattern="_R1_001.fastq.gz"))
fnRs <- sort(list.files(miseq_path, pattern="_R2_001.fastq.gz"))

##### ##### ##### ##### 
##### had issue here, the way #s reading in not in order with metadata sheet...
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sampleNames.p2 <- sapply(strsplit(fnFs, "_"), `[`, 2) # extract sample names
sampleNames.p3 <- sapply(strsplit(fnFs, "_"), `[`, 3) # extract the run # sample
sampleNames<-paste(sampleNames.p2,sampleNames.p3) # compile
sampleNames<-gsub(" ", "_", sampleNames) # remove space and add an underscore

SN.FQ<-as.data.frame(sampleNames)

# if need to merge by a column, say if sequences not all in a single run or separated for some reason...
run.metaD.run1 <- merge(run.metaD, SN.FQ, by="sampleNames")

write.csv(run.metaD.run1, "output/run.metaD.edit.csv")

################################ Specify the full path to the fnFs and fnRs
fnFs <- file.path(miseq_path, fnFs)
fnRs <- file.path(miseq_path, fnRs)
#fnFs[1:3]

########################################
# create pdf of quality profiles for forward samples
# the aggregate = TRUE gives quality plot for ALL reads

pdf("output/qual_profiles_F.pdf")
plotQualityProfile(fnFs, aggregate = TRUE)
dev.off()

# create pdf of quality profiles for reverse samples
pdf("output/qual_profiles_R.pdf")
plotQualityProfile(fnRs, aggregate = TRUE)
dev.off()

# final output is two pdfs of quality score profiles - forward and reverse
#################################

Sys.time() - start_time

