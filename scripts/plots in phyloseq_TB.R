########### ########### investigate the final 16S and 18S library ########### ########### 

# load packages
if (!require("pacman")) install.packages("pacman") # for rapid install if not in library

devtools::install_github("benjjneb/dada2", ref="v1.20") # update to most recent dada2
devtools::install_github("zdk123/SpiecEasi")
remotes::install_github("microbiome/microbiome")


# use pacman to load CRAN packages missing
pacman::p_load('knitr', 'microbiome', 'phyloseq', 'tidyr', 'tidyverse', 'knitr', 'magrittr', 'effects', 'devtools',
               'stringi', 'dplyr', "ggplot2", "gridExtra", "dada2", "phyloseq", "vegan", "cowplot", 'doBy', 'ecodist',
               'glue', 'geosphere', 'data.table', 'patchwork', 'car', 'ggcorrplot', 'FactoMineR', 'devtools', 'reshape',
               'lattice',  'plyr', 'magrittr', 'factoextra', 'multcompView', 'decontam', 'factoextra', 'car', 'mia', 
               'ade4', 'fossil', 'picante', 'reshape', 'readr', 'corrr', 'Hmisc', 
               "decontam","BiocManager", 'ggpubr', 'ggmap', "ggordiplots", "fossil", "SpiecEasi", "igraph", "huge", "MASS")


############################ What are all the ID column meanings?
#### Water microbiome

# T0-prestock:  "background water" before any plankton, plants, or sticks added; nutrients had been added
# T0.5 = 7/4/22
# T1 = 7/19/22
# T2 = 7/29/22
# Inoc. Exp. start = 8/5/22, this was the water at the start of the inoculation experiment
# T3 = 8/17/22, end of the experiment



#### Bottle experiment

# T0, POOLED STOCK = pooled stock of plankton that went into the bottles for the bottle experiment (pooled from all +plankton tanks -- this is what went into all the bottles). The IDs are L+, H+, LP+ from 8/3. Just daphnia.
# T1 = *plankton samples taken in each +plankton treatment tank, represent the plankton community in each tank. This was from the un-pooled stock, later pooled  as the stock for plankton that went into microbiome bottle experiemnt. Same as on 8/3 but this is the tank-specific samples. The tube IDs are "tank number", "1" (for T1), "MB" for microbiome. ie., 1.1_MB, 2.1_MB. Just daphnia.
# T3 = * plankton in the tanks at the end of the experiment. Tank zoop communities, these are tubes with the picked community of plankton from each tank. Just daphnia.
# zoop.bottles = the zooplankton at the end of the bottle incubation.


################### ################### ################### ################### 
################### start with the final PS object out from dada2

# ps.prune is the final phyloseq object from pipeline, can load it in here...
### or load in the raw data and re-assemble

ps.prune<-readRDS("output/TSCC_output/ps.prune.RDS")
# export metadata and play
metaD<-microbiome::meta(ps.prune)

metaD$Experiment<- recode_factor(metaD$Experiment, "zoops bottles" = "zoop bottles")
metaD$Experiment<-factor(metaD$Experiment, levels= c("mesocosm water", "zoop mesocosms", "zoop bottles"))

metaD$Inoc.treatment<-factor(metaD$Inoc.treatment, levels= c("Low", "L+", "L-", "High", "H+", "H-", "Plant", "LP+", "LP-", "stream"))

metaD$Time<- recode_factor(metaD$Time, "Inoc. Exp. start" = "Inoc.Exp.start")
metaD$Time<-factor(metaD$Time, levels= c("T0-prestock", "T0.5", "T1", "T2", "T3", "Inoc.Exp.start"))


# to replace an NA in the factor level here...
facna <- addNA(metaD$Inoc.plank.nutr.source)
levels(facna) <- c(levels(metaD$Inoc.plank.nutr.source), "stream")
metaD$Inoc.plank.nutr.source<-facna
metaD$Inoc.plank.nutr.source<-factor(metaD$Inoc.plank.nutr.source, levels= c("Low", "L+", "High", "H+", "Plant.low", "LP+", "stream"))
metaD$Inoc.plank.nutr.source<- recode_factor(metaD$Inoc.plank.nutr.source, "Plant.low" = "Plant")

# to replace an NA in the factor level here...
facna2 <- addNA(metaD$Inoc.treatment)
levels(facna2) <- c(levels(metaD$Inoc.treatment), "stream")
metaD$Inoc.treatment<-facna2

# make a new factor to reconcile some of the name
metaD$Meso.trt<-metaD$Inoc.Bottle.color
metaD$Meso.trt<- recode_factor(metaD$Meso.trt, 
                               "blue" = "bot.H-",
                               "green" = "bot.H+",
                               "orange" = "bot.LP-",
                               "red"= "bot.L+",
                               "white" = "bot.LP+",
                               "yellow"= "bot.L-")

metaD$Meso.trt<-factor(metaD$Meso.trt, levels= 
                         c("L+", "L-", "H+", "H-", "LP+", "LP-",
                         "bot.L+", "bot.L-", "bot.H+", "bot.H-", "bot.LP+", "bot.LP-", "stream"))
                         
# Make a data frame with a column for the read counts of each sample
sample_sum_df<-as.data.frame(sample_sums(ps.prune))
colnames(sample_sum_df)<-"read.sum"
sample_sum_df$sampleNames<-rownames(sample_sum_df)

# are row names the same? if so, add them in
identical(rownames(sample_sum_df), rownames(metaD))
metaD$read.sum<-sample_sum_df$read.sum

sample_data(ps.prune)<-metaD

ps.fin<-ps.prune # rename and export
saveRDS(ps.fin, "output/TSCC_output/ps.fin.RDS")
########## ########## ########## ########## 

# END METADATA SMITHING until later with final PS object (rare-hell)

########## ########## ########## ########## 

########## EXPLORE PLOTS
##### 16S water and zoops
#### inspect plots
Inoc.pl.nut.src.reads<-ggplot(metaD, aes(x=Experiment, y=read.sum, color=Meso.trt)) + geom_boxplot() +
  theme_bw() + ylim(0,60000)+
  ggtitle("16S: Reads by TRT")

Sample.type.reads<-ggplot(metaD, aes(x=Experiment, y=read.sum, color=Sample.type)) + geom_boxplot() +
  theme_bw() + ylim(0,60000)+
  ggtitle("16S: Reads by year")

# subset just zoops to see if read sum effected by $ of individuals
zoop.df<-metaD[(metaD$Sample.type=="zooplankton"),]
zoop.df$Number.of.individuals.or.ml<- as.numeric(zoop.df$Number.of.individuals.or.ml)
mean(zoop.df$read.sum) # = 29,000
less.30<-zoop.df[(zoop.df$Number.of.individuals.or.ml<30),]; mean(less.30$read.sum); sd(less.30$read.sum) # = 17,341 +/- 14k
more.30<-zoop.df[(zoop.df$Number.of.individuals.or.ml>=30),]; mean(more.30$read.sum); sd(more.30$read.sum) # 30,000 +/- 14k
# plot
ggplot(zoop.df, aes(y = read.sum, x = Number.of.individuals.or.ml)) + geom_point()


#############
# Histogram of sample read counts 16S
hist.depth.16S <- ggplot(metaD, aes(read.sum)) + 
  geom_histogram(color="gold4", fill= "gold2", binwidth = 2000) +
  xlab("Read counts") + ylim(0,16)+
  ggtitle("16S: Sequencing Depth") + theme_classic()

# a different view
read_depth_violin.16S <- ggplot(metaD, aes(x=1, y=read.sum)) +
  geom_violin(color="gold4", fill= "gold2", alpha=0.2) + scale_y_log10(limits = c(10e2,10e4)) +
  xlab("Samples") + ylab("Read Depth") + geom_boxplot(width=0.1)+ theme_classic() +
  ggtitle("16S: Distribution of sample sequencing depth")


####################
# Rarefaction curve 16S
count_table_filt.16S <- as.data.frame(otu_table(ps.fin))

pdf(file="figures/rarefac.16S.pdf", height=6, width=8)
rarecurve(count_table_filt.16S, col="gold4", step=50, cex=0.5, ylim=c(0,1500), xlim=c(0,6000), ylab = "16S ASVs", label=FALSE)
#abline(v = 500, lty = "dotted", col="red", lwd=2)
dev.off()


########## combine and export plots ################
inspect.raw.16S<- plot_grid(Inoc.pl.nut.src.reads, Sample.type.reads, hist.depth.16S, read_depth_violin.16S,
                                ncol=4, rel_widths = c(7,5,3,3))


inspect.raw.16S
dev.copy(pdf, "figures/inspect.raw.16S.pdf", height=5, width=17)
dev.off()



################### 16S ################### ################### ################### 

###### subset the PS objects for needs...
# final and NOT rarefied data: 
ps.fin

# remove samples with < 5000 reads, 219 samples, if not in 1%, remove
ps.fin.og <- prune_samples(sample_sums(ps.fin) > 5000, ps.fin)
ps.fin.og<- prune_taxa(taxa_sums(ps.fin.og) > 2, ps.fin.og) 
ps.fin.og # = 219 samples, 7936 taxa

# rarified to 6500 reads, which is about 1.1*min(sample_sums(ps.fin.og))
PS.rar = rarefy_even_depth(ps.fin.og, rngseed=111, sample.size=6500, replace=F) 
table(sample_data(PS.rar)$Sample.type) 
# 219 samples, 6477 taxa
# 262 OTUs were removed
# 99 water, 120 zoops

######################
# rarify, and hellinger
PS.rar.hell = transform_sample_counts(PS.rar, function(x) x^0.5)
table(tax_table(PS.rar.hell)[, "Kingdom"], exclude = NULL)
# 45 Archaea, 6434 Bacteria

###### ###### ###### 
# Compare full data vs. rarefied data beta diversity
####### full data
ord.full <- ordinate(ps.fin.og, method = 'PCoA', distance = 'bray')
plot_ordination(ps.fin.og, ord.full, color = "Sample.type")

# what samples have low reads
plot_ordination(ps.fin.og, ord.full, shape = "Sample.type") +
  geom_point(aes(color = read.sum > 6500))

Rare.culling<-plot_ordination(ps.fin.og, ord.full, color = "Sample.type") +
  geom_point(aes(shape = read.sum > 6500, size= read.sum > 6500)) +
  scale_size_manual(values= c(6,2)) +
  scale_shape_manual(values=c(1, 16))+
  ggtitle("Samples culled by rarefaction") + theme_classic()

# plot the PCoA for non-rarified
Non.rar.NMDS<-plot_ordination(
  physeq = ps.fin.og,                                                   
  ordination = ord.full) +                                                
  geom_point(aes(color = Sample.type), size = 2) +    
  stat_ellipse(level=0.9, linetype = 2, aes(color=Sample.type)) +
  ggtitle("ALL taxa-all data") +
  theme_classic() 

###### if rarefied
ord.rar <- ordinate(PS.rar, method = 'PCoA', distance = 'bray')
plot_ordination(PS.rar, ord.rar, color = "Sample.type")

Rar.NMDS<- plot_ordination(
  physeq = PS.rar,                                                   
  ordination = ord.rar) +                                                
  geom_point(aes(color = Sample.type), size = 2) +    
  stat_ellipse(level=0.9, linetype = 2, aes(color=Sample.type)) +
  ggtitle("ALL taxa- 6500 rar") +
  theme_classic() 

#### rarefied, but hellinger transformed
ord.hell <- ordinate(PS.rar.hell, method = 'PCoA', distance = 'bray')

Hell.NMDS<- plot_ordination(
  physeq=PS.rar.hell, 
  ordination = ord.hell) +
  geom_point(aes(color = Sample.type), size = 2) +    
  stat_ellipse(level=0.9, linetype = 2, aes(color=Sample.type)) +
  ggtitle("ALL taxa - hellinger-6500 rar") +
  theme_classic() 


NMDS.tests<-plot_grid(Rare.culling + guides(color="none"), 
                      Non.rar.NMDS +  guides(color="none"), 
                      Rar.NMDS +  guides(color="none"), 
                      Hell.NMDS, ncol=4, rel_widths=c(4,4,4,5.5))
NMDS.tests
dev.copy(pdf, "figures/16S.NMDS.tests.pdf", height=6, width=20)
dev.off()
#######

#export final metadata
metaD.final<-microbiome::meta(PS.rar.hell)
write.csv(metaD.final, "output/metaD.final.csv")


###########
# subset for bottles or water exp
PS.bottles = subset_samples(PS.rar.hell, Experiment=="zoop bottles") # 40 samples
PS.mesocosm = subset_samples(PS.rar.hell, Experiment=="zoop mesocosms" | Experiment=="mesocosm water")  # 179 samples
PS.water = subset_samples(PS.rar.hell, Experiment=="mesocosm water")  # 99 samples
PS.mesozoop = subset_samples(PS.rar.hell, Experiment=="zoop mesocosms")  # 99 samples
#########


# test ordination
ord.bottles <- ordinate(PS.bottles, method = 'PCoA', distance = 'bray')

# do zooplankton group by the treatments they CAME from
bottle.1<- plot_ordination(
  physeq=PS.bottles, 
  ordination = ord.bottles) +
  scale_color_manual(values=c("lightskyblue", "palegreen3", "lightsalmon2")) +
  geom_point(aes(color = Inoc.plank.nutr.source), size = 2) +    
  stat_ellipse(level=0.9, linetype = 2, aes(color=Inoc.plank.nutr.source)) +
  ggtitle("Zoop source trt") +
  theme_classic() 

# do zooplankton group by the treatments they WENT TO
bottle.2<- plot_ordination(
  physeq=PS.bottles, 
  ordination = ord.bottles) +
  scale_color_manual(values=c("lightskyblue","dodgerblue3", "palegreen3", "forestgreen", "lightsalmon2", "firebrick2")) +
  geom_point(aes(color = Inoc.treatment), size = 2) +    
  stat_ellipse(level=0.9, linetype = 2, aes(color=Inoc.treatment)) +
  ggtitle("Water trt") +
  theme_classic() 

# combine the 2 together with shapes and colors (shape=source, color = treatment water)
bottle.3<- plot_ordination(
  physeq=PS.bottles, 
  ordination = ord.bottles) +
  scale_color_manual(values=c("lightskyblue","dodgerblue3", "palegreen3", "forestgreen", "lightsalmon2", "firebrick2")) +
  geom_point(aes(shape = Inoc.plank.nutr.source, color=Inoc.treatment), size = 3) +    
  stat_ellipse(level=0.9, linetype = 2, aes(color=Inoc.treatment)) +
  ggtitle("Water trt") +
  theme_classic() 


########################
# test ordination
ord.water <- ordinate(PS.water, method = 'PCoA', distance = 'bray')

water.time<- plot_ordination(
  physeq=PS.mesocosm, 
  ordination = ord.water) +
  geom_point(aes(color = Inoc.treatment, shape=Time), size = 2.5) +    
  scale_color_manual(values=c("lightskyblue", "palegreen3", "lightsalmon2", "gray40")) +
  scale_shape_manual(values=c(3,16,17,15,16,1))+
  stat_ellipse(level=0.9, linetype = 2, aes(color=Inoc.treatment)) +
  ggtitle("Zoop source trt") +
  theme_classic() 

water.time$layers<-water.time$layers[-1]
water.time

####################











###### some richness plots
#richness by Lake
richness.plot_inoc.trt<-plot_richness(ps.prune, x="Inoc.treatment", measures=c("Observed", "Shannon"))  + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
richness.plot_inoc.trt
dev.copy(pdf, "figures/rrichness.plot_inoc.trt.pdf", height=4, width=10)
dev.off() 


#richness by organism (or sample type i.e, water)
richness.plot_organism<-plot_richness(ps.prune, x="Inoc.plank.nutr.source", measures=c("Observed", "Shannon")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
richness.plot_organism
dev.copy(pdf, "figures/richness.plot.organism.pdf", height=4, width=12)
dev.off() 


#richness by Time.point
richness.plot_Time<-plot_richness(ps.prune, x="Time", measures=c("Observed", "Shannon")) + theme_bw()
richness.plot_Time
dev.copy(pdf, "figures/richness.plot_Time.pdf", height=4, width=10)
dev.off() 


#richness by Experiment
richness.plot_Exp<-plot_richness(ps.prune, x="Experiment", measures=c("Observed", "Shannon")) + theme_bw()
richness.plot_Exp
dev.copy(pdf, "figures/richness.plot_Exp.pdf", height=4, width=10)
dev.off() 


#richness by Functional Group
richness.plot_Func<-plot_richness(ps.prune, x="Functional.group", measures=c("Observed", "Shannon")) + theme_bw()
richness.plot_Func
dev.copy(pdf, "figures/richness.plot_Func.pdf", height=4, width=10)
dev.off() 
###### 


##### Read counts and rarefaction curves
# Make a data frame with a column for the read counts of each sample
sample_sum_df<-as.data.frame(sample_sums(PS.fin))
colnames(sample_sum_df)<-"read.sum"
sample_sum_df$sampleNames<-rownames(sample_sum_df)


# merge in the reads
run.metaD<-merge(M, sample_sum_df, by="sampleNames", all.y=TRUE)
write.csv(run.metaD, "output/run.metaD.final.csv")


# Histogram of sample read counts
hist.depth<-ggplot(sample_sum_df, aes(x = read.sum)) + 
  geom_histogram(color = "black", fill = "gold2", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank()) + theme_classic() + geom_vline(xintercept=1000, lty=2)

hist.depth
dev.copy(pdf, "figures/hist.depth.pdf", height=4, width=5)
dev.off() 


##################### ###################
##################### #################

OTU <- otu_table(PS.fin)
class(OTU) <- "matrix" # as.matrix() will do nothing
## you get a warning here, but this is what we need to have

OTU <- t(OTU) # transpose observations to rows

# version 1
pdf(file="figures/rare.raw.pdf", height=6, width=10)
rarecurve(OTU, step=50, cex=0.5, xlim=c(0,10000), ylim=(c(0,300)), label=FALSE)
#abline(v = 5000, lty = "dotted", col="red", lwd=2)
dev.off() 



#### more plots
pdf(file="figures/read.by.Inoc.trt.pdf", height=4, width=12)
boxplot(run.metaD$read.sum~run.metaD$Inoc.treatment, cex.axis=0.5, cex.lab=0.8)
dev.off() 

pdf(file="figures/read.by.sampletype.pdf", height=6, width=5)
boxplot(run.metaD$read.sum~run.metaD$Sample.type, cex.axis=0.5, cex.lab=0.8)
dev.off() 

pdf(file="figures/reads.sample.pdf", height=4, width=7)
ggplot(run.metaD, aes(x=sample_control, y=read.sum, color=Functional.group)) + geom_boxplot()
dev.off() 
##### 



################### ################### ###################
################### ################### ###################

### NMDS
############ NMDS
sample_variables(PS.fin)

# make colors for sites
library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))


###
Time.water<- subset_samples(PS.fin, Experiment=="mesocosm water")

Tpre.phy<-subset_samples(Time.water, Time=="T0-prestock")
T0.5.phy<-subset_samples(Time.water, Time=="T0.5")
T1.phy<-subset_samples(Time.water, Time=="T1")
T2.phy<-subset_samples(Time.water, Time=="T2")
T3.phy<-subset_samples(Time.water, Time=="T3")
 
Zoop.TB.phy<- subset_samples(PS.fin, Experiment=="zoop mesocosms")

Zoop.stock.phy<- subset_samples(PS.fin, Tank=="pooled.stock")

Zoop.bottle.phy<- subset_samples(PS.fin, Experiment=="zoops bottles")

###############
ORD.zoop <- ordinate(Zoop.TB.phy, method='MDS', distance='bray')

#### zoops in tanks
# T0= pooled stock
# T1 = plankton in each tank a few days later
# T3 = zoops at end of experiment

NMDS.zoops.tanks<-plot_ordination(
  physeq = Zoop.TB.phy,                                                   
  ordination = ORD.zoop) +                                                
  geom_point(aes(color=Time), size = 3) +    
  stat_ellipse(level=0.9, linetype = 2, aes(color=Time)) +
  ggtitle("TB end tank Daph") + theme_classic()
  
#####
ORD.zoop.stock <- ordinate(Zoop.stock.phy, method='MDS', distance='bray')

#### zoops in tanks, stocks pooled by treatment for bottle incubations
NMDS.zoops.stock<-plot_ordination(
  physeq = Zoop.stock.phy,                                                   
  ordination = ORD.zoop.stock) +                                                
  geom_point(aes(color=Inoc.plank.nutr.source), size = 3) +    
  geom_polygon(aes(fill=Inoc.Bottle.color), alpha=0.5) +
  guides(color = "none")+
  ggtitle("TB Zoop bottle stock") + theme_classic()

NMDS.zoops.stock

###############
ORD.zoop.bot <- ordinate(Zoop.bottle.phy, method='MDS', distance='bray')
  
#### zoops in tanks
NMDS.zoops.bott.1<-plot_ordination(
    physeq = Zoop.bottle.phy,                                                   
    ordination = ORD.zoop.bot) +                                                
    geom_point(aes(color=Inoc.plank.nutr.source), size = 3) +    
    stat_ellipse(level=0.9, linetype = 2, aes(color=Inoc.plank.nutr.source)) +
    ggtitle("TB Zoopbot source") + theme_classic()

NMDS.zoops.bott.2<-plot_ordination(
  physeq = Zoop.bottle.phy,                                                   
  ordination = ORD.zoop.bot) +                                                
  geom_point(aes(color=Inoc.treatment), size = 3) +    
  stat_ellipse(level=0.9, linetype = 2, aes(color=Inoc.treatment)) +
  ggtitle("TB zoops Inoc treatm") + theme_classic()


plot_grid(
  NMDS.zoops.tanks, 
  NMDS.zoops.stock, 
  NMDS.zoops.bott.1, 
  NMDS.zoops.bott.2,
  rel_widths = c(8,8,8,8), ncol=4)

ggsave("figures/TBzoop.NMDS.pdf", height=5, width=13)


###############
ORD.t.water <- ordinate(Time.water, method='MDS', distance='bray')

#### just time change across all tanks
NMDS.ord.time.water<-plot_ordination(
  physeq = Time.water,                                                   
  ordination = ORD.t.water) +                                                
  geom_point(aes(color=Time), size = 3) +    
  stat_ellipse(level=0.9, linetype = 2, aes(color=Time)) +
  theme_classic()   

NMDS.ord.time.water
dev.copy(pdf, "figures/NMDS.ord.time.water.pdf", height=6, width=9)
dev.off() 


### just water source change, pooled across time
NMDS.ord.TRT<-plot_ordination(
  physeq = Time.water,                                                   
  ordination = ORD.t.water) +                                                
  geom_point(aes(color=Inoc.Bottle.color), size = 3) +    
  stat_ellipse(level=0.9, linetype = 2, aes(color=Inoc.Bottle.color)) +
  #scale_color_brewer(palette = "Dark2") +
  #scale_fill_brewer(palette = "Dark2") + ggtitle("Time1") +
  theme_classic()   

NMDS.ord.TRT
dev.copy(pdf, "figures/NMDS.ord.TRT.pdf", height=6, width=9)
dev.off() 

########## ########## 
########## ########## 
########## ########## 

# Tpre

ORD.Tpre <- ordinate(Tpre.phy, method='MDS', distance='bray')


NMDS.Tpre<-plot_ordination(
  physeq = Tpre.phy,                                                   
  ordination = ORD.Tpre) +                                                
  geom_point(aes(color=Inoc.Bottle.color), size = 3) +    
  geom_polygon(aes(fill=Inoc.Bottle.color), alpha=0.5) +
  guides(color = "none")+
  ggtitle("Time-prestock")+
  theme_classic()    

NMDS.Tpre
########## ########## 

# T0.5

ORD.T0.5 <- ordinate(T0.5.phy, method='MDS', distance='bray')

#### just time change across all tanks
NMDS.T0.5<-plot_ordination(
  physeq = T0.5.phy,                                                   
  ordination = ORD.T0.5) +                                                
  geom_point(aes(color=Inoc.Bottle.color), size = 3) +    
  geom_polygon(aes(fill=Inoc.Bottle.color), alpha=0.5) +
  guides(color = "none")+
  ggtitle("Time0.5")+
  theme_classic()    

NMDS.T0.5

### Time 1
ORD.T1 <- ordinate(T1.phy, method='MDS', distance='bray')

NMDS.T1<-plot_ordination(
  physeq = T1.phy,                                                   
  ordination = ORD.T1) +                                                
  geom_point(aes(color=Inoc.Bottle.color), size = 3) +    
  geom_polygon(aes(fill=Inoc.Bottle.color), alpha=0.5) +
  guides(color = "none")+
  ggtitle("Time1")+
  theme_classic()  

NMDS.T1

### Time 2
ORD.T2 <- ordinate(T2.phy, method='MDS', distance='bray')

NMDS.T2<-plot_ordination(
  physeq = T2.phy,                                                   
  ordination = ORD.T2) +                                                
  geom_point(aes(color=Inoc.Bottle.color), size = 3) +    
  geom_polygon(aes(fill=Inoc.Bottle.color), alpha=0.5) +
  guides(color = "none")+
  ggtitle("Time2")+
  theme_classic()    

NMDS.T2

### Time 3
ORD.T3 <- ordinate(T3.phy, method='MDS', distance='bray')

NMDS.T3<-plot_ordination(
  physeq = T3.phy,                                                   
  ordination = ORD.T3) +                                                
  geom_point(aes(color=Inoc.Bottle.color), size = 3) +    
  geom_polygon(aes(fill=Inoc.Bottle.color), alpha=0.5) +
  guides(color = "none")+
  ggtitle("Time3")+
  theme_classic()   

NMDS.T3

####

extract.legend <- get_legend(
  # create some space to the left of the legend
  NMDS.Tpre + theme(legend.box.margin = margin(0, 0, 0, 10)))


plot_grid(
  NMDS.Tpre+ theme(legend.position = "none"), 
  NMDS.T0.5+ theme(legend.position = "none"), 
  NMDS.T1+ theme(legend.position = "none"), 
  NMDS.T2+ theme(legend.position = "none"), 
  NMDS.T3+ theme(legend.position = "none"),
  extract.legend, 
  rel_widths = c(8,8,8,8,8,3), ncol=6)

ggsave("figures/TBwater.NMDS.pdf", height=5, width=15)


# T1

# T2

# T3


########## 

NMDS.ord.functional<-plot_ordination(
  physeq = PS.fin,                                                   
  ordination = ORD) +                                                
  geom_point(aes(color = Functional.group, shape=Year), size = 3) +    
  stat_ellipse(level=0.9, linetype = 2, aes(color=Functional.group)) +
  #scale_color_brewer(palette = "Dark2") +
  #scale_fill_brewer(palette = "Dark2") + ggtitle("Time1") +
  theme_classic()   

NMDS.ord.functional
dev.copy(pdf, "figures/NMDS.ord.func.group.pdf", height=6, width=9)
dev.off() 

########## 

NMDS.ord.types<-plot_ordination(
  physeq = PS.fin,                                                   
  ordination = ORD) +                                                
  geom_point(aes(color = Sample.type, shape=Year), size = 3) +    
  stat_ellipse(level=0.9, linetype = 2, aes(color=Sample.type)) +
  #scale_color_brewer(palette = "Dark2") +
  #scale_fill_brewer(palette = "Dark2") + ggtitle("Time1") +
  theme_classic()   

NMDS.ord.types
dev.copy(pdf, "figures/NMDS.ord.types.pdf", height=6, width=9)
dev.off() 

########## 

NMDS.ord.lake<-plot_ordination(
  physeq = PS.fin,                                                   
  ordination = ORD) +                                                
  geom_point(aes(color = Lake, shape=Year), size = 3) +    
  #scale_color_brewer(palette = "Dark2") +
  #scale_fill_brewer(palette = "Dark2") + ggtitle("Time1") +
  theme_classic()   

NMDS.ord.lake
dev.copy(pdf, "figures/NMDS.ord.lake.pdf", height=6, width=9)
dev.off() 

### stacked bar plot
pdf(file="figures/stack.bar1.pdf", height=5, width=12)
plot_bar(PS.fin, x="Sample.type", fill="Phylum")
dev.off() 


## other
pdf(file="figures/stack.bar1.pdf", height=5, width=12)
p = plot_bar(PS.fin, "Phylum", fill="Phylum", facet_grid=Year~Organism)
dev.off() 