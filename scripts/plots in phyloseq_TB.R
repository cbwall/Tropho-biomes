## investigate the final 16S library


remotes::install_github("microbiome/microbiome")
library("microbiome")
library("vegan")
library("tidyr")
library("ggplot2")
library("phyloseq")


# ps.prune is the final phyloseq object from pipeline, can load it in here...
### or load in the raw data and re-assemble

ps.prune<-readRDS("output/TSCC_output/ps.prune.RDS")


######### an error in IDing one of the treatments, read in the correct metadata below and add this to the phyloseq object

# read in metadata and re-format so it has the str necessary
M<-read.csv("output/TSCC_output/sam_data.csv")

# set row names
row.names(M)<-M$X 

# remove junk column
M<-M[-1]

make.fac<-c("Location", "Experiment", "Time", "Tank", 
            "Inoc.plank.nutr.source", "Inoc.treatment", "Zoop.tube_Inoc.tube.ID", "Inoc.Bottle.color",
            "Inoc.Bottle.number", "rep",  "Sample.type", "Organism", 
            "Functional.group", "sample_control", "Miseq.ANL")

run.metaD[make.fac] <- lapply(run.metaD[make.fac], factor) # make all these factors



sample_data(ps.prune)<-M

## skip this
#### assemble the phyloseq object from raw files
PS.fin<-
  read_phyloseq(
    otu.file = "output/TSCC_output/otu_table.csv",
    taxonomy.file = "output/TSCC_output/tax_table.csv",
    metadata.file = "output/TSCC_output/sam_data.csv",
    sep = ","
  )

#replace sample data with reformatted metadata
sample_data(PS.fin)<-M



################ Now on to plotting

PS.fin<-ps.prune


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


###### some richness plots
#richness by Lake
richness.plot_inoc.trt<-plot_richness(PS.fin, x="Inoc.treatment", measures=c("Observed", "Shannon"))  + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
richness.plot_inoc.trt
dev.copy(pdf, "figures/rrichness.plot_inoc.trt.pdf", height=4, width=10)
dev.off() 


#richness by organism (or sample type i.e, water)
richness.plot_organism<-plot_richness(PS.fin, x="Inoc.plank.nutr.source", measures=c("Observed", "Shannon")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
richness.plot_organism
dev.copy(pdf, "figures/richness.plot.organism.pdf", height=4, width=12)
dev.off() 


#richness by Time.point
richness.plot_Time<-plot_richness(PS.fin, x="Time", measures=c("Observed", "Shannon")) + theme_bw()
richness.plot_Time
dev.copy(pdf, "figures/richness.plot_Time.pdf", height=4, width=10)
dev.off() 


#richness by Experiment
richness.plot_Exp<-plot_richness(PS.fin, x="Experiment", measures=c("Observed", "Shannon")) + theme_bw()
richness.plot_Exp
dev.copy(pdf, "figures/richness.plot_Exp.pdf", height=4, width=10)
dev.off() 


#richness by Functional Group
richness.plot_Func<-plot_richness(PS.fin, x="Functional.group", measures=c("Observed", "Shannon")) + theme_bw()
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