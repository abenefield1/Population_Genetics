###################################################
######## OUTFLANK LAB ########
###################################################

# in this lab we will identify loci that are highly differentiated among populations
# these loci are candidates for being under selection

###################################################
######## INSTALLING AND LOADING PACKAGES ########
###################################################

### IF YOU ARE ON YOUR OWN LAPTOP:

# load outFLANK
##  install the packages "devtools","ggplot2", "qqman" and "gtools" from the packages tab if you do not already have them
## use the library command to load the packages you need
library("ggplot2")
library("qqman")
library("gtools")
library("devtools")
library("qvalue")
library("OutFLANK")

###################################################
######## LOADING DATA INTO R ########
###################################################

## we are going to use the same datasets as we used with TESS last week. We just need to format them a little differently

## Outflank datasets:
# SNP are in columns: each row contains one value at each SNP locus 
# Entries in cells correspond to the number of reference alleles: 
# 0 means zero copies of the reference allele, 1 means one copy of the reference allele 
# (i.e. heterozygote), and 2 means two copies of the reference allele.
# reference allele is arbitrary- just remember that 0 and 2 are homozygotes and 1 is heterozygote

##########
## READ IN DATASETS
######
## you can  change your input datset here for datasets 1, 2, and 3
dataset.OF<-read.table("dataset3_OF.txt", header=TRUE)

### once you load your input dataset, you can run the rest of the code without changing anything


###################################################
######## FORMATTING AND CHECKING OUR DATA ########
###################################################

### first, we need to do a little reformatting of the dataset:
head(dataset.OF)[1:10]

dataset.OF[is.na(dataset.OF)]<-9 ## convert NAs to 9
## the is.na command finds all the NAs (missing data) in your dataset
## the <-9 changes the NAs to 9's, becaue that's what outFlank wants
## this takes a minute

head(dataset.OF)[1:10] ## check your data- note that the NAs have turned in to 9s
dim(dataset.OF)

## make a vector with your population information
## we will use this to tell OutFLANK how many populations there are
pop.names<-dataset.OF$pop
unique(pop.names) ## tells you the name of each population
nlevels(pop.names) ## tells you how many factor levels in the vector. There is one factor level for each population, so this tells you how many pops

## create a dataframe with just the genotypes
genotypes<-droplevels(dataset.OF[,-1]) # the -1 removes the population column
head(genotypes[1:10])
## the droplevels command here is a weird R thing- don't worry about it for now

## COMPLETE THE COMMAND
# If samples are rows, and loci are columns, how would we check for the 
# number of each of these?
length(genotypes$locus1)
dim(genotypes)


## COMPLETE THE COMMAND
## create a vector that has the names of your loci- these are the column names of the genotype dataframe
# try checking out the "colnames" command - very handy!
locus.names<-colnames(genotypes)


###################################################
######## CALCULATE FST FOR EACH LOCUS ########
###################################################

# first we'll create a matrix of all Fsts for each locus
# you can use ?MakeDiploidFSTmat to see what this is doing

Fsts<-MakeDiploidFSTMat(genotypes, locus.names, pop.names ) #define your data, the locus names, and the populations

## how would we look at just the first few rows from the above command
head(Fsts)

### look at the distribution of Fst values across your loci by making a histogram with the hist function
hist(Fsts$FSTNoCorr, xlab="Fst")

### Check out the relationship between Fst and heterozygosity
plot(Fsts$He, Fsts$FST, xlab= "heterozygosity", ylab="Fst")
abline(lm(Fsts$FST~Fsts$FST), col="red")
summary(lm(Fsts$FST~Fsts$FST))

## lets look at the maximum Fst
max(Fsts$FSTNoCorr,na.rm=T)

## look at the mean Fst of the FSTNoCorr column
mean(Fsts$FSTNoCorr,na.rm=T)

###################################################
######## RUNNING OUTFLANK TO IDENTIFY OUTLIERS ########
###################################################

### Ok, now we're going to calculate the number of outlier loci (those that might be under selection)
##you can use ?OutFLANK for more information and to read about the specifications in the command

Fst.outliers<-OutFLANK(Fsts,LeftTrimFraction=0.05,
                        RightTrimFraction=0.05, Hmin=0.1, NumberOfSamples = nlevels(pop.names),
                        qthreshold=0.05)


## plot the results of the outlier scan
OutFLANKResultsPlotter(Fst.outliers, withOutliers = TRUE, NoCorr = TRUE, Hmin = 0.1, 
                       binwidth = 0.005, Zoom = FALSE, RightZoomFraction = 0.05, 
                       titletext = NULL)


### same plot as above, but zoomed in on the tail of the distribution
# you can see we changed Zoom=TRUE- everything else is the same
OutFLANKResultsPlotter(Fst.outliers, withOutliers = TRUE, NoCorr = TRUE, Hmin = 0.1, 
                       binwidth = 0.005, Zoom = TRUE, RightZoomFraction = 0.05, 
                       titletext = NULL)


## now we can make a dataframe containing all the outlier loci
## we do this by pulling out all the results in our outlier list where Oultier= TRUE
Fst.outlier.list<-Fst.outliers$results[Fst.outliers$results$OutlierFlag==TRUE,]
dim(Fst.outlier.list) ## how many outliers?
print("The number of outliers is")
dim(Fst.outlier.list)
head(Fst.outlier.list)

mean(Fst.outlier.list$FST, na.rm=TRUE)

###################################################
######## PLOTTING OUTLIERS ########
###################################################

### Let's make a manhattan plot and highlight the outlier loci

## calculate pvalues and qvalues for all loci
ManPlots <- pOutlierFinderChiSqNoCorr(Fsts, Fstbar = Fst.outliers$FSTNoCorrbar, 
                                      dfInferred = Fst.outliers$dfInferred, qthreshold = 0.05, Hmin=0.1)

## here we are adding a column to our output to make plotting look nicer
ManPlots$map.name<-seq(1:nrow(ManPlots))

## make the manhattan plot
plot(ManPlots$map.name[ManPlots$He>0.1], ManPlots$FST[ManPlots$He>0.1],
     xlab="locus", ylab="FST", col=rgb(0,0,0,0.2), ylim=c(0,1),main="Dataset 3")

## create a vector of outliers
outliers<- ManPlots$OutlierFlag==TRUE

## color the outlier points on the plot
points(ManPlots$map.name[outliers], ManPlots$FST[outliers], col="magenta", pch=20)   












