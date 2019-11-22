library("adegenet")
library("hierfstat")
library("pegas")
data<-read.table("Dataset3.genind.txt", sep="", header=TRUE)
head(data[1:10])
loci<-data[,c(-1,-2)]
head(loci[1:10])
str(loci)

# labels of the individuals
index <- as.character(data$sample)
# labels of the populations
pops <- as.character(data$pop)

data_genind <- df2genind(loci, ploidy = 2, ind.names = index, 
                          pop = pops, sep = "/") 

data_genind
#View(data_genind)
data_genind@tab
nAll(data_genind) 


hw<-hw.test(data_genind, B=0)


pvals<-hw[,3]
num_loci<-subset(pvals,pvals<0.05)
#View(num_loci)
print("Number of loci outside of HWE:") 
length(num_loci)

# MCMC Permutation
HW<-hw.test(data_genind, B = 1000)
PVALS<-HW[,4]
NUM_loci<-subset(PVALS,PVALS<0.05)
print("Number of loci outside of HWE:") 
length(NUM_loci)

##################################
###### DIVERSITY STATS  ##########
##################################
## Now that everything is formatted properly we can calculate genetic diversity (observed and expected heterozygosity). First, we'll do this with adegenet:

## Lets create a new object called "div" that summarizes the information in our pines_genind object
DIV <- summary(data_genind)
DIV ## this gives us a summary of the information in the div object
names(DIV)## this tells us the names of each slot of the div object
DIV$n
DIV$n.by.pop
DIV$loc.n.all
DIV$pop.n.all
DIV$NA.perc
DIV$Hobs
DIV$Hexp

## look at proportion of heterozygosity per locus
plot(DIV$Hobs, xlab="Locus ID", ylab="Observed Heterozygosity", 
     main="Observed heterozygosity per locus")

## we also want to compare observed and expected heterozygosity
plot(DIV$Hobs,DIV$Hexp, xlab="Hobs", ylab="Hexp", 
     main="Expected vs.observed heterozygosity per locus")
lm<-lm(DIV$Hexp~DIV$Hobs)
summary(lm)
abline(lm, col="red", lwd=3)

## Now we can test if there a significant difference between observed and expected heterozygosity
# Check the variances are equal
bartlett.test(list(DIV$Hexp, DIV$Hobs))
### VARIANCE ISN'T EQUAL: p-value < 2.2e-16: Probably need to do some data transformation

t.test(DIV$Hexp, DIV$Hobs, paired = T, alternative = "greater", var.equal = TRUE)

#########################################
###### BASIC STATS WITH HIERFSTAT  ######
#########################################
## Now, let's do basic statistics with hierfstat. we need to convert genind object to a hierfstat object, because sometimes analysis in R is dumb and requires lots of reformatting. this format is needed to analyze pairwise FST

data_hier <- genind2hierfstat(data_genind) # Create hierfstat object. Populations are sites

# The function basic.stats() provides the observed heterozygosity (Ho), mean gene diversities within population (Hs), Fis, and Fst.

# The function boot.ppfis() provides confidence interval for Fis. 

## calculate basic stats on dataset
basicstat_final <- basic.stats(data_hier, diploid = TRUE, digits = 2)  
## just typing in basic stat gives a summary of the information contained in the slots (the perloc and overall options)
basicstat_final

## let's first look at the relationship between sample size and heterozygosity. we're going to use the "apply" function to summarize the data tables in different slots of our basicstat_final object. creates an object where we the average sample size for each locus in each popthis is the mean of each column in the basicstat_final$n.ind.samp slot
sample.size<-apply(basicstat_final$n.ind.samp, 2 ,mean) 

ho<-apply(basicstat_final$Ho, 2 ,mean) ## average Ho for each locus in each pop. this is the mean of each column in the basicstat_final$Ho slot

## now we plot those two values agains each other:
plot(sample.size, ho) ## plot Ho against sample size

## Ok- now let's plot Hs (the average heterozygosity of subdivided pops) against Ht (the pooled heterozygosity). remember, we discussed these in lecture! you'll see that these values are columns in the basicstat_final$perloc slot- we canjust plot those columns against each other

## plot Hs against Ht
plot(basicstat_final$perloc$Hs, basicstat_final$perloc$Ht, xlab="Hs in subdivided populations", 
     ylab="Ht in pooled popultions")
plot(basicstat_final$perloc$Fst)

# Weir & Cockerham's Fst
wc(data_genind)
# $FST
# [1] 0.2154904
# $FIS
# [1] 0.04863081

### pairwise Fst- this gives us the fst between each pair of populations using weir and cockerham's estimator. This is a measure of genetic distance bewteen each pair of populations. this calculation might take a minute
fst.dist<-genet.dist(data_genind, method = "WC84")
fst.dist
max(fst.dist) ## maximum pairwise fst 0.3015741
min(fst.dist) ## minimum pairwise fst 0.008603233
mean(fst.dist) ## mean fst across all pairs 0.187665

#########################################
###### PAIRWISE FST AND IBD  ##########
#########################################
## some more advanced stuff! let's look for evidence of isolation by distance.
## if you get a message asking to restart R, click "no"

#install.packages(c("geonames","data.table","purrr","geosphere"))
library(geonames)
library(data.table)
library(purrr)
library(geosphere)

# Read in site coordinates:
site.coords<-read.csv("site.coords.csv",header=TRUE)
head(site.coords)
str(site.coords)

## Now, we convert the coordinates from each site into a pairwise geographic distance matrix.T his will give us the distance in km between each site. We're obviously oversimplifying a bit here, but it gets us close.
longlats<-site.coords[,c(3,2)] # Pull out just longitude and latitude
longlats
site.dists<-distm(longlats) # This converts the coordinates into a distance matrix
rownames(site.dists)<-site.coords$name # name the rows
colnames(site.dists)<-site.coords$name # name the columns
site.dists<-site.dists/1000 ## make kilometers
site.dists

# Ok! now we have geographic and genetic distance between sites. Lets plot them. First, convert to correct format for plotting. We have to reorder our Fst distance matrix to be in the same order as the sites - we'll alphabetize them both.
ordering <- sort(attr(fst.dist, "Labels"))
fst.mat <- as.matrix(fst.dist)[ordering, ordering]
fst.dist<-as.dist(fst.mat, diag=TRUE)
fst.dist
site.dists<-as.dist(site.dists, diag=TRUE)
site.dists

# Now plot the genetic distance (fst.dist) agains the geographic distance (site.dists)
plot(site.dists, fst.dist, pch=16, cex=1.5, col="blue", xlab="Geographic Distance (km)", 
     ylab="Genetic Distance (pairwise WC Fst)")

