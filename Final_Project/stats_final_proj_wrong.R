### calculating basic population genetic stats with adegenet
## code adapted from the NESCent population genomics program here: 
## http://popgen.nescent.org/StartSNP.html

library("adegenet")
library("hierfstat")
library("pegas")

## and that your data are space delimited (sep="") 
pines<-read.table("pine_data.txt", sep="", header=TRUE)

# Check the structure of the "state" variable in our pines data set
str(pines$state)
# This tells us that "state" is a "Factor" with "11 levels" - put another way - a 
# categorical variable, with 11 possible categories
# We also could set up this type of structure manually, or add possible categories 
# using a couple options:
names<-sort(unique(pines$state))
names<-c("alabama","arkansas","florida",  "georgia",  "louisiana","mississippi",
         "northcarolina", "oklahoma", "southcarolina" ,"texas", "virginia" )
levels(pines$state)<-c("alabama","arkansas","florida",  "georgia",  "louisiana",
                       "mississippi","northcarolina", "oklahoma", "southcarolina" ,
                       "texas", "virginia" )
levels(pines$state)<-names
#########################################
###### SETTING UP YOUR DATA    ##########
#########################################

#check the size of the file- the dim command gives you the dimensions
dim(pines) 
## should be 550 rows x 3086 columns

### make sure the file looks ok. the head command gives you the first 5 rows 
### of the file and all the columns. 
## Since the file is 3086 rows, we add the [1:10] to the command. This tells R 
## to only show us the first 10 columns
head(pines)[1:10]

## first 4 columns are "metadata"- information about the sample. Next 3084 
## columns are the genetic data
## each row is one individual tree (these are pine trees)
## each column is one locus (one SNP)
## the individual cells show the genotype of the individual at that locus

## we can look at individual rows, columns, and cells:
head(pines)[1:10]
pines[1,1] #this is the data in row one column one- tree id 1066
pines[1,3] #this is the data in row one column three
pines[3,1] #row 3 column 1
pines[,2] #this is EVERYTHING in column 2- we didn't specify a row before the comma, so it shows us all rows
pines[3,] #this is everything in row 3- you'll see its so big we can't even fit it all on our screen
pines[1:3,1:5] ## this is rows 1,2,3 and columns 1,2,3,4,5

## we can also identify columns by their name using the $ operator
pines$tree_id ## this is column 1. same as saying pines[,1]
pines$state ## column 3- same as saying pines[,3]
pines[3] ## this also calls column 3
pines[1:10] ## columns 1:10, all rows

## check the structure of the object pines
str(pines)
## this is a dataframe where each each column is a factor
## each locus has 3 levels- this is good! it means that there are 2 homozygotes 
## and 1 heterozygote (3 different possible values)
## if we found a locus with 4 or 5 levels, we would know there's an error
## NA is missing data

### first we need to get the data formatted correctly for analysis
# To work with the data, we need to convert the R object returned by read.table() to 
# a “genind” object. 
# To achieve this, we will create a dataframe with only genotypes, and keep 
# only a subset of the first 15 SNP loci (to make calculations faster). 
# The small dataset can then be converted to a “genind” object (for package adegenet). 
# The “genind” object can then easily be converted into a “hierfstat” (package hierfstat) object.

## keep only the first 15 loci, cut out the rest of the data. 
## we create a new object, "locus", so we don't overwrite "pines"
## this cuts off the first 4 columns (the metadata), and then all the columns 
## from 20:3086 (the rest of the loci), but keeps all the rows (all individuals)
locus <- pines[, -c(1, 2, 3, 4, 20:3086)]   
head(locus)
dim(locus)

# the column names are confusing (each one is a separate SNP locus), so let's 
# just rename them
head(locus)
# this command tells R to change all the column names to "SNPlocus", followed 
# by a number to the end of locus name, so that each column has a unique identifier
colnames(locus)<-paste0("SNPlocus", "_", seq(1:ncol(locus))) 

## use the "head" command to see the new locus names- same data, just easier 
## column names
head(locus)

## now we're going to create new objects that hold all our metadata- information 
## about sampling location, etc
## we'll store this information in an object separate from the genetic data, and 
## then put it all back together later
## all this data is still accessible in the pines object- that hasn't gone 
## anywhere. It's just hanging out in our workspace waiting for us to use it again

ind <- as.character(pines$tree_id) # labels of the individuals
population <- as.character(pines$state) # labels of the populations

############ Now we're going to format everything for the pop genetics commands

?df2genind
## this command pulls up the help page (lower right screen) for the command 
## that creates the "genind" object, which stores all our genetic information.
## read through the "usage" section- this tells you what each argument within 
## the () does
## X is always your data, and then you enter additional information to tell the 
## df2fgenind command what you want it to do
## in the "Arguments" section, you get more detail about each piece of 
## information you need to tell df2genind
## If you scroll all the way down, there are examples of usage (which are also 
## in the next part of this script)

## Now put your data into the format required for analyzing population genetics- 
## that's what the "df2genind" command does
## in this command, we have a series of arguments telling R that the our genetic 
## data is in the object "locus", 
## we have diploid individuals (ploidy=2),
## our individual ids are in the object "ind", our population information is in 
## the object "population"
## "sep" means that our alleles are not separated by any character- they are 
## AA, AC, AT, etc. 
## If alleles were formatted as "A/A" or "A.A", we would say that sep= "/" or 
## sep=".", respectively

pines_genind <- df2genind(locus, ploidy = 2, ind.names = ind, 
                          pop = population, sep = "") ## create the genind object
# pines_genind@pop<-levels(pines_genind@pop, levels=c("alabama","arkansas","florida",  "georgia",  "louisiana","mississippi",
#                                                        "northcarolina", "oklahoma", "southcarolina" ,"texas", "virginia" ))
## look at the data in your genind object
## we can use the @ command to look at the information in the different 
## "slots" of the pines_genind object (this is the same as the $ operator we 
## used before)
## note that RStudio will tab-complete these for you
pines_genind ## summary of the object
pines_genind@tab
pines_genind@loc.fac
pines_genind@loc.n.all
pines_genind@all.names
pines_genind@ploidy
pines_genind@type
pines_genind@call
pines_genind@pop

nAll(pines_genind) # Number of alleles per locus


#########################################
###### HARDY WEINBERG EQUILIBRIUM  ######
#########################################

### Now that everything is formatted, lets test for HWE:
## we'll call the help command for hw.test
## we see that this will do a traditional chi-square test for deviations from 
## expected values under HWE
?hw.test

hw.test(pines_genind, B=0)

# we also an option to use a different method for testing the significane of 
# deviations from HWE- a Monte Carlo permutation procedure
hw.test(pines_genind, B = 1000)

##################################
###### DIVERSITY STATS  ##########
##################################
## Now that everything is formatted properly we can calculate genetic diversity 
## (observed and expected heterozygosity)
# First, we'll do this with adegenet:
  
## Lets create a new object called "div" that summarizes the information in 
## our pines_genind object
div <- summary(pines_genind)  

div ## this gives us a summary of the information in the div object

names(div)## this tells us the names of each slot of the div object

## for techinal reasons, we select slots in the div object using $ instead of @ 
## (this is because we are working with both S3 and S4 object types- don't stress 
## about this detail)
## use the $ to look at the information contained in each slot
## note that RStudio will tab-complete these for you
div$n
div$n.by.pop
div$loc.n.all
div$pop.n.all
div$NA.perc
div$Hobs
div$Hexp

## Now we can plot the information from the div object using the "plot" command
?plot
# read through the arguments of plot and make sure you understand the code below

## look at proportion of heterozygosity per locus
plot(div$Hobs, xlab="Locus ID", ylab="Observed Heterozygosity", 
     main="Observed heterozygosity per locus")

## the main command here is that we are plotting the Hobs (observed heterozygosity). 
# The rest of the information- xlab, ylab, main- is how we tell R what to label 
# the axes of our plot

## we also want to compare observed and expected heterozygosity
plot(div$Hobs,div$Hexp, xlab="Hobs", ylab="Hexp", 
     main="Expected heterozygosity as a function of observed heterozygosity per locus")

## Now we can test if there a significant difference between observed and expected 
# heterozygosity

# Check the variances are equal
?bartlett.test
bartlett.test(list(div$Hexp, div$Hobs))

?t.test
t.test(div$Hexp, div$Hobs, paired = T, alternative = "greater", var.equal = TRUE)

#########################################
###### BASIC STATS WITH HIERFSTAT  ######
#########################################

########################
## Now, let's do basic statistics with hierfstat  
## we need to convert genind object to a hierfstat object, because sometimes 
## analysis in R is dumb and requires lots of reformatting
## this format is needed to analyze pairwise FST

pines_hier <- genind2hierfstat(pines_genind) # Create hierfstat object

#### 
# Populations are states (e.g. Alabama, Arkansas, etc are the pops) 
# The function basic.stats() provides the observed heterozygosity (Ho), mean 
# gene diversities within population (Hs), Fis, and Fst.
# The function boot.ppfis() provides confidence interval for Fis. 

?basic.stats
## the help page explains how all the different variables are calculated. 
## Some of these should look familiar!

## calculate basic stats on pines dataset
basicstat <- basic.stats(pines_hier, diploid = TRUE, digits = 2)  
basicstat$n.ind.samp
basicstat$pop.freq
basicstat$Ho
basicstat$Hs
basicstat$perloc
basicstat$overall

## just typing in basic stat gives a summary of the information contained in the 
## slots (the perloc and overall options)
basicstat

## we can pull lots of information out of this basic stat object, and compare 
## different measures to each other
## we use the $ operator to select pieces of information
## we can store these pieces of information in their own objects, or plot them 
## against each other

## let's first look at the relationship between sample size and heterozygosity
## we're going to use the "apply" function to summarize the data tables in 
## different slots of our basicstat object

## creates an object where we the average sample size for each locus in each pop 
#this is the mean of each column in the basicstat$n.ind.samp slot
sample.size<-apply(basicstat$n.ind.samp, 2 ,mean) 


ho<-apply(basicstat$Ho, 2 ,mean) ## average Ho for each locus in each pop
#this is the mean of each column in the basicstat$Ho slot

## now we plot those two values agains each other:
plot(sample.size, ho) ## plot Ho against sample size

## Ok- now let's plot Hs (the average heterozygosity of subdivided pops) against Ht (the pooled heterozygosity). remember, we discussed these in lecture! you'll see that these values are columns in the basicstat$perloc slot- we canjust plot those columns against each other

## plot Hs against Ht
plot(basicstat$perloc$Hs, basicstat$perloc$Ht, xlab="Hs in subdivided populations", 
     ylab="Ht in pooled popultions")

plot(basicstat$perloc$Fst)

### different ways of calculating FSt
## weir and cockerham's estimate- corrects for population size
## we wont' really spend time on this in class, but it is the most common way to 
## estimate Fst with real data, so you'll see it in papers
## you may want to use Weir & Cockerham's Fst for your final projects, so here's 
## how to calculate it:

wc(pines_genind)

### pairwise Fst- this gives us the fst between each pair of populations using 
### weir and cockerham's estimator
## this is a measure of genetic distance bewteen each pair of populations
## this calculation might take a minute

fst.dist<-genet.dist(pines_genind, method = "WC84")
fst.dist
max(fst.dist) ## maximum pairwise fst
min(fst.dist) ## minimum pairwise fst
mean(fst.dist) ## mean fst across all pairs

#########################################
###### PAIRWISE FST AND IBD  ##########
#########################################

##############
## some more advanced stuff! let's look for evidence of isolation by distance...
## this code is definitely more sophisticated, so don't worry about trying to 
## understand all the pieces...for now (unless you want to)
## first, we load the right packages
## if you get a message asking to restart R, click "no"

install.packages(c("geonames","data.table","purrr","geosphere"))
library(geonames)
library(data.table)
library(purrr)
library(geosphere)

## create an object that contains the name of each state where sampling was done
## this pulls out each unique state from our pines dataset and puts them in a vector
states<-sort(unique(pines$state)) 
levels(states)[9]<-"south carolina" # put a space in "south carolina"
states

## use online database to download geographic coordinates for each state
## this is basically telling R to take the list of our state names, search an online 
## database for the gps coordinates for those states, 
## and then download the coordinates into R. Pretty cool!

options(geonamesUsername = "escordato")
get_coords <- function(name, country) {
  res <- GNsearch(name = name, state=name, country = "US")  
  res<-res[1,]
  res$lng<-as.numeric(res$lng)
  res$lat<-as.numeric(res$lat)
  out <- data.frame(name = res$adminName1, lat = res$lat, lon = res$lng)
  return(out)  
}

## this condenses the geographic coordinates for each state into a nice dataframe
state.coords <- states %>% 
  map(get_coords, country = "US") %>% 
  rbindlist()
state.coords<-as.data.frame(state.coords)
state.coords ## ta da! a list of geographic coordinates for each state

## now, we convert the coordinates from each state into a pairwise geographic distance matrix
## this will give us the distance in km between each state
## we're obviously oversimplifying a bit here, but it gets us close

longlats<-state.coords[,c(3,2)] ## pull out just longitude and latitude
longlats
state.dists<-distm(longlats) ## this converts the coordinates into a distance matrix
rownames(state.dists)<-state.coords$name # name the rows
colnames(state.dists)<-state.coords$name # name the columns
state.dists<-state.dists/1000 ## make kilometers
state.dists

## ok! now we have geographic and genetic distance between sites. lets plot them...
## first, convert to correct format for plotting
## we have to re-order our fst distance matrix to be in the same order as the states- 
## we'll alphabetize them both
ordering <- sort(attr(fst.dist, "Labels"))
fst.mat <- as.matrix(fst.dist)[ordering, ordering]
fst.dist<-as.dist(fst.mat, diag=TRUE)
fst.dist

state.dists<-as.dist(state.dists, diag=TRUE)
state.dists

## now plot the genetic distance (fst.dist) agains the geographic distance (state.dists)
plot(state.dists, fst.dist, pch=16, cex=1.5, col="blue", xlab="geographic distance (km)", 
     ylab="genetic distance (pairwise WC Fst)")
