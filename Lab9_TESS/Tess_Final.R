####### How to run structure-type analyses in R
# Download packages
#install.packages(c("fields","RColorBrewer","mapplots", "dplyr","ggplot2", "devtools","rworldmap", "maps","raster")) ## install packages
#install.packages(c("fields","mapplots", "devtools","rworldmap","raster")) ## install packages

# if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("LEA")
#devtools::install_github("bcm-uga/TESS3_encho_sen")

#devtools::install("/Users/Amy 1/Desktop/Population_Genetics/Lab/R/Lab9/TESS3_encho_sen-master")
library("devtools")

#BiocManager::install("qvalue")
#install_github("whitlock/OutFLANK")
library(OutFLANK)

### load packages
library(LEA)
library(tess3r)
library(mapplots)
library(fields)
library(dplyr)
library(ggplot2)
library(raster)
library(rworldmap)
library(maps)
source("http://membres-timc.imag.fr/Olivier.Francois/Conversion.R") ## supplemental functions
source("http://membres-timc.imag.fr/Olivier.Francois/POPSutilities.R") ## supplemental functions

#### using TESS:
## https://bcm-uga.github.io/TESS3_encho_sen/articles/main-vignette.html
## this is the TESS website, which has a nice guide for using the program

######################################################################
## set your working directory to where the files are located

#read in the data
dataset1<-read.table("dataset1.geno")
dim(dataset1)
head(dataset1)[1:10] # 226 individuals: 15000 SNPs
coords1<-read.table("dataset1.coord", header=FALSE, sep="")
head(coords1)
coords1<-as.matrix(coords1) ## convert to a matrix for the program to work properly

## this is the main clustering command
#dataset1.TESS<- tess3(X = dataset1,
#                      coord = coords1,
#                      ploidy = 2,
#                      K = 1:6,
#                      rep = 1,
#                      max.iteration=1000)
# Write to file:
#saveRDS(dataset1.TESS, file = "dataset1.TESS.rds")
dataset1.TESS<-readRDS(file="dataset1.TESS.rds")

# we tell it the dataset (X) = dataset1
# the coord = the coords1 object
# ploidy is diploid (2)
# fit K from 1-6
# do 1 repetition per K
# iterate the program 1000 times

# now we plot the cross validation to help choose K
#The best choice for the K value is when the cross-validation curve exhibits a plateau or starts increasing.
plot(dataset1.TESS, pch = 19, col = "blue",
     xlab = "Number of ancestral populations",
     ylab = "Cross-validation score")

# OK, let's make structure-style barplots now. The q matrix is the object TESS created that has the assignment proportions for each individual. The number of columns in this object corresponds to the number of clusters (k)- so K=3 has 3 columns. Each row is an individual. The values for each column for each individual (row) are the proportions of its genotype assigned to each cluster. So to make the barplot, we're going to plot the values for each individual as one bar. You can envision this as the rows of the q matrix being tipped on their sides, so each row is turned into one bar on the bar plot. The amount of each color in each bar corresponds to the proportions of ancestry- so an individual with 0.5 in each cluster will have a bar that is half one color, half the other color. You can change the value of K here to make plots for different K values:
q.matrix.dataset1 <- qmatrix(dataset1.TESS, K = 2) 
head(q.matrix.dataset1)

#####
# make the STRUCTURE-like barplot for the Q-matrix 
barplot(q.matrix.dataset1, border=NA, space = 0,
        xlab = "Individuals", ylab = "Ancestry proportions", 
        main = "Ancestry matrix") ->bp

## this command makes a bar plot from the q.matrix
# border = "gray50" is the color of the borders of the bars. You can try changing this, or set border = NA
# space = 0 means no spaces between bars-  try deleting this and re-doing your plot and see what happens

########## Now, let's analyze our second dataset
#read in the data
dataset2<-read.table("dataset2.geno")
dim(dataset2)
head(dataset2)[1:10]
coords2<-read.table("dataset2.coord", header=FALSE, sep="")
head(coords2)
coords2<-as.matrix(coords2)

## do the clustering analysis, using the same parameters as before
# dataset2.TESS<- tess3(X = dataset2,
#                       coord = coords2,
#                       ploidy = 2,
#                       K = 1:6,
#                       rep = 1,
#                       max.iteration = 10000)

plot(dataset2.TESS, pch = 19, col = "blue",
     xlab = "Number of ancestral populations",
     ylab = "Cross-validation score")

## you can change the value of K here to make plots for different K values
q.matrix.dataset2 <- qmatrix(dataset2.TESS, K = 3) 
head(q.matrix.dataset2)

# make the STRUCTURE-like barplot for the Q-matrix 
barplot(q.matrix.dataset2, space = 0, border=NA,
        xlab = "Individuals", ylab = "Ancestry proportions", 
        main = "Ancestry matrix") -> bp


## NOW! Do this all on your own for dataset 3. type the code in here and turn this script in with your lab assignment
# yes, you may copy and paste the code from above and just modify it appropriately for dataset 3. 
## You might want to consider setting border = NA in your barplot
dataset3<-read.table("dataset3.geno")
dim(dataset3)
head(dataset3)[1:10]
coords3<-read.table("dataset3.coord", header=FALSE, sep="")
head(coords3)
coords3<-as.matrix(coords3)

## do the clustering analysis, using the same parameters as before
dataset3.TESS<- tess3(X = dataset3,
                      coord = coords3,
                      ploidy = 2,
                      K = 1:6,
                      rep = 1,
                      max.iteration = 10000)

plot(dataset3.TESS, pch = 19, col = "blue",
     xlab = "Number of ancestral populations",
     ylab = "Cross-validation score",
     main="Dataset 3: 826 Individuals, 300 SNPs")

## you can change the value of K here to make plots for different K values
q.matrix.dataset3 <- qmatrix(dataset3.TESS, K = 4) 
head(q.matrix.dataset3)

# make the STRUCTURE-like barplot for the Q-matrix 
barplot(q.matrix.dataset3, space = 0, border=NA,
        xlab = "Individuals", ylab = "Ancestry proportions", 
        main = "Ancestry matrix") -> bp

### A cool thing about TESS is that it is a spatially explicit model
## that is, we can incorporate information about geographic sampling location into our models
## here, let's plot how our measures of ancestry (assignment to different genetic clusters) vary across the landscape
## we'll use the coords file to tell TESS where each sample came from (its geographic coordinates)

## here, you can change which q.matrix and coord objects are entered
## so, if you want to plot the q.matrix from dataset3, you can change this to q.matrix.map<-q.matrix3

## just make sure the q.matrix and coords match- so q.matrix from dataset3 and coords3 etc
## you can go back and change the value of K in your q-matrix, and then see how your map changes

q.matrix.map<-q.matrix.dataset3 ## change this value for different datasets
coords<-coords3 ## change this value for different coordinates

## load a map of the world
#install.packages("rworldxtra")
library(rworldxtra)

map.polygon <- getMap(resolution = "li")
quartz()
plot(map.polygon) 

## set colors
my.colors <- c("tomato", "orange", "lightblue", "wheat","olivedrab","mediumorchid2","darkgray")
my.palette <- CreatePalette(my.colors, 9)
# here I create a different combination of colors for each dataset at the K=2 model
palette.1<-my.palette[c(3,4)]
palette.2<-my.palette[c(1,5,6,7)]
palette.3<-my.palette[c(1,2)]

### make plots
##
## if you want to be fancy, change the col.palette to the the corresponding number of your dataset
## so if you're plotting datset 3, change to col.palette=palette.3
## these only work with K=2
## if you are plotting a q-matrix with a larger k than 2, set col.palette = my.palette

# this function connects the your genetic clusters and coordinates with the map you downloaded
q.matrix.map<-q.matrix.dataset3 ## change this value for different datasets
coords<-coords3
# it may take a few seconds to run
ancestry.coeffs <- ggtess3Q(q.matrix.map, coords, map.polygon = map.polygon, col.palette = palette.2) 
# now we plot- everything after "ancestry.coeffs" is about aesthetics of the map
ancestry.coeffs +
  geom_path(data = map.polygon, aes(x = long, y = lat, group = group)) +
  xlim(30, 140) + 
  ylim(20, 75) + 
  coord_equal() + 
  geom_point(data = as.data.frame(coords), aes(x = V1, y = V2), size = 1) + 
  xlab("Longitude") +
  ylab("Latitude") + 
  theme_bw()+
  ggtitle("Dataset 3, k=4")

### Dataset2
map.polygon <- getMap(resolution = "li")
quartz()
plot(map.polygon) 
q.matrix.map<-q.matrix.dataset2 ## change this value for different datasets
coords<-coords2
ancestry.coeffs <- ggtess3Q(q.matrix.map, coords, map.polygon = map.polygon, col.palette = palette.2) 
ancestry.coeffs +
  geom_path(data = map.polygon, aes(x = long, y = lat, group = group)) +
  xlim(30, 140) + 
  ylim(20, 75) + 
  coord_equal() + 
  geom_point(data = as.data.frame(coords), aes(x = V1, y = V2), size = 3) + 
  xlab("Longitude") +
  ylab("Latitude") + 
  theme_bw()+
  ggtitle("Dataset2, k=3")



# Dataset1
q.matrix.map<-q.matrix.dataset1 ## change this value for different datasets
coords<-coords1 
map.polygon <- getMap(resolution = "li")
quartz()
plot(map.polygon) 
ancestry.coeffs <- ggtess3Q(q.matrix.map, coords, map.polygon = map.polygon, col.palette = palette.1) 
ancestry.coeffs +
  geom_path(data = map.polygon, aes(x = long, y = lat, group = group)) +
  xlim(30, 140) + 
  ylim(20, 75) + 
  coord_equal() + 
  geom_point(data = as.data.frame(coords), aes(x = V1, y = V2), size = 3) + 
  xlab("Longitude") +
  ylab("Latitude") + 
  theme_bw()+
  ggtitle("Dataset1, k=2")
