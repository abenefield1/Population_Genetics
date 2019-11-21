### today we're going to be looking at population structure in our pines data
##  the first part of the lab will involve the same setup as last week- 
##  we will read in our data and create a genind object

### calculating basic population genetic stats with adegenet
## code adapted from the NESCent population genomics program here: 
## http://popgen.nescent.org/StartSNP.html

## install the required packages
## you can also do this "point and click" style on the lower-right 
## panel under "packages" --> "install.packages
## if you already installed the packages last week, you don't need to 
## re-install them. You just need to load them with the library commands
install.packages(c("adegenet","hierfstat","pegas", "dplyr"))

library("adegenet")
library("hierfstat")
library("pegas")
library("dplyr")

## set your working directory here. This is the folder on your computer 
## you can also do this in the lower-right window of RStudio by going to files, 
## navigating to your folder, and clicking more--> set as working directory
## you will need to edit this line to your correct directory

pines<-read.table("pine_data.txt", sep="", header=TRUE)

## keep only the first 15 loci, cut out the rest of the data. 
## we create a new object, "locus", so we don't overwrite "pines"
## this cuts off the first 4 columns (the metadata), and then all the columns 
## from 20:3086 (the rest of the loci), but keeps all the rows (all individuals)
#locus <- pines[, 5:19]   
locus <- pines[, c(-1,-2,-3,-4)]   
dim(locus)

# the column names are confusing (each one is a separate SNP locus), 
# so let's just rename them
colnames(locus)<-paste0("SNPlocus", "_", seq(1:ncol(locus))) 
head(locus)[1:10]

## now we're going to create new objects that hold all our metadata- information 
## about sampling location, etc
## we'll store this information in object separate from the genetic data, 
## and then put it all back together later
## all this data is still accessible in the pines object- that hasn't gone anywhere. 
## It's just hanging out in our workspace waiting for us to use it again

ind <- as.character(pines$tree_id) # labels of the individuals
population <- as.character(pines$state) # labels of the populations

############ Now we format everything for the pop genetics commands in adegenet and hierfstat
## you can use the help command ?df2genind to remind yourself what 
## this function is doing
## create the genind object
pines_genind <- df2genind(locus, ploidy = 2, ind.names = ind, 
                          pop = population, sep = "") 
pines_genind

pines_hier <- genind2hierfstat(pines_genind) # Create hierfstat object

######################################################
######### Let's look at some simple patterns of clustering using indpca, 
######### find.clusters, and dapc
###############################################################

# The function indpca() does a PCA on the centered matrix of individuals’ 
## allele frequencies.
## we can look at genetic relationships based on allele frquency differences 
## among INDIVIDUALS using this method
## read about the analysis here:
?indpca

## do the PCA
pca_pines<- indpca(pines_hier) 
## plot results 
plot(pca_pines, cex = 0.7)
## Okay that sucks. Lets make it better
# We'll pull out population names
indlabels<-pines_hier$pop
# Then the first two pc axes from our pca_pines object, and plot them by color
plot(pca_pines$ipca$l1$RS2~pca_pines$ipca$l1$RS1,pch=19,col=indlabels, xlab="PC1", ylab="PC2",main="Allele Freq. Differences Among Individuals")
legend("bottomleft",legend=unique(indlabels),col=1:length(indlabels),pch=19, cex=0.8)

#### Now lets use unsupervised clustering to try to find the number of genetic clusters that best describes our data
## since we don't know how many genetic clusters there are, we will let the clustering algorithm find the best number for us
## We use the function find.clusters() with the maximum number of PCA axes, which is the same number as our loci. 
## The idea here is to first apply a PCA to our data, which reduces the number of dimensions
## we then use a discriminant function to look at how variation in our data is spread across different numbers of genetic clusters (k)
## Different statistics are used to then choose the "best" (most likely) number of genetic clusters in the data
## the statistic used is chosen with the "criterion" argument- here we will use the default
## our goal is to choose a number of clusters that maximizes the variation among groups while minimizing variation within groups
## we will also get an output of how confident we are in the number of clusters in our data- that is, a membership probability for each individual to each clusters
## In this example, we used choose.n.clust = FALSE. However, if you use the option TRUE and then you will be able to choose the number of clusters.

## check out the help page here- you can alter the arguments to this 
## function in lots of different ways. feel free to explore!
?find.clusters

set.seed(20160308) # Setting a seed for a consistent result
## cluster selection procedure
clusters <- find.clusters(pines_genind, n.pca=15, max.n.clust = 10, 
                          choose.n.clust = FALSE) 
# explore the information contained in the clusters object
names(clusters) 
clusters$Kstat
clusters$stat
clusters$grp
clusters$size

### now we're going to plot the number of genetic clusters
dapc_results <- dapc(pines_genind, clusters$grp, n.pca = 10, n.da = 5) 
## the n.pca and n.da arguments tell the function how many PCs and 
## linear discriminants to use. you can change these numbers and see 
## how your results change

dapc_results
summary(dapc_results) ## summarizes your dapc
scatter(dapc_results, col=c("darkblue","purple","green","orange","red","blue")) # plot of the dapc results
## the ellipses are summaries of the cloud of points
## the lines connect each point to the center of its cluster
## the col argument allows you to assign a color to each cluster 
## (this works for all plots, not just dapc)


### let's now plot the same data, but in a different way- 
### make the points bigger and transparent, remove the ellipses, 
### add a legend
### visualizing data in a such a way that it's easy to pick 
### out the patterns is an important skill!
### feel free to play around with the graphical parameters 
### and see how they affect your plot

colors<-c("darkblue","purple","green","orange","red","blue")## define the colors you use in the plot
scatter(dapc_results, scree.da=FALSE,  pch=20, cell=0, cstar=0, col=colors, solid=.4,
cex=2,clab=0, leg=TRUE, txt.leg=paste("Cluster",1:5)) ## plot data

### now we can add a minimum spanning tree between these clusters
## here the lines show the distance between clusters in PC space- 
## longer lines indicate clusers that are more differentiated
## these lines do NOT translate to geographic distances- 
## they are distances between clusters in PC space

## add the minimum spanning tree with the argument mstree=TRUE
scatter(dapc_results, cstar=0, mstree= TRUE) 

## make the plot look nicer:
scatter(dapc_results, ratio.pca=0.3, pch=20, cell=0,
        cstar=0, col=colors, solid=.4, cex=2, clab=0,
        mstree=TRUE, scree.da=FALSE,
        leg=TRUE, lwd=3, txt.leg=paste("Cluster",1:3))
points(dapc_results$grp.coord[,1], dapc_results$grp.coord[,2], pch=4,
       cex=3, lwd=8, col="black") ## add "x' to group center
points(dapc_results$grp.coord[,1], dapc_results$grp.coord[,2], pch=4,
       cex=3, lwd=2, col=colors)


################################
#### Now let's try looking at our data like tree using the nj commands
####################################

## we'll just make a very simple neighbor-joining tree
## so that the tree is readable, we'll randomly sample 50 individual 
## trees from the dataset, but use all the loci:
## the "sample_n" command is telling R to randomly 
## pull 50 rows from the pines dataset and put them in a new 
## object called nj.pines
nj.pines<-sample_n(pines, 50)

## now we'll format our "nj.pines" dataset into a genind object exactly like we did above- 
## we can just copy-paste the same code but change the object names so we don't overwrite our previous analysis
nj.locus <- nj.pines[, -c(1, 2, 3, 4)]   
colnames(nj.locus)<-paste0("SNPlocus", "_", seq(1:ncol(locus))) 
head(nj.locus)[1:10]
nj.ind <- as.character(nj.pines$tree_id) # labels of the individuals
nj.population <- as.character(nj.pines$state) # labels of the populations

nj.pines_genind <- df2genind(nj.locus, ploidy = 2, ind.names = nj.ind, pop = nj.population, sep = "") ## create the genind object
## if you get a weird error message here about duplicated labels, just ignore it

### now that we have our 50 individuals formatted, let's make a tree
?nj
?plot.phylo

## make the tree
quartz()
nj.tree<-nj(dist(as.matrix(nj.pines_genind))) ## create the tree based on a genetic distance matrix
nj.tree$tip.label<-clusters$grp ## change tip labels to population of orign
plot(nj.tree, typ = "fan", use.edge.length = FALSE, cex=0.7) ## plot the tree

tiplabels(pch=20, col=nj.pines_genind@pop, cex=2) ## color each tip by population, using the colors we just created
title("Neighbor-joining tree of pine data")

### you can change the "typ" of tree to view it in different formats
#try getting rid of the "typ" argument and looking at the default: 
plot(nj.tree, use.edge.length = FALSE, cex=0.6) ## plot the tree
tiplabels(pch=20, col=nj.pines_genind@pop, cex=1) ## color each tip by population, using the colors we just created
title("Neighbor-joining tree of pine data")


