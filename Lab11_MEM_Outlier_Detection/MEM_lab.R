### Moran spectral outlier detection to identify outlier loci
### Moran spectral randomizaiton to associate outlier loci with environmental variables

## reminder- if you are on your own laptop, you don't need to install packages if you have already installed them- just call them with the library command
## if you are on the lab desktop you need to install packages because they are erased when you logout
# either run the code below or click on the "packages" tab and install each one manually
install.packages(c("adespatial"), dependencies=TRUE)

# Load packages
# -------------
library(spdep)
library(adespatial)
library(raster)
library(viridis)


## set your working directory to where your data and habitat files are located
#setwd("")

## COMPLETE THE COMMAND ##
## load the file "MEM_data.csv" using the read.csv function. Save this to the variable 'data'
data   <- read.csv("MEM_data.csv")

##look at your data
head(data)[1:10]
dim(data)
length(data$L1)
View(data)

## subset data into dataframes for analysis

Coord <- data.matrix(data[,1:2])   # The UTM coordinates of the 500 individuals
Env   <- as.data.frame(data[,3])   # Habitat factor for the 500 individuals
Loci  <- data.matrix(data[,4:103]) # Genotypes at 100 loci for the 500 individuals

# Calculate Moran eigenvector maps, vectors, and values
# --------------------------------
############### Make a map of the spatial relationships between individuals in our dataset
# the first couple steps here tell R how close our different points are to eachother
# this "maps out" the spatial structure in our data

# first, we make a Gabriel graph that maps out all our individuals
# Gabriel graphs show the spatial relationship among points in a dataset on a plane
# look at the ?gabrielneigh help file for more info
# here, we are going to define neighboring points
ind.map      <- graph2nb(gabrielneigh(Coord), sym=TRUE)  # Gabriel graph: neighbor definition

## look at what our Gabriel graph looks like
plot(ind.map, coords=Coord, col=1, pch=16, cex=0.8)

## next, we calculate distances between neighbors with these two commands
disttri <- nbdists(ind.map, Coord)                       # Calculate distances
fdist   <- lapply(disttri, function(x) x^(-1))      # Use inverse distance weights

## now, we are going to weight our Gabriel graph by the distances between neighbors we just calculated
listW   <- nb2listw(ind.map, glist=fdist, style="W")     # Revised spatial weights matrix



############### Habitat maps
# now we're going to look at how our individuals are distributed across space
# load the habitat file
# The habitat file has two different environments, shown here in  yellow and purple
# you can think of this as the same as our beach mouse example in class: 
# imagine yellow is the beach and green is the forest

Habitat <- raster("habitat_data.asc")     # The Habitat ascii file

# plot the habitat data
# these take some time- be patient
plot(Habitat, axes=F, legend=F, box=F,                 # Plot the habitat map
     col=c("#29ae80ff","#FDE725FF"), alpha=0.7) 

# Add the our map of individuals to the habitat map
plot(ind.map, coords=Coord, col="grey44", pch=16, cex=0.6, add=T)


############### Now we conduct MSOD (Moran spectral outlier detection) on all our loci
# first, we compute a Moran's eigenvector map (MEM) for each locus
# MEMs are orthogonal synthetic variables that provide a decomposition of the spatial relationships among individuals based on a spatial weighting matrix
# each MEM axis describes a component of spatial variation among individuals
eigenanalysis     <- scores.listw(listW, MEM.autocor = "all") # MEM eigenanalysis

# next, compile the results of our MEM analysis into a list
mem     <- list(vectors = as.matrix(eigenanalysis), values = attr(eigenanalysis, "values")) # MEM eigenvectors and eigenvalues
mem$values <- mem$values / abs(sum(mem$values))     # Rescale eigenvalues to Moran's I


###### Now that we have determined spatial distribution of individuals, 
# we need to use those spatial relationships to find outlier loci
# to do that, we will look for correlations between loci and our MEM axes 
# the power spectrum (which is the squared correlation) tells us how the allele 
# frequency at one SNP is distributed in space.
# space is defined by the MEMs (they are spatial eigenvectors)
# so we calculate a separate power spectrum for each locus
# if we see a lot of variation in the power spectrum acros MEMs, then 
# the frequency of the allele varies spatially

# Calculate cor.mem, which contains for each locus the vector of its correlations with all MEM axes. 
cor.mem <- cor(Loci, mem$vectors, use="pairwise.complete.obs") # Correlation between allele frequency at each locus and MEM axes     

## we also can calculate the average power spectrum across all loci for each MEM
mean.power.spectrum    <- apply(cor.mem^2, 2, mean)     # power.spectrum = is the average of the squared correlations for each MEM


# Plot power spectra
par(mfrow=c(1,2))
barplot((cor.mem^2)[1,], ylim=c(0, 0.12), main="power spectrum locus 1")  # power spectrum for locus 1
barplot((cor.mem^2)[2,], ylim=c(0, 0.12),main="power spectrum locus 2")  # power spectrum for locus 2
dev.off()

par(mfrow=c(1,2))
barplot(mean.power.spectrum, ylim=c(0, 0.12))  # Average power spectrum (all 100 loci)
barplot(mean.power.spectrum, xlim=c(0,50), ylim=c(0, 0.12))  # Average power spectrum (all 100 loci)
dev.off()

### now that we have a power spectrum for each locus, we want to know which loci show significant spatial varation
# we therefore calculate z-scores for power spectra and assess which ones are significant
# the z scores tell us how much an individual locus deviates from the average power spectrum
# loci that deviate MORE are likely to be outliers

### next few lines calculate the z scores
# If you don't know what a Z score is, don't stress about these
Dev <- sweep(cor.mem^2, 2, mean.power.spectrum, "/") - 1  # Subtract average power spectrum from each locus.
Dev[Dev > 0] <- 0                              # Set positive deviations to zero.
Dev <- apply(Dev, 1, sum)                      # Sum of negative deviations
z <- scale(Dev)                                # Standardize

# Plot z-scores for power spectra
# -------------------------------

## define cutoffs for significance
# here we are going to put three lines on our plot, indicating three different significance levels
dev.off()
plot.new()
cutoffs <- abs(qnorm(c(0.05, 0.01, 0.001)/2))  # Cutoffs (can be modified! But we won't bother here)

## plot z scores
plot(z, ylim=c(-7,5), main= "Z scores for each locus")                     # Plot the z-scores
for(h in 1:length(cutoffs)) {                    # Add lines for the three significance cutoffs
  lines(c(0,100), rep(cutoffs[h],2), lty=h)
  lines(c(0,100), rep(-cutoffs[h],2), lty=h)
}

points(1, z[1], pch=16, col="red")             # color in candidate outlier loci
text(6, -6, "locus1", cex=0.8 )
points(90, z[90], pch=16, col="red")
text(90, -3.5, "locus90", cex=0.8 )
points(95, z[95], pch=16, col="red")
text(97, -3, "locus95", cex=0.8)

### make a list of candidate outlier loci
cutoff.msod <- cutoffs[1]                              # we'll use the lest stringent cutoff
Candidates.msod <- c(1:100)[abs(z)>cutoff.msod]        # makes a list of the loci that are significant at this significance level
Candidates.msod                                        # list of possible outlier loci


###### Moran spectral randomization (MSR) 
# Now that we have identified potential outlier loci, we want to see if those loci are associated with environmental variables
# we will look for correlations betwen allele frequencies at each locus and habitat
# we will control for isolation by distance when we do this
# we control for IBD by removing correlations between individuals that are close to each other in space

# calculate correlations between our MEMs and habitat 
# also calculate correlations between MEMs and x and y coordinates (to control for IBD)
# -----------------------------------------------------------------
Env.locus.cor <- cor(Env, mem$vectors) ## remember, we defined Env at the beginning of the script
IBD.xcoord <- cor(Coord[,1], mem$vectors)
IBD.ycoord <- cor(Coord[,2], mem$vectors)

## we are going to use a permutation test to examine our correlations
# this means that we are going to calcualte OBSERVED correlations between loci and habitat
# we then compare those observed correlations to lots of randomly generated correlations in our dataset

cutoff.msr <- 0.05    # Set a significance cutoff for our permutation test
nPerm <- 199          # Set number of permutations (may choose e.g. 499 or 999)

#### we write a function to perform MSR test
# don't stress about this- this is how we tell R to do the permutations

get.pvalue.msr <- function(r.XV=R.XV, cor.mem=cor.mem, nPerm=199){
  R.XV.rand <- matrix(r.XV, nPerm, ncol(r.XV), byrow=TRUE) 
  R.XV.rand <- R.XV.rand * sample(c(-1,1), length(R.XV.rand), replace=TRUE)
  Cor.obs <- abs(as.vector(cor.mem %*% t(r.XV)))
  Cor.rand <- abs(cor.mem %*% t(R.XV.rand))
  P.values.MSR <- apply((cbind(Cor.obs,Cor.rand) >= Cor.obs), 1, mean)
  P.values.MSR
}

#### MSR test for candidate outlier loci associated with habitat
# this is where we run the  function we made above to do our permutation tests
# we test for correlations between our loci and our environmental variables
# we also test for correlations between loci and distance on the x and y axes (IBD.x and IBD.y)

Env.cor.results <- get.pvalue.msr(r.XV=Env.locus.cor, cor.mem=cor.mem[Candidates.msod,], nPerm=nPerm)
IBD.x.results <- get.pvalue.msr(r.XV=IBD.xcoord, cor.mem=cor.mem[Candidates.msod,], nPerm=nPerm)
IBD.y.results <- get.pvalue.msr(r.XV=IBD.ycoord, cor.mem=cor.mem[Candidates.msod,], nPerm=nPerm)

## we can now print out which loci were significantly associated with each varible
print(paste("Loci significantly associated with environment:", names(Env.cor.results)[Env.cor.results < cutoff.msr]))

print(paste("Loci with IBD on x-axis:", names(IBD.x.results)[IBD.x.results < cutoff.msr]))

print(paste("Loci with IBD on y-axis::", names(IBD.y.results)[IBD.y.results < cutoff.msr]))

###### Now! We can plot the distribution of allele frequencies at differnet loci across our habitat
# these take some time to load- be patient
# run and save each plot individually

# first, define our allele colors
# the homozygotes will be yellow and green, and the heterozygote will be purple

allele_colors <- c("#29ae80ff", "#433880ff" ,"#FDE725FF")          # Set allele colors

### Let's plot our habitat map with our differences in allel frequencies at outlier loci
# first- plot the habitat map
plot(Habitat, axes=F, legend=F, box=F, col=c("#29ae80ff","#FDE725FF"), alpha=0.7) 

# add on the Gabriel graph we made showing proximity between individuals & distribution on landscape
plot(ind.map, coords=Coord, add=T) 

# now, we can color each individual by its allele frequency at Locus 1, our first candidate outlier locus
points(Coord, col = allele_colors[Loci[,1] + 1], pch=20)
title("Allele frequency distribution for Locus 1")

## Ok- now we can plot our second candidate outlier locus
plot(Habitat, axes=F, legend=F, box=F, col=c("#29ae80ff","#FDE725FF"), alpha=0.7)
plot(ind.map, coords=Coord, add=T)
points(Coord, col = allele_colors[Loci[,90] + 1], pch=20)
title("Allele frequency distribution for Locus 90") 

## finally, for comparison, we will plot the spatial distribution of a neutral locus
plot(Habitat, axes=F, legend=F, box=F, col=c("#29ae80ff","#FDE725FF"), alpha=0.7)
plot(ind.map, coords=Coord, add=T)
points(Coord, col = allele_colors[Loci[,2] + 1], pch=20)
title("Allele frequency distribution for Locus 2 (neutral locus)")
