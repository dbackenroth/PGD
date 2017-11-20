library(data.table)
library(dplyr)
library(tidyr)
library(emdbook)


cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

options(stringsAsFactors=F)

data.dir <- "/Users/db2175/NonDData/Data/Carmi/"

genotypes <- c("AA", "BB", "AB")   

array.error.matrix <- matrix(
  c(0.999, 0.00006, 0.00094,
    0.00006,  0.999, 0.00094,
    0.0105,  0.0105,  0.979),
  byrow=T, nrow=3, ncol=3,
  dimnames=list(genotypes, genotypes)
)

ref.allele.probability <- c(AA=0.998, BB=0.011, AB=0.54)

bb.params <- c(prob=0.554, theta=0.473)

mydnbinomlog <- function(x, mean, sd){
  size <- mean^2 - (sd^2 - mean)
  dnbinom(x=x, size=size, mu=mean, log=T)
}

# pd should be a data.table or data.frame with columns Type (Sequencing or Array)
# and GenData
# returns 
util.GenerateDataProbabilities <- function(pd){
  pd <- as.data.table(pd)
  seq.data <- pd[Type=="Sequencing" & !is.na(GenData), ]
  array.data <- pd[Type=="Array" & !is.na(GenData), ]
  na.data <- pd[is.na(GenData), ]
  # want to calculate probability of observed data given
  # each genotype
  for (true.genotype in genotypes){
    array.genotype.probs <- array.error.matrix[true.genotype, ]
    array.data[, (true.genotype):=array.genotype.probs[array.data$GenData]]
    spl <- strsplit(seq.data$GenData, split="/", fixed=T)
    ref <- as.numeric(unlist(spl)[c(T,F)])
    alt <- as.numeric(unlist(spl)[c(F,T)])
    depth <- ref + alt
    
    seq.data[, (true.genotype):=dbetabinom(x=ref, size=depth, 
                                             #prob=bb.params['prob'], 
                                             prob=ref.allele.probability[true.genotype],
                                             theta=bb.params['theta'])]
    #} else {
    #  seq.data[, (true.genotype):=dbinom(ref, depth, prob=ref.allele.probability[true.genotype])]
    #}
  }
  array.data[is.na(AA), `:=`(AA=1, AB=1, BB=1)]
  na.data[, `:=`(AA=1, AB=1, BB=1)]
  pd <- rbind(seq.data, array.data, na.data)
}


# input and output are both in this form:
#CHROM      POS Father Mother Child1 Child2
# chr21 15206086     AB     AB     AB    3/0
# chr21 15214708     AB     AA     AB    0/0
util.Downsample <- function(observations, truncate.limit){
  pd <- gather(observations, Sample, GenData, -CHROM, -POS) %>%
    mutate(Type=ifelse(GenData %in% genotypes, "Array", "Sequencing")) %>%
    filter(!is.na(GenData)) %>%
    as.data.table()
  if (nrow(pd[Type=="Sequencing", ])>0){
    pd[Type=="Sequencing", `:=`(
      Ref.Count=as.numeric(
        unlist(strsplit(GenData, "/", fixed=T))[c(T,F)]), 
      Alt.Count=as.numeric(
        unlist(strsplit(GenData, "/", fixed=T))[c(F,T)]))]
    pd[, `:=`(Total.Count=Ref.Count+Alt.Count)]
    if (!is.na(truncate.limit)){
      pd[Total.Count > truncate.limit, 
       `:=`(New.Ref.Count=rhyper(nn=.N, m=Ref.Count, 
                                 n=Alt.Count, k=truncate.limit))]
      pd[Total.Count > truncate.limit, 
       `:=`(New.Alt.Count = truncate.limit - New.Ref.Count)]
      pd[Total.Count > truncate.limit, 
       `:=`(GenData=paste0(New.Ref.Count, "/", New.Alt.Count))]
      pd[, `:=`(Ref.Count=NULL, Alt.Count=NULL, Total.Count=NULL, 
              New.Ref.Count=NULL, New.Alt.Count=NULL)]
    }
  }
  observations <- select(pd, CHROM, POS, Sample, GenData) %>% 
    spread(Sample, GenData) %>% 
    arrange(POS)
}

util.GetRecombinationProbability <- function(cM){
  cm.diffs <- cM[2:length(cM)] - 
    cM[1:(length(cM) - 1)]
  probs <- (1 - exp(-2*cm.diffs/100))/2
}

# positions should be a numeric vector
# returns probabilities for males and female of recombination between
# each successive pair of positions
# in list with elements male.recomb.probs and female.recomb.probs
util.GetRecombinationProbabilities <- function(chr, positions){
  gmap.male <- read.table(paste0("/Users/db2175/Dropbox/Daniel/Biostatistics/Carmi/Refined_genetic_map_b37/male_", chr, ".txt"), header=T, stringsAsFactors=F)
  gmap.female <- read.table(paste0("/Users/db2175/Dropbox/Daniel/Biostatistics/Carmi/Refined_genetic_map_b37/female_", chr, ".txt"), header=T, stringsAsFactors=F)
  male.interp <- approx(gmap.male$pos, gmap.male$cM, xout=positions)
  female.interp <- approx(gmap.female$pos, gmap.female$cM, xout=positions)
  male.recomb.probs <- util.GetRecombinationProbability(male.interp$y)
  female.recomb.probs <- util.GetRecombinationProbability(female.interp$y)
  male.recomb.probs[is.na(male.recomb.probs)] <- min(male.recomb.probs, na.rm=T)
  female.recomb.probs[is.na(female.recomb.probs)] <- min(female.recomb.probs, na.rm=T)
  return(list(male.recomb.probs=male.recomb.probs, 
              female.recomb.probs=female.recomb.probs))
}

# from in rows
# chr is 1,2,3...
# pos is numeric
GetGenomicCoordinates <- function(chr, pos, buffer=25e+06){
  chrom.sizes <- read.table("hg19.chrom.sizes", 
                            col.names=c("chr", "length")) %>%
    filter(chr %in% paste0("chr", 1:22)) %>% 
    mutate(chr=as.numeric(gsub("chr", "", chr))) %>%
    arrange(chr)
  chrom.sizes$start <- c(0, cumsum(chrom.sizes$length[1:21] + buffer))
  coord <- chrom.sizes$start[chr] + pos
}

GetMidChromosomalCoordinates <- function(buffer=25e+06){
  chrom.sizes <- read.table("hg19.chrom.sizes", 
                            col.names=c("chr", "length")) %>%
    filter(chr %in% paste0("chr", 1:22)) %>% 
    mutate(chr=as.numeric(gsub("chr", "", chr))) %>%
    arrange(chr)
  GetGenomicCoordinates(chr=1:22, pos=chrom.sizes$length / 2, buffer=buffer)
}