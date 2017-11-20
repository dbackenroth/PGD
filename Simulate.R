library(R.utils)
library(dplyr)
options(stringsAsFactors=F)
source("Utilities.R")


# generate genotypes
# then generate probabilistic genotypes
# Hom_Ref, Hom_Alt, Het
# 

# probabilistic data is three columns
# Prob(data|AA(hom reference)), 
# Prob(data|BB(hom alternative)), 
# Prob(data|AB(heterozygous))


# type is data frame with columns
# Name, Type (Array or Sequencing)
# data should be Father Mother Child1 Child2 Child3
#                  AA    BB      AB    AB      AB
# A is reference, B is alternate


# Prob(observed genotype XX | actual genotype YY)  =0.00006
# Prob(observed genotype XY | genotype genotype XX)=0.00094
# Prob(observed genotype XX or YY | actual genotype XY)  =0.0105

# rows has true genotype, columns has observed genotype

genotypes <- c("AA", "BB", "AB")            

array.error.matrix <- matrix(
  c(0.999, 0.00006, 0.00094,
    0.00006,  0.999, 0.00094,
    0.0105,  0.0105,  0.979),
  byrow=T, nrow=3, ncol=3,
  dimnames=list(genotypes, genotypes)
)

ref.allele.probability <- c(AA=0.998, BB=0.011, AB=0.54)

# for array generates genotypes (AA, BB, AB)
# for sequencing generates counts   num.ref/num.alt
GenerateProbabilisticData <- function(data, subject.data, sequencing.coverage){
  subjects <- subject.data$Name
  types <- subject.data$Type
  obs.data <- matrix(NA, nrow=nrow(data), ncol=length(subjects))
  colnames(obs.data) <- subjects
  genotypes <- rownames(array.error.matrix)
  for (i in 1:length(subjects)){
    subj.type <- types[i]
    subj.name <- subjects[i]
    subj.dat <- data[, subj.name]
    for (true.genotype in genotypes){
      genotype.locs <- which(subj.dat==true.genotype)
      num.variants <- length(genotype.locs)
      if (subj.type=="Array"){
        genotype.probs <- array.error.matrix[true.genotype, ]
        obs.genotypes <- sample(genotypes, 
                                size=num.variants, 
                                prob=genotype.probs, 
                                replace=T)
        obs.data[genotype.locs, subj.name] <- obs.genotypes
      } else if (subj.type=="Sequencing"){
        ref.prob <- ref.allele.probability[true.genotype]
        depth <- rpois(num.variants, sequencing.coverage)
        num.ref <- rbinom(num.variants, size=depth, 
                            prob=ref.prob)
        num.alt <- depth - num.ref
        obs.data[genotype.locs, subj.name] <- 
          paste0(num.ref, "/", num.alt)
      }
    }
  }
  as.data.frame(obs.data)
}

# Recombine should either use constant meiosis.prob or 
# genetic maps
MakeFamily <- function(n.children){
  set.seed(1)
  snp.array.data <- read.table(paste0(data.dir, "SNPArrayData.txt"), header=T)
  #snp.array.data <- snp.array.data[sort(sample(nrow(snp.array.data), 2000)), ]
  positions <- snp.array.data$POS
  recomb <- GetRecombinationProbabilities(positions)
  male.recomb.probs <- recomb$male.recomb.probs
  female.recomb.probs <- recomb$female.recomb.probs
  chr.length <- length(positions)
  #meiosis.prob <- 0.005
  f1 <- SimulateHaplotype(chr.length)
  f2 <- SimulateHaplotype(chr.length)
  m1 <- SimulateHaplotype(chr.length)
  m2 <- SimulateHaplotype(chr.length)
  child.haplotypes <- list()
  for (i in 1:n.children){
    child.haplotypes[[paste0("c", i, ".f")]] <- Recombine(f1, f2, female.recomb.probs)
    child.haplotypes[[paste0("c", i, ".m")]] <- Recombine(m1, m2, male.recomb.probs)
  }
  df <- data.frame(Father=FlipBA(paste0(f1, f2)), 
                   Mother=FlipBA(paste0(m1, m2)))
  for (i in 1:n.children){
    df[, paste0("Child", i)] <- 
      FlipBA(paste0(
        child.haplotypes[[paste0("c", i, ".f")]]$haplotype, 
        child.haplotypes[[paste0("c", i, ".m")]]$haplotype
    ))
  }
  ground.truth <- vector('list', n.children*2)
  for (i in 1:n.children){
    ground.truth[[paste0("C", i, "F")]] <- child.haplotypes[[paste0("c", i, ".f")]]$which.parental.haplotype
    ground.truth[[paste0("C", i, "M")]] <- child.haplotypes[[paste0("c", i, ".m")]]$which.parental.haplotype
  }
  ground.truth <- as.data.frame(do.call('cbind', ground.truth))
  states <- character(chr.length)
  for (i in 2:n.children){
    state <- rep("ni", chr.length)
    father.identical <- ground.truth$C1F==ground.truth[, paste0("C", i, "F")]
    mother.identical <- ground.truth$C1M==ground.truth[, paste0("C", i, "M")]
    state[father.identical] <- "hp"
    state[mother.identical] <- "hm"
    state[mother.identical & father.identical] <- "id"
    if (i==2){
      states <- state
    } else {
      states <- paste0(states, ";", state)
    }
  }
  ground.truth$STATE <- states
  ground.truth$POS <- positions
  df$POS <- positions
  list(ground.truth=ground.truth, data=df, male.recomb.probs=male.recomb.probs, 
       female.recomb.probs=female.recomb.probs)
}

SimulateHaplotype <- function(length){
  sample(c("A", "B"), length, replace=T)
}

Recombine <- function(hap1, hap2, probs){
  l <- length(hap1)
  meioses <- rbinom(l - 1, 1, probs)
  which.chrom <- c(0, cumsum(meioses)) %% 2 + 1
  if (rbinom(1, 1, 0.5)==1){
    which.chrom <- 3 - which.chrom
  }
  haplotypes <- rbind(hap1, hap2)
  index.mat <- cbind(which.chrom, 1:l)
  haplotype <- haplotypes[index.mat]
  list(haplotype=haplotype, which.parental.haplotype=which.chrom)
}