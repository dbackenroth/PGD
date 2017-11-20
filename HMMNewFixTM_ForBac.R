library(R.utils)
library(dplyr)
library(tidyr)
options(stringsAsFactors=F)

genotypes <- c("AA", "AB", "BB")

source("Utilities.R")
source("PrepareData.R")
source("Simulate.R")

basic.states <- c("hm", "hp", "id", "ni")

CombineTogether <- function(){
  for (chrom in paste0("chr", 1:22)){
    f.csv <- paste0(data.dir, "/HmmResultsChildrenTogether.bb/", chrom, "allres.csv")
    chr.data <- LoadGroundTruth(chrom)[, c("CHROM", "POS", "Father", "Mother", "Child1")]
    f <- paste0(data.dir, "/HmmResultsChildrenTogether.bb/", chrom, ".Rdata")
    load(f)
    r <- res$best.genotype[, c("Child2", "Child3", "Child4", "POS")]
    chr.data <- merge(chr.data, r)
    write.csv(chr.data, f.csv, row.names=F)
  }
}


# returns a list of three data frames
#      observations
#       CHROM POS Child1 Child2 [Child...] Father Mother
#       samples with array genotypes have AA, AB, or BB
#       samples with sequencing have Ref/Alt
#      predictions
#       CHROM POS Father Mother Child1 Child 2 [Child ...] PredState
#       each sample has a predicted genotype for each SNP
#       PredState is like hp;hm;hm
#      if forward.backward, marginal
#       CHROM POS hm;hm;hm....ni;ni;ni;Child2SamePat, Child2SameMat
#        first the posterior probability of being in each state
#        then the marginal probability for each child of sharing the same maternal
#        and paternal haplotype as Child1
# children should be Child2, ... (but consecutively starting at Child2)
HMMRealData <- function(children="Child2", 
                        chr="chr21", 
                        forward.backward=F){
  all.data <- GetAllRealData(chr=chr)
  # keep the right children, label data as sequencing or array
  all.data <- filter(all.data, Sample %in% c("Mother", "Father", "Child1", children)) %>%
    mutate(Type=ifelse(Gen_Call %in% genotypes, "Array", "Sequencing"))
  # get rid of sequencing data for Child1
  all.data <- filter(all.data, 
                     !(Sample=="Child1" & Type=="Sequencing"))
  # get rid of array data for other children
  all.data <- filter(all.data, 
                     !(Sample %in% children & Type=="Array"))
  # only keep array positions
  array.positions <- filter(all.data, Sample=="Father")$POS
  all.data <- filter(all.data, POS %in% array.positions)
  
  observed.data <- all.data %>% select(-Type) %>% spread(Sample, Gen_Call)
  hmm.res <- DoHMM(observed.data, 
                   forward.backward=forward.backward)
  hmm.res
}

# observations should be a data frame with columns
#   CHROM, POS, Father, Mother, Child1, Child2 and optionally other
#   consecutively numbered children
#   CHROM      POS Child1 Child4 Father Mother
#   chr21 15206086     AB    0/0     AB     AB
#   chr21 15291929     AB    1/4     AA     AB
#   1 is the number of reference alleles, 4 is the number of alternate alleles
#   there should be data for only one chromosome in observations

DoHMM <- function(observations, forward.backward=F){
  # convert observations to long format
  obs.long <- gather(observations, Sample, GenData, -CHROM, -POS) %>%
    mutate(Type=ifelse(GenData %in% genotypes, "Array", "Sequencing"))
  
  if (length(unique(observations$CHROM))>1){
    stop("DoHMM error: Only one chromosome allowed at a time")
  }
  observations <- arrange(observations, POS)
  recomb <- util.GetRecombinationProbabilities(observations$CHROM[1], observations$POS)
  obs.prob <- util.GenerateDataProbabilities(obs.long)
  
  n.children <- length(grep("Child", colnames(observations)))
  samples <- c("Mother", "Father", paste0("Child", 1:n.children))
  if (!all(samples %in% colnames(observations))){
    stop("DoHMM error:All samples", samples, "should be observations data.frame\n")
  }
  prob.data <- list()
  for (sample in samples){
    s.data <- obs.prob[Sample==sample, ]
    setnames(s.data, genotypes, paste0(sample, "_", genotypes))
    s.data <- as.data.frame(s.data) %>% arrange(POS)
    prob.data[[sample]] <- s.data[, c("POS", paste0(sample, "_", genotypes))]
  }
  prob.data <- Reduce(merge, prob.data)
  hmm.results <- RunHMM(probabilities=prob.data, 
                        male.recomb.probs=recomb$male.recomb.probs, 
                        female.recomb.probs=recomb$female.recomb.probs, 
                        n.children=n.children, forward.backward=forward.backward)
  gen <- hmm.results$best.genotype
  
  PredState <- hmm.results$best.state
  res <- list(observations=observations, 
              predictions=cbind(observations[, c("CHROM", "POS")], gen, PredState))
  if (forward.backward){
    marginal.probs <- hmm.results$forward.m + hmm.results$backward.m
    marginal.probs <- t(apply(marginal.probs, 1, function(x){x-SumProbabilities(x)}))
    marginal.probs <- exp(marginal.probs)
    mp <- as.data.frame(marginal.probs)
    cns <- colnames(mp)
    for (child in 2:n.children){
      child.state <- unlist(lapply(as.list(cns), function(x){
        unlist(strsplit(x, ";"))[child - 1]}))
      pats <- which(child.state %in% c("id", "hp"))
      mp[, paste0("Child", child, "SamePat")] <- apply(mp[, pats], 1, sum)
      mats <- which(child.state %in% c("id", "hm"))
      mp[, paste0("Child", child, "SameMat")] <- apply(mp[, mats], 1, sum)
    }
    res$marginal <- cbind(observations[, c("CHROM", "POS")], mp)
  }
  res
}

# given father's and mother's genotype and number of children
# generates table of all possible inheritance states
# and child genotypes
# along with conditional probability (probability of observing
# given child genotypes conditional on parent genotypes and inheritance state)
MakeCombinations <- function(father="AB", mother="AB", n.children=2){
  segregations <- intToBin(0:(2^n.children - 1))
  # generate data frames of all possible segregations of genotypes to children
  # each segregation is a string like 0010; so first haplotype goes to first, 
  # second and fourth children and second haplotype to third child
  grid <- expand.grid(segregations, segregations, stringsAsFactors=F)
  colnames(grid) <- c("Father", "Mother")
  # add inheritance state labels
  n.minus.1 <- n.children - 1
  for (j in 1:nrow(grid)){
    father.identical <- rep(NA, n.minus.1)
    mother.identical <- rep(NA, n.minus.1)
    for (i in 1:n.minus.1){
      father.identical[i] <- substr(grid$Father[j], 1, 1) == substr(grid$Father[j], i+1, i+1)
      mother.identical[i] <- substr(grid$Mother[j], 1, 1) == substr(grid$Mother[j], i+1, i+1)
    }
    state <- NULL
    for (i in 1:n.minus.1){
      child.state <- "ni"
      if (mother.identical[i]) child.state <- "hm"
      if (father.identical[i]) child.state <- "hp"
      if (mother.identical[i] & father.identical[i]) child.state <- "id"
      state <- c(state, child.state)
    }
    grid$State[j] <- paste0(state, collapse=";")
  }
  
  father.spl <- unlist(strsplit(father, ""))
  mother.spl <- unlist(strsplit(mother, ""))
  for (i in 1:n.children){
    ChildFromFather <- father.spl[as.numeric(substr(grid$Father, i, i)) + 1]
    ChildFromMother <- mother.spl[as.numeric(substr(grid$Mother, i, i)) + 1]
    grid[, paste0("Child", i)] <- paste0(ChildFromFather, ChildFromMother) %>% FlipBA()
  }
  grid <- grid[, c("State", paste0("Child", 1:n.children))]
  # Columns you want to group by
  grp_cols <- c("State", paste0("Child", 1:n.children))
  dots <- lapply(grp_cols, as.symbol)
  u.grid <- group_by_(grid, .dots=dots) %>%
    summarize(n=n()) %>% ungroup() %>% group_by(State) %>%
    mutate(CondProb=n/sum(n))
  u.grid$Father <- father
  u.grid$Mother <- mother
  u.grid
}

# 145 same as 1

AddCode <- function(x){
  x <- mutate(x, Code=paste0(State, ":", "Father", "_", Father, ":", "Mother", "_", Mother)) %>% as.data.frame()
  n.children <- length(grep("Child", colnames(x)))
  for (i in 1:n.children){
    child.num <- paste0("Child", i)
    x$Code <- paste0(x$Code, ":", child.num, "_", x[, child.num])
  }
  x
}

MakeEmissionsProbabilityDictionary <- function(n.children=2, states){
  grid <- expand.grid(genotypes, genotypes, stringsAsFactors=F)
  colnames(grid) <- c("Father", "Mother")
  r <- vector('list', nrow(grid))
  for (i in 1:nrow(grid)){
    r[[i]] <- MakeCombinations(father=grid$Father[i], mother=grid$Mother[i], 
                               n.children=n.children)
  }
  all <- bind_rows(r) %>% AddCode()
  dict.all <- all$CondProb
  names(dict.all) <- all$Code
  # fill in zeros
  grid.list <- list(State=states, Father=genotypes, Mother=genotypes)
  for (i in 1:n.children){
    grid.list[[paste0("Child", i)]] <- genotypes
  }
  f <- expand.grid(grid.list, stringsAsFactors=F) %>% 
    AddCode()
  zero.codes <- f$Code[!f$Code %in% names(dict.all)]
  dict.zero <- rep(0, length(zero.codes))
  names(dict.zero) <- zero.codes
  dict <- c(dict.all, dict.zero)
}

FlipBA <- function(x){
  x[x=="BA"] <- "AB"
  x
}

GetStatesList <- function(n.children){
  basic.states <- c("hm", "hp", "id", "ni")
  states.l <- list()
  for (i in 1:(n.children - 1)){
    states.l[[i]] <- basic.states
  }
  states.grid <- expand.grid(states.l, stringsAsFactors=F)
  for (i in 2:n.children){
    if (i==2){
      states <- states.grid[, i-1]
    } else {
      states <- paste0(states, ";", states.grid[, i-1])
    }
  }
  states
}

TestHMM <- function(n.children=2, n.with.sequencing=1){
  family.dat <- MakeFamily(n.children)
  ground.truth <- family.dat$ground.truth
  dat <- family.dat$data
  subject.data <- data.frame(Name=setdiff(colnames(dat), c("CHROM", "POS")), 
                             Type=c(rep("Array", 2 + n.children - n.with.sequencing), 
                                    rep("Sequencing", n.with.sequencing)))
  samples <- subject.data$Name
  observed.data <- GenerateProbabilisticData(dat, 
                                             subject.data, 
                                             sequencing.coverage=1.8)
  observed.data$CHROM <- 1
  observed.data$POS <- 1:nrow(observed.data)
  sequenced.samples <- subject.data$Name[subject.data$Type=="Sequencing"]
  obs.long <- gather(observed.data, Sample, GenData, -CHROM, -POS) %>%
    mutate(Type=ifelse(Sample %in% sequenced.samples, "Sequencing", "Array")) %>%
    as.data.table()
  obs.prob <- util.GenerateDataProbabilities(obs.long)
  prob.data <- list()
  for (sample in samples){
    s.data <- obs.prob[Sample==sample, ]
    setnames(s.data, genotypes, paste0(sample, "_", genotypes))
    s.data <- as.data.frame(s.data) %>% arrange(POS)
    prob.data[[sample]] <- s.data[, c("POS", paste0(sample, "_", genotypes))]
  }
  prob.data <- Reduce(merge, prob.data)
  hmm.results <- RunHMM(probabilities=prob.data, sample.names=samples, 
                        male.recomb.probs=family.dat$male.recomb.probs, 
                        female.recomb.probs=family.dat$female.recomb.probs, 
                        n.children=n.children)
  state <- hmm.results$best.state
  tt <- table(state, ground.truth$STATE)
  bg <- hmm.results$best.genotype
  for (i in 1:n.children){
    bg.this <- bg[, paste0("Child", i)]
    gt <- dat[, paste0("Child", i)]
    tt <- table(bg.this, gt)
    prop.correct <- sum(diag(tt))/sum(tt)
    cat("Prop.correct child ", i, ":", prop.correct, "\n")
  }
}

# probabilities has columns POS Mother_AA Mother_BB Mother_AB 
#                               Father_AA Father_BB, Father_AB etc...
RunHMM <- function(probabilities, male.recomb.probs, 
                   female.recomb.probs, n.children,
                   forward.backward=F){
  start.time <- proc.time()
  # states is "hm;hm;hm" "hp;hm;hm" etc...   
  states <- GetStatesList(n.children)
  n.states <- length(states)
  n.variants <- nrow(probabilities)
  # returns a dictionary with probabilities
  # hm;hp;hp:Father_AA:Mother_AA:Child1_AA:Child2_AA:Child3_AA:Child4_AA 
  # 1 
  # probability of observing given child genotypes given the parent genotypes and the
  # inheritance state
  emissions.prob.dict <- MakeEmissionsProbabilityDictionary(n.children=n.children, 
                                                            states=states)
  emissions.probs <- matrix(NA, nrow=n.variants, ncol=n.states)
  colnames(emissions.probs) <- states
  start.pos <- gregexpr(pattern ='Father', names(emissions.prob.dict)[1])[[1]][1]
  # all combinations of genotypes
  combinations <- unique(substr(names(emissions.prob.dict), start.pos, nchar(names(emissions.prob.dict))))
  n.combinations <- length(combinations)
  comb.probabilities <- matrix(NA, nrow=n.variants, ncol=n.combinations)
  colnames(comb.probabilities) <- combinations
  # for each combination multiply the probabilities of all the
  # component genotypes
  # now we know for each combination of genotypes that probability of the
  # data given that genotype
  for (comb in combinations){
    probs.to.multiply <- unlist(strsplit(comb, ":", fixed=T))
    ptm <- probabilities[, probs.to.multiply]
    comb.probabilities[, comb] <- apply(ptm, 1, prod)
  }
  for (i in 1:n.states){
    state <- states[i]
    state.emissions.prob.dict <- emissions.prob.dict[substr(names(emissions.prob.dict), 1, start.pos-2)==state]
    state.emissions.prob.dict <- state.emissions.prob.dict[state.emissions.prob.dict > 0]
    state.probs <- rep(0, n.variants)
    for (j in 1:length(state.emissions.prob.dict)){
      state.code <- names(state.emissions.prob.dict)[j]
      genotype.code <- substr(state.code, start.pos, nchar(state.code))
      genotype.probs <- comb.probabilities[, genotype.code]
      state.probs <- state.probs + genotype.probs * state.emissions.prob.dict[j]
    }
    emissions.probs[, i] <- log(state.probs)
  }
  if (F){
    save(emissions.probs, probabilities, 
         file=paste0(data.dir, "/Temp/temp.Rdata"))
  }
  initial.state.prob <- log(rep(0.25, length(states)))
  cat("Starting Viterbi\n")
  viterbi <- DoViterbi(emissions.probs, male.recomb.probs=male.recomb.probs,
                       female.recomb.probs=female.recomb.probs,
                       initial.state.prob=initial.state.prob, states=states, 
                       n.children=n.children)
  if (forward.backward){
    cat("Starting forward pass\n")
    forward.m <- GetForwardMatrix(emission.probs=emissions.probs, 
                                  male.recomb.probs=male.recomb.probs, 
                                  female.recomb.probs=female.recomb.probs, 
                                  initial.state.prob=initial.state.prob, 
                                  states=states, n.children=n.children)
    cat("Starting backward pass\n")
    backward.m <- GetBackwardMatrix(emission.probs=emissions.probs, 
                                    male.recomb.probs=male.recomb.probs, 
                                    female.recomb.probs=female.recomb.probs, 
                                    initial.state.prob=initial.state.prob, 
                                    states=states, n.children=n.children)
    colnames(forward.m) <- colnames(backward.m) <- colnames(emissions.probs)
  } else {
    forward.m <- NULL
    backward.m <- NULL
  }
  state <- states[viterbi$viterbi]
  # get most likely genotypes
  # use comb.probabilities
  # 
  u.states <- unique(state)
  emissions.positive.combinations <- names(emissions.prob.dict[emissions.prob.dict > 0])
  best.genotype <- rep(NA, n.variants)
  for (s in u.states){
    which.variants <- which(state==s)
    state.combinations <- emissions.positive.combinations[substr(emissions.positive.combinations, 1, start.pos - 2)==s]
    s.c.genotypes <- substr(state.combinations, start.pos, nchar(state.combinations))
    s.comb.probabilities <- comb.probabilities[which.variants, s.c.genotypes]
    s.comb.probs.colnames <- colnames(s.comb.probabilities)
    temp.max <- apply(s.comb.probabilities, 1, which.max)
    max.state <- s.comb.probs.colnames[temp.max]
    best.genotype[which.variants] <- max.state
  }
  sample.names <- unlist(strsplit(unlist(strsplit(best.genotype[1], ":")), "_"))[c(T,F)]
  best.genotype <- data.frame(Temp=best.genotype) %>% separate_("Temp", sample.names, sep=":")
  best.genotype <- apply(best.genotype, 2, function(x){unlist(strsplit(x, "_"))[c(F,T)]}) %>% as.data.frame()
  return(list(best.genotype=best.genotype, best.state=state, 
              forward.m=forward.m, backward.m=backward.m))
}

TransitionM.PowerMatrices <- function(states){
  N.states <- length(states)
  mother.power.m <- father.power.m <- matrix(1, nrow=N.states, ncol=N.states, 
                                             dimnames=list(states, states))
  diag(mother.power.m) <- diag(father.power.m) <- NA
  for (i in 1:N.states){
    for (j in 1:N.states){
      if (i==j) next
      father.power <- 0
      mother.power <- 0
      from <- unlist(strsplit(states[i], ";"))
      to <- unlist(strsplit(states[j], ";"))
      prob <- 1
      for (k in 1:length(from)){
        if (from[k]=="hm" & to[k]=="hp" | 
            from[k]=="hp" & to[k]=="hm" | 
            from[k]=="id" & to[k]=="ni" | 
            from[k]=="ni" & to[k]=="id"){
          father.power <- father.power + 1
          mother.power <- mother.power + 1
        }
        if (from[k]=="hm" & to[k]=="ni" |
            from[k]=="ni" & to[k]=="hm" |
            from[k]=="hp" & to[k]=="id" |
            from[k]=="id" & to[k]=="hp"){
          mother.power <- mother.power + 1
        }
        if (from[k]=="hp" & to[k]=="ni" |
            from[k]=="ni" & to[k]=="hp" |
            from[k]=="hm" & to[k]=="id" |
            from[k]=="id" & to[k]=="hm"){
          father.power <- father.power + 1
        }
      }
      mother.power.m[i, j] <- mother.power
      father.power.m[i, j] <- father.power
    }
  }
  return(list(father=father.power.m, mother=mother.power.m))
}

TransitionM2 <- function(father.matrix, mother.matrix, male.recomb.prob, 
                         female.recomb.prob, n.children){
  male.matrix <- 
    male.recomb.prob^father.matrix * (1 - male.recomb.prob)^(n.children - father.matrix) + 
    male.recomb.prob^(n.children - father.matrix) * (1 - male.recomb.prob)^father.matrix
  female.matrix <- 
    female.recomb.prob^mother.matrix * (1 - female.recomb.prob)^(n.children - mother.matrix) + 
    female.recomb.prob^(n.children - mother.matrix) * (1 - female.recomb.prob)^mother.matrix
  tm <- male.matrix * female.matrix
  for (i in 1:ncol(tm)){
    tm[i, i] <- 1 - sum(tm[i, ], na.rm=T)
  }
  return(log(tm))
}

DoViterbi <- function(emission.probs.m, male.recomb.probs, female.recomb.probs, 
                      initial.state.prob, states, n.children){
  tm.matrices <- TransitionM.PowerMatrices(states)
  NUM.STATES <- length(initial.state.prob)
  N <- nrow(emission.probs.m)
  viterbi.matrix <- matrix(NA, nrow=N, ncol=NUM.STATES)
  viterbi.pointers <- matrix(NA, nrow=N, ncol=NUM.STATES)
  viterbi.matrix[1, ] <- initial.state.prob + emission.probs.m[1, ]
  for (i in 2:N) {
    temp.matrix <- viterbi.matrix[i - 1, ] + 
      TransitionM2(father.matrix=tm.matrices$father, 
                   mother.matrix=tm.matrices$mother,
      #TransitionM(states=states, 
                  male.recomb.prob=male.recomb.probs[i-1], 
                  female.recomb.prob=female.recomb.probs[i-1], 
      n.children=n.children)
    viterbi.matrix[i, ] <- apply(temp.matrix, 2, max)
    emission.probs <- c(emission.probs.m[i,])
    dim(emission.probs) <- c(NUM.STATES, 1)
    viterbi.matrix[i, ] <- viterbi.matrix[i, ] + emission.probs
    viterbi.pointers[i, ] <- apply(temp.matrix, 2, which.max)
  }
  viterbi.states = vector(length = N)
  viterbi.states[N] = which.max(viterbi.matrix[N, ])
  max.prob <- max(viterbi.matrix[N, ])
  for (i in (N - 1):1) {
    viterbi.states[i] <- viterbi.pointers[i + 1, viterbi.states[i + 1]]
  }
  print(max.prob)
  return(list(mp=max.prob, viterbi=viterbi.states))
}

# adds two log-space probabilities using the identity
# log (p1 + p2) = log p1 + log(1 + exp(log p2 - log p1))
AddTwoProbabilities <- function(x, y){
  if (is.infinite(x)) return (y)
  if (is.infinite(y)) return (x)
  sum.probs <- max(x, y) + log1p(exp(-abs(x - y)))
}

# adds multiple log-space probabilities
SumProbabilities <- function(x){
  sum.probs <- x[1]
  for (i in 2:length(x)){
    sum.probs <- AddTwoProbabilities(sum.probs, x[i])
  }
  return(sum.probs)
}

# get the forward probabilities
GetForwardMatrix <- function(emission.probs.m, male.recomb.probs, 
                             female.recomb.probs, 
                             initial.state.prob, states, n.children){
  tm.matrices <- TransitionM.PowerMatrices(states)
  NUM.STATES <- length(initial.state.prob)
  N <- nrow(emission.probs.m)
  forward.matrix <- matrix(NA, nrow=N, ncol=NUM.STATES)   # matrix to hold forward probabilities
  forward.matrix[1, ] <- initial.state.prob + emission.probs.m[1, ]
  for (i in 2:N){
    # compute matrix with probability we were in state j and are now in state i
    # in temp.matrix[j, i] (ignoring emission of current token)
    temp.matrix <- forward.matrix[i - 1, ] + 
      TransitionM2(father.matrix=tm.matrices$father, 
                   mother.matrix=tm.matrices$mother,
                   male.recomb.prob=male.recomb.probs[i-1], 
                   female.recomb.prob=female.recomb.probs[i-1], 
                   n.children=n.children)
    # find the probability that we are in each of the three states
    sum.probs <- apply(temp.matrix, 2, SumProbabilities)
    forward.matrix[i, ] <- sum.probs + emission.probs.m[i, ]
  }  
  return(forward.matrix)  
}

# get the backward probabilities
GetBackwardMatrix <- function(emission.probs.m, male.recomb.probs, 
                              female.recomb.probs, 
                              initial.state.prob, states, n.children){
  tm.matrices <- TransitionM.PowerMatrices(states)
  NUM.STATES <- length(initial.state.prob)
  N <- nrow(emission.probs.m)
  
  backward.matrix <- matrix(NA, nrow=N, ncol=NUM.STATES)   # matrix to hold backward probabilities
  
  backward.matrix[N, ] <- rep(0, NUM.STATES)
  for (i in (N - 1):1){
    temp.matrix <- TransitionM2(father.matrix=tm.matrices$father, 
                                mother.matrix=tm.matrices$mother,
                                male.recomb.prob=male.recomb.probs[i], 
                                female.recomb.prob=female.recomb.probs[i], 
                                n.children=n.children) + 
      matrix(backward.matrix[i + 1, ], NUM.STATES, NUM.STATES, byrow=T) +
      matrix(emission.probs.m[i+1, ], NUM.STATES, NUM.STATES, byrow=T)
    backward.matrix[i, ] <- apply(temp.matrix, 1, SumProbabilities)
  }  
  final.prob <- backward.matrix[1, ] + emission.probs.m[1, ] + initial.state.prob
  return(backward.matrix)  
}




# CombineAll <- function(){
#   for (chrom in paste0("chr", 1)){
#     f.csv <- paste0(data.dir, "/HmmResults.bb/", chrom, "allres.csv")
#     chr.data <- LoadGroundTruth(chrom)[, c("CHROM", "POS", "Father", "Mother", "Child1")]
#     for (child in paste0("Child", 2:4)){
#       f <- paste0(data.dir, "/HmmResults.bb/", chrom, child, ".Rdata")
#       load(f)
#       r <- res$best.genotype[, c(child, "POS", "CHROM")]
#       chr.data <- merge(chr.data, r)
#     }
#     write.csv(chr.data, f.csv, row.names=F)
#   }
# }