library(dplyr)
library(ggplot2)
library(gridExtra)

source("PrepareReadCountData.R")    # provides GCCorrect()
source("Utilities.R")               # provides mydnbinomlog

PlotSamples <- function(){
  samples <- paste0("Child", 1:4)
  override.factors <- c(NA, NA, NA, NA)
  hp <- list()
  for (i in 1:4){
    l <- RDRunTwice(sample=samples[i], override.factor=override.factors[i])
    hp[[i]] <- l$plot
  }
  pdf("Results/FourChildren.pdf", height=5, width=8)
  grid.arrange(hp[[1]], hp[[2]], 
               hp[[3]], hp[[4]], 
               nrow=4)
  dev.off()
}

# runs the HMM twice
# the second time, calculates an adjustment factor based on the median
# ratio of read count for inferred normal regions to the expected 
# read count for those regions
RDRunTwice <- function(sample="Child3", override.factor=NA){
  if (!is.na(override.factor)){
    mult.factor <- override.factor
  } else mult.factor <- 1
  # GetGCCorrect() 
  dat <- GetGCCorrect() %>% filter(Normal | Sample==sample)
  r1 <- RDRunHMM(dat=dat, sample=sample, mult.factor=mult.factor)
  normal <- r1$l %>% filter(State=="TwoCopies")
  adj.factor <- median(normal$NormGC/(normal$WindowMedian/mult.factor))
  print(adj.factor)
  r2 <- RDRunHMM(dat=dat, mult.factor=adj.factor, sample=sample)
  p <- PlotHMMResults(r2$l)
  list(hmm.res=r2, plot=p)
}

RDRunHMM <- function(dat, sample, mult.factor){
  
  dat <- select(dat, chr, start, end, Sample, NormGC, Normal)
  stopifnot(sample %in% dat$Sample)
  stopifnot(is.numeric(mult.factor) & length(mult.factor)==1)
  
  # calculate median and SD for each window, using only normal samples
  window.data <- filter(dat, Normal) %>% 
    group_by(chr, start) %>%
    summarize(WindowMedian=median(NormGC), 
              WindowSD=sd(NormGC))
  dat <- merge(dat, window.data)
  dat$WindowMedian <- dat$WindowMedian * mult.factor
  #sd.mean.fit1 <- lm(WindowSD ~ WindowMedian + 0, data=window.data)
  dat <- dat %>% filter(Sample==sample) %>%
    mutate(Deviation=abs(NormGC-WindowMedian))
  sd.mean.fit2 <- lm(Deviation ~ WindowMedian + 0, data=dat)
  emps <- GetEmissionProbabilities(dat %>% filter(Sample==sample), sd.mean.fit2)
  tm <- RDTransitionMatrix()
  ip <- RDInitialStateProb()
  l <- list()
  all.emps <- list()
  for (Chr in 1:22){
    chr.dat <- filter(emps, chr==Chr) %>% arrange(start)
    chr.emps <- select(chr.dat, Del2, Del1, Normal, Dup1, Dup2) %>% as.matrix()
    vit <- rd.DoViterbi(emission.probs.m=chr.emps, initial.state.prob=ip, transition.m=tm)
    State <- c("ZeroCopies", "OneCopy", "TwoCopies", "ThreeCopies", "FourCopies")[vit$viterbi]
    l[[Chr]] <- cbind(chr.dat, State)
    all.emps[[Chr]] <- 
      select(chr.dat, chr, start, end, Sample, Del2, Del1, Normal, Dup1, Dup2)
  }
  list(l=bind_rows(l), emps=bind_rows(all.emps))
}

PlotHMMResults <- function(l){
  
  l <- select(l, chr, start, end, State, NormGC, WindowMedian, Sample)
  
  l <- l %>% mutate(gen.start=GetGenomicCoordinates(chr, start), 
               gen.end=GetGenomicCoordinates(chr, end))
  midcoords <- GetMidChromosomalCoordinates()
  cnvs.combined <- group_by(l, chr) %>%
    select(chr, gen.start, gen.end, State) %>%
    arrange(gen.start) %>%
    mutate(DiffState=c(0, !State[2:n()]==State[1:(n()-1)])) %>%
    mutate(Segment=cumsum(DiffState)) %>%
    group_by(chr, Segment) %>%
    summarize(seg.start=min(gen.start), 
              seg.end=max(gen.end), 
              State=State[1], 
              NumStates=length(unique(State))) %>%
    filter(!State=="TwoCopies")
  print(as.data.frame(cnvs.combined))
  p <- ggplot(NULL) + 
    geom_point(data=l, 
               aes(x=gen.start, y=NormGC/WindowMedian, group=chr, col=State), size=0.25) +
    theme_bw() + 
    geom_rect(data=cnvs.combined, aes(xmin=seg.start, xmax=seg.end, ymin=0, ymax=2, fill=State), alpha=0.25) + 
    geom_hline(yintercept=1, alpha=0.1) + 
    scale_y_continuous(breaks=c(0.5, 1, 1.5), limit=c(0, 2)) + 
    ylab("RC/median") + 
    scale_x_continuous(breaks=midcoords, labels=1:22) + 
    facet_wrap(~Sample, ncol=1) + 
    scale_fill_manual("", guide=F, values=c(ZeroCopies="yellow", 
                                   OneCopy="red", 
                                   ThreeCopies="blue", 
                                   FourCopies="purple")) + 
    scale_color_manual("", guide=F, values=c(ZeroCopies="yellow", 
                                             OneCopy="red", 
                                             TwoCopies="black",
                                             ThreeCopies="blue", 
                                             FourCopies="purple")) + 
    xlab("") + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
}


GetEmissionProbabilities <- function(dat, sd.mean.fit){
  dat <- select(dat, Sample, chr, start, end, NormGC, WindowMedian)
  CalcLogProbs <- function(x, mean, fit){
    mydnbinomlog(x=round(x), mean=mean, sd=predict(fit, newdata=data.frame(WindowMedian=mean)))
  }
  dat$Normal <- CalcLogProbs(x=dat$NormGC, mean=dat$WindowMedian, fit=sd.mean.fit)
  dat$Del1 <- CalcLogProbs(x=dat$NormGC, mean=dat$WindowMedian / 2, fit=sd.mean.fit)
  dat$Del2 <- CalcLogProbs(x=dat$NormGC, mean=dat$WindowMedian/10, fit=sd.mean.fit)
  dat$Dup1 <- CalcLogProbs(x=dat$NormGC, mean=dat$WindowMedian * 3 / 2, fit=sd.mean.fit)
  dat$Dup2 <- CalcLogProbs(x=dat$NormGC, mean=dat$WindowMedian * 2, fit=sd.mean.fit)
  dat
}

RDTransitionMatrix <- function(){
  transition.m <- matrix(c(NA, 1/300, 1/300, 0, 0,
                           1/8000, NA, 1/300, 0, 0, 
                           0, 1/8000, NA, 1/8000, 0, 
                           0, 0, 1/300, NA, 1/8000,
                           0, 0, 1/300, 1/300, NA), 
                          byrow=T, nrow=5)
  for (i in 1:5){
     transition.m[i,i] <- 1 - sum(transition.m[i, ], na.rm=T)
  } 
  log(transition.m)
}

RDInitialStateProb <- function(){
  log(c(1/4000000, 1/2000, 1998/2000, 1/2000, 1/4000000))
}


rd.DoViterbi <- function(emission.probs.m,
                      initial.state.prob, 
                      transition.m){
  NUM.STATES <- length(initial.state.prob)
  N <- nrow(emission.probs.m)
  viterbi.matrix <- matrix(NA, nrow=N, ncol=NUM.STATES)
  viterbi.pointers <- matrix(NA, nrow=N, ncol=NUM.STATES)
  viterbi.matrix[1, ] <- initial.state.prob + emission.probs.m[1, ]
  for (i in 2:N) {
    temp.matrix <- viterbi.matrix[i - 1, ] + 
      transition.m
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
  return(list(mp=max.prob, viterbi=viterbi.states))
}

# TestGGC <- function(){
#   ll <- RunAllChildren()
#   ll$gen.start <- GetGenomicCoordinates(ll$chr, ll$start)
#   ll$gen.end <- GetGenomicCoordinates(ll$chr, ll$end)
#   midcoords <- GetMidChromosomalCoordinates()
#   cnvs <- filter(ll, State %in% c("Deletion", "Duplication"))
#   p <- ggplot(NULL) + 
#     geom_point(data=ll, 
#                aes(x=gen.start, y=NormDepth/WindowMedian, group=chr), size=0.25) +
#     theme_bw() + 
#     geom_rect(data=cnvs, aes(xmin=gen.start, xmax=gen.end, ymin=0.3, ymax=1.7, fill=State), alpha=0.45) + 
#     geom_hline(yintercept=1, alpha=0.1) + 
#     scale_y_continuous(breaks=c(0.5, 1, 1.5), limit=c(0.3, 1.7)) + 
#     ylab("Depth/median") + 
#     scale_x_continuous(breaks=midcoords, labels=1:22) + 
#     facet_wrap(~Sample, ncol=1) + 
#     scale_fill_manual("", values=c(Deletion="red", Duplication="blue")) + 
#     xlab("Chromosome") + 
#     theme(panel.grid.major = element_blank(), 
#           panel.grid.minor = element_blank())
#   ggsave("ReadDepthHMM.jpeg", width=6, height=6)
# }

# RunAllChildren <- function(){
#   l <- list()
#   for (child in paste0("Child", 1:4)){
#     l[[child]] <- RDRunHMM(child=child)
#   }
#   l <- bind_rows(l)
# }
# 
# a <- function(){
#   mult.factors <- seq(0.7, 1.3, by=0.01)
#   mps <- rep(NA, length(mult.factors))
#   for (i in 1:length(mult.factors)){
#     mps[i] <- RDRunHMM(mult.factor=mult.factors[i])
#   }
#   browser()
# }


# grab_grob <- function(){
#   grid.echo()
#   grid.grab()
# }