source("HMMNewFixTM_ForBac.R")

options(stringsAsFactors=F)
dir <- "/Users/db2175/NonDData/Data/Carmi/ResultsNovember20/"

# get haplotype blocks for ground truth and for HMM
# are mistakes close to edges

DoWithFB <- function(chr){
  res <- HMMRealData(chr=chr, children=paste0("Child", 2:4), forward.backward=T)
  save(res, file=paste0(dir, "/", chr, ".FB.Rdata"))
}

CompareBlocks <- function(){
  gt <- read.csv(paste0(dir, "/GroundTruthBlocks.csv")) %>% mutate(Type="GroundTruth")
  hmm <- read.csv(paste0(dir, "/HMMBlocks.csv")) %>% mutate(Type="Sequencing")
  both <- bind_rows(gt, hmm) %>% as.data.table()
  res <- list()
  for (child in 2:4){
    dat.subset <- both[Child==child, ]
    setkey(dat.subset, Chromosome, Start, End)
    points <- dat.subset[, .(Points=sort(unique(c(Start, End)))), 
                         by=Chromosome]
    intervals <- points[, .(Start=sort(Points)[1:(length(Points) - 1)], 
                            End=sort(Points)[2:length(Points)]), 
                        by=Chromosome]
    setkey(intervals, Chromosome, Start, End)
    intervals$End <- intervals$End - 1
    overlaps <- foverlaps(intervals, dat.subset)
    overlaps <- overlaps[, .(Chromosome, i.Start, i.End, Value, Child, Parent, 
                             Type)] %>% as.data.frame()
    res[[child]] <- overlaps
  }
  
  all.res <- bind_rows(res) %>% 
    group_by(Child, Chromosome, Parent, i.Start, i.End) %>%
    summarize(Correct=length(unique(Value))==1) %>% ungroup() %>%
    mutate(Length=i.End-i.Start)
  
  # analyze distribution of lengths of right and wrong segments
  Right.Wrong.Segments <- group_by(all.res, Child, Parent, Chromosome) %>%
    arrange(i.Start) %>%
    mutate(Segment=cumsum(c(0, abs(diff(Correct))))) %>%
    group_by(Segment, add=T) %>%
    summarize(Correct=Correct[1], Start=min(i.Start), End=max(i.End)) %>%
    ungroup() %>%
    mutate(Length=End - Start + 1) %>%
    arrange(Child, Parent, Chromosome, Start)
  
  incorrect <- filter(Right.Wrong.Segments, !Correct) %>% arrange(Length)

  summ <- group_by(all.res, Child, Parent) %>%
    summarize(Correct=sum(Length[Correct]), 
              Total=sum(Length)) %>%
    ungroup() %>%
    arrange(Child, Parent) %>%
    mutate(Accuracy=round(Correct/Total, 3))
  
  browser()
  
}

HaplotypeBlocks <- function(){
  l <- list()
  for (type in c("HMM", "GroundTruth")){
    for (child in 2:4){
      for (chr in paste0("chr", 1:22)){
        if (type=="GroundTruth"){
          chr.data <- LoadGroundTruth(chr)
        } else {
          load(paste0(dir, "/", chr, ".noFB.Rdata"))
          chr.data <- res$predictions
        }
        for (maternal in c(T,F)){
          l[[paste0(child, chr, maternal)]] <- 
            BlocksInfo(chr.data=chr.data, child=child, 
                       maternal=maternal)
          l[[paste0(child, chr, maternal)]]$Chromosome <- chr
        }
      }
    }
    all <- bind_rows(l) %>% arrange(Parent, Child, Chromosome, Start)
    if (type=="GroundTruth"){
      f <- paste0(dir, "/GroundTruthBlocks.csv")
    } else {
      f <- paste0(dir, "/HMMBlocks.csv")
    }
    write.csv(all, f, row.names=F)
  }
}

# chr should be in form "chr4" etc...
# rle.exclude specifies the biggest block (measured in number of 
# informative SNPs) that will be excluded
# informative SNPs are SNPs where the applicable parent is heterozygous
# and at least one of the child and the parent is homozygous

# returns data.frame with columns Start, End, NumInfSNPs, Child
#     Parent ('paternal' or 'maternal') and Chromosome ('chr1', etc...)
BlocksInfo <- function(chr.data, maternal=T, child=4, rle.exclude=4){
  all.haps <- HaplotypesAllChildren(chr.data)
  par.type <- ifelse(maternal, "M", "P")
  c1.col <- paste0("Child1", par.type)
  c2.col <- paste0("Child", child, par.type)
  all.haps <- all.haps[!is.na(all.haps[, c1.col]) & !is.na(all.haps[, c2.col]), ]
  all.haps$Comp <- all.haps[, c1.col]==all.haps[, c2.col]
  hets <- all.haps[all.haps[, ifelse(maternal, "Mother", "Father")]=="AB", ]
  rle.r <- rle(cumsum(abs(diff(hets$Comp))))
  hets$Run.Length <- c(1, rep(rle.r$lengths, times=rle.r$lengths))
  filt.hets <- filter(hets, Run.Length>rle.exclude)
  blocks.info <- select(filt.hets, POS, Comp) %>% 
    arrange(POS) %>%
    mutate(BlockID=c(0, cumsum(abs(diff(Comp)))))
  blocks.summ <- group_by(blocks.info, BlockID) %>%
    summarize(Start=min(POS), End=max(POS), NumInfSNPs=n(), 
              Value=Comp[1] + 0) %>%
    mutate(Child=child, 
           Parent=ifelse(maternal, "maternal", "paternal")) %>% select(-BlockID)
  as.data.frame(blocks.summ)
}

# inferred haplotypes for any location where one of the children 
# has an informative site
# if a child isn't informative, NA
# includes columns POS, Child1P, Child1M, ..., CHROM, Child1, Child2, ..., Father
# Mother
# Child1P is A or B (or NA), Child1 is AA, AB or BB
HaplotypesAllChildren <- function(chr.data){
  haps.list <- list()
  for (child in paste0("Child", 1:4)){
    haps.list[[child]] <- HaplotypeOneChild(child, chr.data)
    colnames(haps.list[[child]]) <- c(paste0(child, "P"), 
                                      paste0(child, "M"), 
                                      "POS")
  }
  m <- function(x, y){merge(x, y, all=T)}
  all.haps <- Reduce("m", haps.list)
  all.haps <- merge(all.haps, chr.data)
}

# child name should be "Child1", "Child2", etc...
# returns for positions where one of child, mother or father is heterozygous
# a data.frame with Paternal, Maternal and POS
# Paternal is the paternally transmitted allele, likewise for Maternal
# POS is the position
HaplotypeOneChild <- function(child.name, chr.data){
  #chr.data <- LoadGroundTruth(chr)
  
  colnames(chr.data)[colnames(chr.data)==child.name] <- "Child"
  chr.data <- filter(chr.data, !(Child=="AB" & Mother=="AB" & Father=="AB")) %>% arrange(POS)
  father <- chr.data$Father
  mother <- chr.data$Mother
  child <- chr.data$Child
  
  het.phase.sites <- child=="AB"
  
  paternal <- rep(NA, length(father))
  maternal <- rep(NA, length(mother))
  paternal[father=="BB" & het.phase.sites] <- "B"
  paternal[mother=="AA" & het.phase.sites] <- "B"
  paternal[father=="AA" & het.phase.sites] <- "A"
  paternal[mother=="BB" & het.phase.sites] <- "A"
  
  maternal[het.phase.sites] <- ifelse(paternal[het.phase.sites]=="A", "B", "A")
  
  paternal[child=="AA"] <- "A"
  paternal[child=="BB"] <- "B"
  maternal[child=="AA"] <- "A"
  maternal[child=="BB"] <- "B"
  
  res <- data.frame(Paternal=paternal, Maternal=maternal, POS=chr.data$POS)
  #res[, child.name] <- ch1$Child
  #res$Father <- father
  #res$Mother <- mother
  return(res)
}

AnnotateSegments <- function(){
  rwsegs <- read.csv(paste0(data.dir, "/ResultsJuly12/bb.RightWrongSegments.csv"))
  l <- list()
  k <- 1
  for (chr in paste0("chr", 1:22)){
    load(paste0("/Users/db2175/NonDData/Data/Carmi/RealData/", chr, ".real.data"))
    array.sub.chr <- filter(array.dat, CHROM==chr & Sample %in% c("Father", "Mother")) %>%
      spread(Sample, Gen_Call) %>%
      filter(!(Father=="AB" & Mother=="AB"))
    for (child in 2:4){
      for (parent in c("maternal", "paternal")){
        segs.sub <- filter(rwsegs, Chromosome==chr & Child==child & Parent==parent)
        if (parent=="maternal"){
          array.sub <- filter(array.sub.chr, Mother=="AB")
        } else {
          array.sub <- filter(array.sub.chr, Father=="AB")
        }
        seq.sub <- filter(seq.dat, 
                          Sample==paste0("Child", child) & 
                            CHROM==chr & 
                            (Ref>0 | Alt > 0) & 
                            POS %in% array.sub$POS)
        for (i in 1:nrow(segs.sub)){
          segs.sub$Num.Inf.SNPs.With.Coverage[i] <- sum(seq.sub$POS >= segs.sub$Start[i] &
                                                          seq.sub$POS <= segs.sub$End[i])
        }
        l[[k]] <- segs.sub
        k <- k + 1
      }
    }
  }
  all <- bind_rows(l)
  all <- arrange(all, Chromosome, Child, Start)
  write.csv(all, "SegmentsWithNumberOfSNPs.csv", row.names=F)
}



DoWithoutFB <- function(){
  for (chr in paste0("chr", 1:22)){
    print(chr)
    res <- HMMRealData(chr=chr, children=paste0("Child", 2:4), forward.backward=F, 
                       beta.binomial=T, truncate.limit=NA)
    save(res, file=paste0(dir, "/", chr, ".noFB.Rdata"))
  }
}



