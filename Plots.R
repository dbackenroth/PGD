library(ggplot2)

dir <- "/Users/db2175/NonDData/Data/Carmi/ResultsNovember20/"

MakeChromosomePlots <- function(){
  for (chr.no in 21){
    chr <- paste0("chr", chr.no)
    PlotAllMarginals(chr=chr, start=0, end=400000000, plot.title=chr, gene.boundaries=F)
    ggsave(paste0("Results/", chr, ".haplotypes.pdf"), width=15-chr.no/3, height=6)
    ggsave(paste0("Results/", chr, ".haplotypes.jpeg"), width=15-chr.no/3, height=6)
  }
}

MakePlots <- function(){
  PlotAllMarginals(chr="chr21", start=1000000, end=44000000, left.buffer=2000000, right.buffer=2000000, plot.title="chr21", gene.boundaries=F)
  ggsave("")
}

PlotAllMarginals <- function(chr="chr13", start=32889617, end=32973809, left.buffer=1000000, right.buffer=1000000, plot.title="BRCA2", gene.boundaries=T){
  dat <- list()
  informatives <- list()
  for (child in paste0("Child", 2:4)){
    m <- MarginalsInfo(chr=chr, child=child, start=start, end=end, left.buffer=left.buffer, right.buffer=right.buffer)
    dat[[child]] <- m$dat
    informatives[[child]] <- m$informatives
  }
  dat <- bind_rows(dat)
  informatives <- bind_rows(informatives)
  ground.truth.blocks <- read.csv(paste0(dir, "/GroundTruthBlocks.csv")) %>% filter(Chromosome==chr, Value==1) %>% mutate(ymin=0, ymax=1) %>% mutate(Parent=ifelse(Parent=="maternal", "Maternal", "Paternal"), Child=paste0("Child", Child))
  p <- ggplot(NULL) + geom_line(data=dat, aes(x=POS, y=Marginal)) + theme_bw() + scale_y_continuous(limits=c(-0.05, 1)) + geom_segment(data=informatives, aes(x=x, xend=xend, y=y, yend=yend)) + facet_wrap(~Parent + Child, nrow=2) + ylab("Marginal probability of sharing of Child 1 haplotype") + xlab("Position") + ggtitle(plot.title) + geom_rect(data=ground.truth.blocks, aes(xmin=Start, xmax=End, ymin=ymin, ymax=ymax), fill="green", alpha=0.2)
  if (gene.boundaries){
    p <- p + geom_vline(xintercept=c(start, end), color="blue")
  }
  print(p)
}


# plot informative SNPs for that parent with read depth
# don't use ground truth
# 
MarginalsInfo <- function(chr, child, start, end, 
                          left.buffer=2000000, right.buffer=2000000){
  left <- start - left.buffer
  right <- end + right.buffer
  load(paste0(dir, "/", chr, ".FB.Rdata"))
  marg <- res$marginal[, c("CHROM", "POS", paste0(child, "SameMat"), paste0(child, "SamePat"))] %>%
    filter(POS >= left & POS <= right)
  obs <- res$observations[, c("CHROM", "POS", "Father", "Mother", child)] %>%
    filter(POS >= left & POS <= right) %>% separate_(child, 
                                                     c("Ref", "Alt"), 
                                                     sep="/", convert=T)
  pred <- res$predictions[, c("CHROM", "POS", child)]
  obs <- merge(obs, pred)
  GetInformativeReads <- function(obs, parent, child){
    InfReads <- rep(0, nrow(obs))
    for (i in 1:nrow(obs)){
      sel.parent.obs <- obs[i, parent]
      other.parent.obs <- obs[i,  ifelse(parent=="Mother", "Father", "Mother")]
      if (obs[i, child]=="AB"){
        if (sel.parent.obs=="AB" & other.parent.obs=="AB"){
          InfReads[i] <- 0
        } else {
          if (sel.parent.obs=="AA"){
            InfReads[i] <- obs[i, "Ref"]
          } else if (sel.parent.obs=="BB"){
            InfReads[i] <- obs[i, "Alt"]
          } else if (sel.parent.obs=="AB"){
            if (other.parent.obs=="AA"){
              InfReads[i] <- obs[i, "Alt"]
            } else if (other.parent.obs=="BB"){
              InfReads[i] <- obs[i, "Ref"]
            }
          }
        }
      } else {
        InfReads[i] <- 0
      }
    }
    InfReads
  }
  obs$MotherInf <- GetInformativeReads(obs, "Mother", child)
  obs$FatherInf <- GetInformativeReads(obs, "Father", child)
  dat <- merge(obs, marg) 
  mip <- data.frame(x=filter(dat, MotherInf>0)$POS, y=-0.05, yend=-0.01) %>% mutate(xend=x) %>% mutate(Parent="Maternal", Child=child)
  fip <- data.frame(x=filter(dat, FatherInf>0)$POS, y=-0.05, yend=-0.01) %>% mutate(xend=x) %>% mutate(Parent="Paternal", Child=child)
  informatives <- bind_rows(mip, fip)
  setnames(dat, paste0(child, c("SameMat", "SamePat")), c("Maternal", "Paternal"))
  dat <- dat %>% 
    select(CHROM, POS, Maternal, Paternal) %>%
    gather(Parent, Marginal, -CHROM, -POS) %>%
    mutate(Child=child) 
  list(informatives=informatives, dat=dat)
}