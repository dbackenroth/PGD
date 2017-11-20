# timing for with and without forward backward for chromosome 21

# timing for all together and separately for without forward backward

# accuracy for all together and separately

FBTiming <- function(){
  ptm <- proc.time()
  HMMRealData(chr="chr1", children=c("Child2", "Child3", "Child4"), forward.backward=F)
  print((proc.time()-ptm)[3])
  ptm <- proc.time()
  HMMRealData(chr="chr1", children=c("Child2", "Child3", "Child4"), forward.backward=T)
  print((proc.time()-ptm)[3])
}

EvaluateHMM <- function(chr, children){
  Evaluate <- function(predictions, truth){
    pred <- gather(predictions, Sample, Predicted, -CHROM, -POS)
    truth <- mutate(truth, NumHets=as.numeric(Mother=="AB") + as.numeric(Father=="AB")) %>%
      gather(Sample, Truth, -NumHets, -CHROM, -POS)
    summary <- merge(pred, truth) %>% group_by(NumHets, Sample) %>% 
      summarize(Accuracy=mean(Truth==Predicted), 
                NumSNPs=n())
  }
  hmm.res <- HMMRealData(chr=chr, children=children, forward.backward=F)
  ground.truth <- LoadGroundTruth(chr=chr)
  Evaluate(predictions=hmm.res$predictions, truth=ground.truth)
}


AllComparison <- function(){
  l <- list()
  i <- 1
  for (chr in paste0("chr", 22:1)){
    print(chr)
    ptm <- proc.time()
    l[[i]] <- EvaluateHMM(chr=chr, children="Child2") %>%
      mutate(Chr=chr, Together=F)
    time <- proc.time() - ptm
    l[[i]]$time <- time[3]
    i <- i + 1
    ptm <- proc.time()
    l[[i]] <- EvaluateHMM(chr=chr, children=c("Child2", "Child3", "Child4")) %>%
      mutate(Chr=chr, Together=T)
    time <- proc.time() - ptm
    l[[i]]$time <- time[3]
    i <- i + 1
  }
  merged <- bind_rows(l)
  save(merged, file="Results/AllChildrenComparison.Rdata")
  res <- filter(merged, Sample %in% c("Child2", "Child3", "Child4")) %>%
    group_by(NumHets, Together, Sample) %>% 
    summarize(Acc=sum(Accuracy * NumSNPs) / sum(NumSNPs), 
              Tot=sum(NumSNPs), 
              Time=sum(time))
  browser()
}
# 
# 
# OneChildComparison <- function(){
#   l <- list()
#   i <- 1
#   for (chr in paste0("chr", 22:1)){
#     l[[i]] <- EvaluateHMM(chr=chr, children="Child2", 
#                           beta.binomial=F, truncate.limit=NA) %>%
#       mutate(Chr=chr, Together=F, BB=F, Truncate=NA)
#     i <- i + 1
#     l[[i]] <- EvaluateHMM(chr=chr, children="Child2", 
#                           beta.binomial=F, truncate.limit=5) %>%
#       mutate(Chr=chr, Together=F, BB=F, Truncate=5)
#     i <- i + 1
#     l[[i]] <- EvaluateHMM(chr=chr, children="Child2", 
#                           beta.binomial=T, truncate.limit=NA) %>%
#       mutate(Chr=chr, Together=F, BB=T, Truncate=NA)
#     i <- i + 1
#     l[[i]] <- EvaluateHMM(chr=chr, children="Child2", 
#                           beta.binomial=T, truncate.limit=5) %>%
#       mutate(Chr=chr, Together=F, BB=T, Truncate=5)
#     i <- i + 1
#   }
#   merged <- bind_rows(l)
#   save(merged, file="Results/OneChildComparison.Rdata")
#   child2.onehet <- filter(merged, NumHets==1 & Sample=="Child2")
#   acc1 <- group_by(child2.onehet, BB, Truncate) %>% summarize(Acc=sum(Accuracy * NumSNPs) / sum(NumSNPs), Tot=sum(NumSNPs))
#   child2.twohet <- filter(merged, NumHets==2 & Sample=="Child2")
#   acc2 <- group_by(child2.twohet, BB, Truncate) %>% summarize(Acc=sum(Accuracy * NumSNPs) / sum(NumSNPs), Tot=sum(NumSNPs))
# }