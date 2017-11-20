options(stringsAsFactors=F)
options(scipen=30)

read.depth.dir <- "/Users/db2175/NonDData/Data/CarmiBAMs/PGDlympho_data/PGDlympho_data/CoverageFiles/"
sequence.stats.dir <- "/Users/db2175/NonDData/Data/CarmiBAMs/PGDlympho_data/PGDlympho_data/SequenceStatFiles/"


d <- "/Users/db2175/NonDData/Data/CarmiBAMs/PGDlympho_data/PGDlympho_data/PGSblind_BAMs/"

runs.dirs <- data.frame(
  RunName=c("R1", "R2", "R1", "R2"), 
  FileList=c(paste0(d, "references_run1/samples.txt"), 
        paste0(d, "references_run2/samples.txt"), 
        paste0(d, "samples_run1/samples.txt"), 
        paste0(d, "samples_run2/samples.txt")), 
  RunCode=c("PGSRef1", "PGSRef2", "PGSSamples1", "PGSSamples2"))



GetGCCorrect <- function(){
  f <- "/Users/db2175/NonDData/Data/CarmiBAMs/PGDlympho_data/PGDlympho_data/AllReadDepth1000000.GC.txt"
  if (!file.exists(f)){
    norm <- GCCorrect()
    save(norm, file=f)
  } else {
    load(f)
  }
  norm
}

# GCCorrect
# normalizes so all samples have the same GC
# curve as the presumed normal samples
# and the same total readcount
# also adds in MeanNormals, SDNormals and MedianNormals

# chr, start, end, PctGC, Sample, ReadCount, Family, Normal, NormGC, 
# MeanNormals, SDNormals, MedianNormals
GCCorrect <- function(){
  rd <- GetAllReadDepth1000000() %>% filter(!Failed) %>% select(-Failed)
  l <- list()
  for (sample in unique(rd$Sample)){
    sub <- filter(rd, Sample==sample)
    sub$Predicted <- predict(gam(ReadCount ~ s(PctGC), family=Gamma(link="log"), data=sub), newdata=data.frame(PctGC=sub$PctGC), type='response')
    l[[sample]] <- sub
  }
  rd <- bind_rows(l)
  norm <- group_by(rd, chr, start, end) %>% 
    mutate(NormGC=ReadCount * mean(Predicted[Normal]) / Predicted) %>%
    ungroup() %>% select(-Predicted)
#   stats <- group_by(norm, chr, start, end) %>%
#     filter(Normal) %>%
#     summarize(MeanNormals=mean(NormGC[Normal]), 
#               SDNormals=sd(NormGC[Normal]), 
#               MedianNormals=median(NormGC[Normal]))
# norm <- merge(norm, stats)
  norm
}

GetAllReadDepth1000000 <- function(){
  fread("/Users/db2175/NonDData/Data/CarmiBAMs/PGDlympho_data/PGDlympho_data/AllReadDepth1000000.txt") %>% as.data.frame()
}

AllReadDepth1000000 <- function(){
  rd <- ReferenceSamplesReadDepth()
  cd <- CombineReadDepth()
  all <- bind_rows(rd, cd)
  write.table(all, "/Users/db2175/NonDData/Data/CarmiBAMs/PGDlympho_data/PGDlympho_data/AllReadDepth1000000.txt", quote=T, sep="\t", row.names=F)
}

# chr start end PctGC Sample ReadCount Family Normal Failed
CombineReadDepth <- function(){
  rd <- list()
  for (chr in 1:22){
    r <- read.table(paste0(read.depth.dir, "/chr", chr, ".txt"))
    colnames(r) <- c("chr", "start", "end", "33866", 
                     "33876", "33882", "33889")
    colnames(r) <- c("chr", "start", "end", paste0("Child", 1:4))
    s <- read.table(paste0(sequence.stats.dir, "/chr", chr, ".stats.txt"), comment.char="", header=T)
    s <- s[, c(1:3, 5, 10)]
    colnames(s) <- c("chr", "start", "end", "PctGC", "NumN")
    m <- merge(r, s)
    m$start <- m$start - 1
    rd[[chr]] <- m
  }
  rd <- bind_rows(rd) %>% 
    filter(NumN==0) %>% 
    select(-NumN) %>% 
    gather(Sample, ReadCount, -chr, -start, -end, -PctGC)
  rd$Family <- "Eth"
  rd$Normal <- rd$Sample %in% c("Child1", "Child2", "Child4")
  rd$Failed <- F
  rd
}

# chr start end PctGC Sample ReadCount Family Normal Failed
ReferenceSamplesReadDepth <- function(){
  rd <- list()
  for (i in 1:nrow(runs.dirs)){
    sample.names <- read.table(runs.dirs$FileList[i], header=F)$V1
    c.sample.names <- paste0(runs.dirs$RunName[i], "_", sample.names)
    c.sample.names <- gsub("-", "_", c.sample.names)
    sample.dict <- sample.names
    names(sample.dict) <- c.sample.names
    sg <- runs.dirs$RunCode[i]
    for (chr in 1:22){
      r <- read.table(paste0(read.depth.dir, "/", sg, "chr", 
                             chr, ".1000000.txt"))
      colnames(r) <- c("chr", "start", "end", c.sample.names)
      r$chr <- gsub("chr", "", r$chr) %>% as.numeric()
      s <- read.table(
        paste0(sequence.stats.dir, "/chr", chr, ".stats.txt"), 
        comment.char="", header=T)
      s <- s[, c(1:3, 5, 10)]
      colnames(s) <- c("chr", "start", "end", "PctGC", "NumN")
      s$start <- s$start - 1
      m <- merge(r, s) %>% gather(Sample, ReadCount, -chr, -start, -end, -PctGC, -NumN)
      m$OrigSampleName <- sample.dict[m$Sample]
      rd[[paste0(sg, chr)]] <- m
    }
  }
  all <- bind_rows(rd) %>% 
    filter(NumN==0) %>% select(-NumN)
  csv1 <- fread("/Users/db2175/NonDData/Data/CarmiBAMs/PGDlympho_data/PGDlympho_data/PGSblind_BAMs/transloseq-blindrun1.csv", skip="Sample_ID", header=T)
  csv2 <- fread("/Users/db2175/NonDData/Data/CarmiBAMs/PGDlympho_data/PGDlympho_data/PGSblind_BAMs/pipeline_sample_sheet_PGS_blindRun2.csv", skip="Sample_ID", header=T)
  csv <- rbind(csv1, csv2)[!Sample_ID=="", .(Sample_ID, Sample_Plate)] %>%
    unique(by=NULL)
  setnames(csv, c("OrigSampleName", "Family"))
  all <- merge(all, csv) %>% select(-OrigSampleName)
  presumed_normal <- c("R2_106L_1", "R2_108L_5", "R2_109L_9", 
                       "R2_110L_10", 
                       "R2_118L_12", "R2_120L_17", "R2_122L_22", 
                       "R2_125L_26", 
                       "R2_2M_3", "R1_55L_6", "R1_56L_11", 
                       "R1_57L_13", "R1_58L_15", 
                       "R1_68L_11", "R1_76L_4", 
                       "R1_XX_135F", "R2_XX_135F", 
                       "R1_XY_98F", "R2_XY_98F")
  failed <- c("R1_46L_3blk", 
    "R1_75L_3", 
    "R1_71L_1", 
    "R2_126L_27", 
    "R1_79L_7")
  all$Normal <- all$Sample %in% presumed_normal
  all$Failed <- all$Sample %in% failed
  all
}

