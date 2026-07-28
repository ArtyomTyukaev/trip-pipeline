library(rtracklayer)
library(dplyr)

names_with_P_list = c('5ss_Expression_rep1_S56', '5ss_Expression_rep2_S57',
                      'CPA_Expression_rep2_S54', 'CPA_Expression_rep3_S55')

# Path to dir with TRIP data (*.combined_counts.txt)
path = "/home/projects/ulitsky/lidiag/TRIP_data/set1_combined_1/"

# Path to file with exp names
namesExp<-read.table(paste0(path,"names_all_clean.txt"))$V1

# Arguments from the terminal
myArgs = commandArgs(trailingOnly = T)
if (length(myArgs)>0)
{
  path = myArgs[1]
  namesExp_file = myArgs[2]
  namesExp<-read.table(paste0(path, namesExp_file))$V1
}
#path="/home/projects/ulitsky/shared/data/human/TRIP/IntronTRIP_run3/"
#namesExp = read.table(paste0(path, "names.exp.txt"))$V1
print(path)
print(namesExp)
# External annotation data
trans2gene = read.table("/home/projects/ulitsky/igoru/data/human/MANE.GRCh38.v1.0.refseq_genomic.trans2gene.txt",sep="\t")
maneExp = read.table("/home/projects/ulitsky/igoru/data/human/ENCODE/K562/K562_1_MANEReadsPerGene.out.tab",sep="\t",header=F)
maneBED = read.table("/home/projects/ulitsky/igoru/data/human/MANE.GRCh38.v1.0.refseq_genomic.bed",sep="\t",header=F)
maneBED$Length = maneBED$V3-maneBED$V2
maneBED$ExonicLength = sapply(strsplit(as.character(maneBED$V11),","),function(x){sum(as.numeric(x))})
intronBED = read.table("/home/projects/ulitsky/lidiag/genomic_data/GCF_000001405.40_GRCh38.p14_MANE_introns_new_format.bed",sep="\t",header=F)
intronBED$Length = intronBED$V3-intronBED$V2
gencodeBED = read.table("/home/projects/ulitsky/lidiag/genomic_data/gencode.v46.annotation.noncoding.bed",sep="\t",header=F)
gencodeBED$Length = gencodeBED$V3-gencodeBED$V2

bed_chrom_hmm <- gzfile("/home/projects/ulitsky/shared/data/human/K562_ChromHMM_Kellis_ENCFF319VXX.bed.gz")
chrom_hmm = import(bed_chrom_hmm)

AnnotTabs = list()
repeats = list()
GcTabs = list()
for (name in namesExp)
{
  print(name)
  # Reading data
  AnnotTabs[[name]] = read.table(paste0(path, name, ".combined_counts.txt"), sep = "\t", header = T)
  table(AnnotTabs[[name]]$exp)
  
  # Adding expression
  AnnotTabs[[name]]$exp.norm = log((1 + AnnotTabs[[name]]$exp) / (1 + AnnotTabs[[name]]$norm), 2)
  AnnotTabs[[name]]$exp.log10 = log((1 + AnnotTabs[[name]]$exp), 10)
  
  # Adding annotation (to obtain tracks use make_tracks.sh)
  maneSense = read.table(paste0(path, name, ".MANE_all.txt"), sep = "\t", header = F,quote="",comment="")
  maneSense[941:943,]
  maneExonic = read.table(paste0(path, name, ".MANE_exonic.txt"), sep = "\t", header = F)
  maneAntisense = read.table(paste0(path, name, ".MANE_all_antisense.txt"), sep = "\t", header = F)
  maneExonicAntisense = read.table(paste0(path, name, ".MANE_exonic_antisense.txt"), sep = "\t", header = F)
  intronSense = read.table(paste0(path, name, ".MANE_introns_sense.txt"), sep = "\t", header = F)
  intronAntisense = read.table(paste0(path, name, ".MANE_introns_antisense.txt"), sep = "\t", header = F)
  gencode_noncoding = read.table(paste0(path, name, ".gencode_noncoding_all.txt"), sep = "\t", header = F)
  
  # Add MANE data to AnnotTabs
  AnnotTabs[[name]]$MANE_all = maneSense[match(rownames(AnnotTabs[[name]]), maneSense$V4), "V10"]
  AnnotTabs[[name]]$MANE_exonic = maneExonic[match(rownames(AnnotTabs[[name]]), maneExonic$V4), "V10"]
  AnnotTabs[[name]]$MANE_antisense = maneAntisense[match(rownames(AnnotTabs[[name]]), maneAntisense$V4), "V10"]
  AnnotTabs[[name]]$MANE_exonic_antisense = maneExonicAntisense[match(rownames(AnnotTabs[[name]]), maneExonicAntisense$V4), "V10"]
  AnnotTabs[[name]]$MANE_intron = intronSense[match(rownames(AnnotTabs[[name]]), intronSense$V4), "V10"]
  AnnotTabs[[name]]$MANE_intron_antisense = intronAntisense[match(rownames(AnnotTabs[[name]]), intronAntisense$V4), "V10"]
  AnnotTabs[[name]]$gencode_noncoding = gencode_noncoding[match(rownames(AnnotTabs[[name]]), gencode_noncoding$V4), "V10"]
  AnnotTabs[[name]]$annot_rna_type = gencode_noncoding[match(rownames(AnnotTabs[[name]]), gencode_noncoding$V4), "V13"]
  
  AnnotTabs[[name]] <- AnnotTabs[[name]] %>%
    mutate(
      annot_rna_type = as.character(annot_rna_type),
      annot_rna_type = case_when(
        is.na(annot_rna_type) & !is.na(MANE_all) ~ "MANE_mRNA",
        TRUE ~ annot_rna_type
      )
    )
  AnnotTabs[[name]]$gencode_noncoding_name = gencode_noncoding[match(rownames(AnnotTabs[[name]]), gencode_noncoding$V4), "V14"]
  
  # Add gene symbols and expression data
  AnnotTabs[[name]]$MANE_all_symbol = trans2gene[match(AnnotTabs[[name]]$MANE_all, trans2gene$V2), "V1"]
  AnnotTabs[[name]]$MANE_all_exp = maneExp[match(trans2gene[match(AnnotTabs[[name]]$MANE_all, trans2gene$V2), "V1"], maneExp$V1), "V2"]
  AnnotTabs[[name]]$MANE_exonic_exp = maneExp[match(trans2gene[match(AnnotTabs[[name]]$MANE_exonic, trans2gene$V2), "V1"], maneExp$V1), "V2"]
  AnnotTabs[[name]]$MANE_antisense_exp = maneExp[match(trans2gene[match(AnnotTabs[[name]]$MANE_antisense, trans2gene$V2), "V1"], maneExp$V1), "V2"]
  AnnotTabs[[name]]$MANE_exonic_antisense_exp = maneExp[match(trans2gene[match(AnnotTabs[[name]]$MANE_exonic_antisense, trans2gene$V2), "V1"], maneExp$V1), "V2"]
  
  # Normalize expression data
  AnnotTabs[[name]]$MANE_all_exp.norm = 1000 * (AnnotTabs[[name]]$MANE_all_exp / maneBED[match(AnnotTabs[[name]]$MANE_all, maneBED$V4), "ExonicLength"])
  AnnotTabs[[name]]$MANE_all_exp.norm.log10 <- log((AnnotTabs[[name]]$MANE_all_exp.norm + 1), 10)
  AnnotTabs[[name]]$MANE_exonic_exp.norm = 1000 * (AnnotTabs[[name]]$MANE_exonic_exp / maneBED[match(AnnotTabs[[name]]$MANE_exonic, maneBED$V4), "ExonicLength"])
  AnnotTabs[[name]]$MANE_antisense_exp.norm = 1000 * (AnnotTabs[[name]]$MANE_antisense_exp / maneBED[match(AnnotTabs[[name]]$MANE_antisense, maneBED$V4), "ExonicLength"])
  AnnotTabs[[name]]$MANE_exonic_antisense_exp.norm = 1000 * (AnnotTabs[[name]]$MANE_exonic_antisense_exp / maneBED[match(AnnotTabs[[name]]$MANE_exonic_antisense, maneBED$V4), "ExonicLength"])
  
  # Process intronic data
  AnnotTabs[[name]]$MANE_intronic = AnnotTabs[[name]]$MANE_all
  AnnotTabs[[name]]$MANE_intronic[!is.na(AnnotTabs[[name]]$MANE_exonic)] = NA
  
  # Calculate distances and relative positions
  AnnotTabs[[name]]$MANElength = maneBED[match(AnnotTabs[[name]]$MANE_all, maneBED$V4), "Length"]
  AnnotTabs[[name]]$MANE_intron_length = intronBED[match(AnnotTabs[[name]]$MANE_intron, intronBED$V4), "Length"]
  AnnotTabs[[name]]$DistToStart = (AnnotTabs[[name]]$Pos - maneBED[match(AnnotTabs[[name]]$MANE_intronic, maneBED$V4), "V2"])
  AnnotTabs[[name]]$DistToEnd = (maneBED[match(AnnotTabs[[name]]$MANE_intronic, maneBED$V4), "V3"] - AnnotTabs[[name]]$Pos)
  AnnotTabs[[name]]$DistToTSS = ifelse(AnnotTabs[[name]]$Strand == "+", AnnotTabs[[name]]$DistToStart, AnnotTabs[[name]]$DistToEnd)
  AnnotTabs[[name]]$DistToTES = ifelse(AnnotTabs[[name]]$Strand == "-", AnnotTabs[[name]]$DistToStart, AnnotTabs[[name]]$DistToEnd)
  AnnotTabs[[name]]$RelPos = AnnotTabs[[name]]$DistToTSS / AnnotTabs[[name]]$MANElength
  
  
  # Calculate distances to intron boundaries
  AnnotTabs[[name]]$DistToIntronStart = (AnnotTabs[[name]]$Pos - intronBED[match(AnnotTabs[[name]]$MANE_intron, intronBED$V4), "V2"])
  AnnotTabs[[name]]$DistToIntronEnd = (intronBED[match(AnnotTabs[[name]]$MANE_intron, intronBED$V4), "V3"] - AnnotTabs[[name]]$Pos)
  AnnotTabs[[name]]$DistToIntron5 = ifelse(AnnotTabs[[name]]$Strand == "+", AnnotTabs[[name]]$DistToIntronStart, AnnotTabs[[name]]$DistToIntronEnd)
  AnnotTabs[[name]]$DistToIntron3 = ifelse(AnnotTabs[[name]]$Strand == "-", AnnotTabs[[name]]$DistToIntronStart, AnnotTabs[[name]]$DistToIntronEnd)
  AnnotTabs[[name]]$RelIntronPos = AnnotTabs[[name]]$DistToIntron5 / AnnotTabs[[name]]$MANE_intron_length
  
  # Adding non-coding genes info 
  AnnotTabs[[name]]$gencodelength = gencodeBED[match(AnnotTabs[[name]]$gencode_noncoding, gencodeBED$V4), "Length"]
  AnnotTabs[[name]]$DistToStart <- ifelse(
    is.na(AnnotTabs[[name]]$DistToStart), 
    AnnotTabs[[name]]$Pos - gencodeBED[match(AnnotTabs[[name]]$gencode_noncoding, gencodeBED$V4), "V2"],
    AnnotTabs[[name]]$DistToStart
  )
  
  AnnotTabs[[name]]$DistToEnd <- ifelse(
    is.na(AnnotTabs[[name]]$DistToEnd), 
    gencodeBED[match(AnnotTabs[[name]]$gencode_noncoding, gencodeBED$V4), "V3"] - AnnotTabs[[name]]$Pos, 
    AnnotTabs[[name]]$DistToEnd
  )
  AnnotTabs[[name]]$DistToTSS <- ifelse(
    is.na(AnnotTabs[[name]]$DistToTSS), 
    ifelse(AnnotTabs[[name]]$Strand == "+", AnnotTabs[[name]]$DistToStart, AnnotTabs[[name]]$DistToEnd), 
    AnnotTabs[[name]]$DistToTSS
  )
  
  AnnotTabs[[name]]$DistToTES <- ifelse(
    is.na(AnnotTabs[[name]]$DistToTES), 
    ifelse(AnnotTabs[[name]]$Strand == "-", AnnotTabs[[name]]$DistToStart, AnnotTabs[[name]]$DistToEnd), 
    AnnotTabs[[name]]$DistToTES
  )
  
  AnnotTabs[[name]]$RelPos <- ifelse(
    is.na(AnnotTabs[[name]]$RelPos), 
    AnnotTabs[[name]]$DistToTSS / AnnotTabs[[name]]$gencodelength, 
    AnnotTabs[[name]]$RelPos
  )
  
  # Adding ChromHMM
  myPoses = GRanges(seqnames = AnnotTabs[[name]]$Chr, ranges = IRanges(start = AnnotTabs[[name]]$Pos, width = 1))
  overlapsChromHmm = findOverlaps(myPoses, chrom_hmm)
  AnnotTabs[[name]][queryHits(overlapsChromHmm), "ChromHMM"] = chrom_hmm[subjectHits(overlapsChromHmm)]$name
  
  #Adding GC content
  GcTabs[[name]] <- read.table(paste0(path, name, "_with_gc.bed"), sep = "\t", header = FALSE,fill=T)
  AnnotTabs[[name]]$GC <- GcTabs[[name]]$V8[match(rownames(AnnotTabs[[name]]), GcTabs[[name]]$V7)]
  
  # Deviding on groups
  gc_mean <- mean(AnnotTabs[[name]]$GC, na.rm = TRUE)
  gc_sd <- sd(AnnotTabs[[name]]$GC, na.rm = TRUE)
  
  low_threshold <- gc_mean - gc_sd
  high_threshold <- gc_mean + gc_sd
  
  AnnotTabs[[name]]$GC_category <- cut(
    AnnotTabs[[name]]$GC,
    breaks = c(-Inf, low_threshold, high_threshold, Inf),
    labels = c("Low GC", "Average GC", "High GC"),
    include.lowest = TRUE
  )
  
  #Adding repeats
  repeats[[name]] <- read.table(paste0(path, name, ".repeats.txt"), sep = "\t", header = FALSE)
  AnnotTabs[[name]]$RepeatGroup <- repeats[[name]]$V13[match(rownames(AnnotTabs[[name]]), repeats[[name]]$V4)]
  
  #Adding exp class (linear regression exp ~ host.exp)
  if (name %in% names_with_P_list) {
    valid_rows <- complete.cases(AnnotTabs[[name]]$exp.norm, AnnotTabs[[name]]$MANE_all_exp.norm.log10)
    
    model <- lm(exp.norm ~ MANE_all_exp.norm.log10, data = AnnotTabs[[name]][valid_rows, ])
    
    AnnotTabs[[name]]$residuals <- NA
    AnnotTabs[[name]]$category <- NA
    
    AnnotTabs[[name]]$residuals[valid_rows] <- residuals(model)
    
    quantile_breaks <- quantile(AnnotTabs[[name]]$residuals[valid_rows], probs = seq(0, 1, length.out = 6))
    AnnotTabs[[name]]$category[valid_rows] <- cut(AnnotTabs[[name]]$residuals[valid_rows], breaks = quantile_breaks, labels = c("Class1", "Class2", "Class3", "Class4", "Class5"), include.lowest = TRUE)
  } else {
    valid_rows <- complete.cases(AnnotTabs[[name]]$exp.log10, AnnotTabs[[name]]$MANE_all_exp.norm.log10)
    max(AnnotTabs[[name]]$exp.log10)
    
    model <- lm(exp.log10 ~ MANE_all_exp.norm.log10, data = AnnotTabs[[name]][valid_rows, ])
    
    AnnotTabs[[name]]$residuals <- NA
    AnnotTabs[[name]]$category_5classes <- NA
    AnnotTabs[[name]]$category_3classes <- NA
    
    AnnotTabs[[name]]$residuals[valid_rows] <- residuals(model)
    
    # Deviding into 5 classes
    quantile_breaks <- quantile(AnnotTabs[[name]]$residuals[valid_rows], probs = seq(0, 1, length.out = 6))
    AnnotTabs[[name]]$category_5classes[valid_rows] <- cut(AnnotTabs[[name]]$residuals[valid_rows], breaks = quantile_breaks, labels = c("Class1", "Class2", "Class3", "Class4", "Class5"), include.lowest = TRUE)
    
    # Deviding into 3 classes
    quantile_breaks <- quantile(AnnotTabs[[name]]$residuals[valid_rows], probs = seq(0, 1, length.out = 4))
    AnnotTabs[[name]]$category_3classes[valid_rows] <- cut(AnnotTabs[[name]]$residuals[valid_rows], breaks = quantile_breaks, labels = c("Class1", "Class2", "Class3"), include.lowest = TRUE)
    
    AnnotTabs[[name]]$category_5classes <- factor(AnnotTabs[[name]]$category_5classes, levels = c(1, 2, 3, 4, 5),
                                                  labels = c("Class1", "Class2", "Class3", "Class4", "Class5"))
    
    AnnotTabs[[name]]$category_3classes <- factor(AnnotTabs[[name]]$category_3classes, levels = c(1, 2, 3),
                                                  labels = c("Class1", "Class2", "Class3"))
    
  }
  
  
  
  # Write final table to file
  write.table(AnnotTabs[[name]], paste0(path, name, ".annotated_tables.txt"), quote = F, col.names = T, row.names = T, sep = "\t")
}
