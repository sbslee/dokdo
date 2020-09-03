library(dada2)

args <- commandArgs(trailingOnly=TRUE)
fastq_dir <- args[[1]]
output_dir <- args[[2]]
tirmF <- as.integer(args[[3]])
trimR <- as.integer(args[[4]])

dir.create(output_dir)

fnFs <- sort(list.files(fastq_dir, pattern="_R1_001.fastq", full.names=TRUE))
fnRs <- sort(list.files(fastq_dir, pattern="_R2_001.fastq", full.names=TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

filtFs <- file.path(output_dir, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(output_dir, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft=c(tirmF,trimR),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

seqtab <- makeSequenceTable(mergers)

table(nchar(getSequences(seqtab)))

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus",
                                    multithread=TRUE, verbose=TRUE)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN),
               sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF",
                     "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names

write.csv(track, paste0(output_dir, "/track.csv"))
write.csv(seqtab.nochim, paste0(output_dir, "/asv.csv"))

pdf(paste0(output_dir, "/qc.pdf"))
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])
plotErrors(errF, nominalQ=TRUE)
barplot(table(nchar(getSequences(seqtab))))
dev.off()