#Sjogren-s-syndrome-Microbiome-Analysis

#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 2) {
  stop("Uso: Rscript dada2_pipeline.R <ruta_fastq> <ruta_referencia_fasta>")
}

fastq_path <- normalizePath(args[1])
ref_fasta <- normalizePath(args[2])

start_time <- Sys.time()
set.seed(531)

library(dada2)

# List fastq files
fastqs <- list.files(fastq_path, pattern = "\\.fastq\\.gz$", full.names = TRUE)
fnFs <- sort(fastqs[grepl("_1\\.fastq\\.gz$", fastqs)])
fnRs <- sort(fastqs[grepl("_2\\.fastq\\.gz$", fastqs)])

# Sample names
sample.names <- sapply(strsplit(basename(fnFs), "_1.fastq.gz"), `[`, 1)

# Filtered output paths (Quality)
filt_path <- file.path(fastq_path, "filtered")
dir.create(filt_path, showWarnings = FALSE)
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

# Verificar y sincronizar el número de lecturas forward y reverse
for (i in seq_along(fnFs)) {
  fq1 <- ShortRead::readFastq(fnFs[i])
  fq2 <- ShortRead::readFastq(fnRs[i])
  
  # Verificar si el número de lecturas difiere
  n <- min(length(fq1), length(fq2))
  
  if (length(fq1) != length(fq2)) {
    message(paste("Lecturas desajustadas en la muestra:", sample.names[i], "- recortando a", n, "lecturas."))
    fq1 <- fq1[1:n]
    fq2 <- fq2[1:n]
    
    # Guardar las lecturas recortadas en nuevos archivos fastq
    ShortRead::writeFastq(fq1, file.path(filt_path, paste0(sample.names[i], '_F_trunc.fastq.gz')))
    ShortRead::writeFastq(fq2, file.path(filt_path, paste0(sample.names[i], '_R_trunc.fastq.gz')))
    
    # Actualizar las rutas de los archivos a los recortados
    fnFs[i] <- file.path(filt_path, paste0(sample.names[i], '_F_trunc.fastq.gz'))
    fnRs[i] <- file.path(filt_path, paste0(sample.names[i], '_R_trunc.fastq.gz'))
  }
}

# Filtering
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN = 0, maxEE = c(5, 5), truncQ = 2, rm.phix = TRUE, compress = TRUE, multithread = TRUE)
DFout <- data.frame(out)

# Keep only samples with enough reads
sample.names0 <- sapply(strsplit(rownames(subset(DFout, reads.out > 50)), "_1.fastq.gz"), `[`, 1)
filtFs0 <- file.path(filt_path, paste0(sample.names0, "_F_filt.fastq.gz"))
filtRs0 <- file.path(filt_path, paste0(sample.names0, "_R_filt.fastq.gz"))

# Error rates
errF <- learnErrors(filtFs0, nbases = 1e5, multithread = TRUE)
errR <- learnErrors(filtRs0, nbases = 1e5, multithread = TRUE)

# Dereplication
derepFs <- derepFastq(filtFs0)
derepRs <- derepFastq(filtRs0)
names(derepFs) <- sample.names0
names(derepRs) <- sample.names0

# Denoising
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)

# Merge pairs
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)

# Sequence table and chimera removal
seqtab <- makeSequenceTable(mergers)
seqtabHS <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE)
saveRDS(seqtabHS, "dd2.seqtab.rds")

# Remove chimeras (redundant step, kept for compatibility)
seqtab_all_no_chimeras <- removeBimeraDenovo(seqtabHS, method = "consensus", multithread = TRUE)
saveRDS(seqtab_all_no_chimeras, "dd2.seqtab_all_no_chimeras.rds")

# Taxonomy assignment
taxa <- assignTaxonomy(seqtab_all_no_chimeras, refFasta = ref_fasta, multithread = TRUE)
saveRDS(taxa, "ASV_seq_rdp_set18.rds")

# Tracking
getN <- function(x) sum(getUniques(x))
track <- cbind(out,
               sapply(dadaFs, getN),
               sapply(dadaRs, getN),
               sapply(mergers, getN),
               rowSums(seqtabHS))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names

# Save results
write.csv(taxa, "taxTab.csv")
write.csv(seqtab_all_no_chimeras, "seqtabNoC.csv")
write.csv(track, "track.csv")

# Time
cat("Pipeline completed in", difftime(Sys.time(), start_time, units = "mins"), "minutes\n")
