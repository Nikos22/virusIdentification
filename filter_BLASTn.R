#!/usr/bin/Rscript

library(mltools)
library(reshape2)
library(Hmisc)
library(tidyr)
library(dplyr)
library(stringr)
library(IRanges)
library(ggplot2)

plotRanges <- function(x, xlim=x, main=deparse(substitute(x)), col="black", sep=0.5, ...)
{
  height <- 1
  if (is(xlim, "IntegerRanges"))
    xlim <- c(min(start(xlim)), max(end(xlim)))
  bins <- disjointBins(IRanges(start(x), end(x) + 1))
  plot.new()
  plot.window(xlim, c(0, max(bins)*(height + sep)))
  ybottom <- bins * (sep + height) - height
  rect(start(x)-0.5, ybottom, end(x)+0.5, ybottom + height, col=col, ...)
  title(main)
  axis(1)
}


args <- commandArgs(TRUE)

df <- read.table("more apple virome/agia_VS_allVirusesViriodsApple_blastn_quote.txt",
                 header = FALSE,
                 sep="\t",
                 as.is = TRUE,
                 stringsAsFactors = FALSE)                                                 
                 
# blastn -query all_agia.fasta -db ~/Documents/ioanna_milia/allVirusesViroids_Malus.fasta -outfmt "6 std slen qlen qcovhsp qcovs stitle" -max_target_seqs 3000 -out agia_VS_allVirusesViriodsApple_blastn.txt  -evalue 1 -num_threads 10 &
                 
colnames(df) <- c("qaccver", "saccver", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "slen", "qlen", "qcovhsp", "qcovs", "stitle")                 

################################
######### filter hits ##########
################################


df_n <- group_by(df, qaccver) %>%
  filter(as.numeric(pident) >= 40, 
         as.numeric(evalue) < 1e-5,
         as.numeric(qcovhsp)>10) %>%
  dplyr::arrange(as.numeric(evalue),
                 desc(as.numeric(pident)),
                 desc(as.numeric(length)),
                 saccver) %>%
  filter(saccver == saccver[1]) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  group_by(saccver) %>%
  summarise(cov = round(100*sum(as.data.frame(reduce(IRanges(pmin(sstart, send), pmax(sstart, send))))$width)/unique(slen),3),
            length = unique(slen),
            no_contigs = length(unique(qaccver)),
            desc = head(stitle,1),
            contigs = toString(unique(qaccver))) %>%
  dplyr::arrange(desc(cov), desc(length)) %>%
  filter(cov>1, 
         !grepl("Malus x domestica cultivar Golden Delicious*|Malus x domestica mitochondrion*", desc))


write.csv(df_n,
			file="Viral_matches_per_virus_AGIA.csv",
			row.names = FALSE)
