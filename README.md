# RandomSequenceGeneratoR
This repository describes a R workflow for the generation of random DNA sequences, based on the approach described in Tourlousse et al. (2017).

**Approach:**
The process for generating the sequences starts from a set of random k-mers that are subsequently concatenated into longer sequences and screened based on user-defined conditions, namely homopolymers, GC content, repeats, palindromes and between-sequence identity.

## Loading required R packages

```
library(dplyr)
library(Biostrings)
library(stringr)
library(stringi)
library(seqinr)
```

## Generating set of k-mers

A set of k-mers is generated using R's `sample` command, using the following set of parameters:

`kmer_seed` = 123 *# seed for reproducibility #*

`kmer_size` = 12 *# k-mer size #*

`kmer_count` = 10000 *# number of k-mers to generate #*

`kmer_baseProbabilitiesACGT` = c(0.25,0.25,0.25,0.25) *# sampling probabilities #*

`kmer_GClow` = 30 *# lower bound on k-mer GC content #*

`kmer_GChigh` = 70 *# upper bound on k-mer GC content #*

`kmer_maximumLengthHompolymers` = 4 *# upper bound on homopolymer length #*

Based on this set of parameters, additional script parameters are set as follows:

```
set.seed(kmer_seed)
seedNumbers <- sample(1:(kmer_count*10000), kmer_count*100, replace = FALSE)
outputCounter <- 1
seedCounter <- 1
```

The k-mers are then generated using the following loop:

```
kmersDF <- NULL
while (outputCounter <= kmer_count) {
  set.seed(seedNumbers[seedCounter])
  tmp_kmersDF <- data.frame("kmerID" = paste0("kmer", seedCounter),
                            "kmer" = paste(Biostrings::DNA_BASES[sample(1:4,
                                                                        kmer_size,
                                                                        rep = TRUE,
                                                                        prob = kmer_baseProbabilitiesACGT)],
                                           collapse=""))
  tmp_kmersDF <- tmp_kmersDF %>% 
    dplyr::filter(!grepl(paste(rep("A", kmer_maximumLengthHompolymers), collapse=""), kmer)) %>% 
    dplyr::filter(!grepl(paste(rep("C", kmer_maximumLengthHompolymers), collapse=""), kmer)) %>% 
    dplyr::filter(!grepl(paste(rep("T", kmer_maximumLengthHompolymers), collapse=""), kmer)) %>% 
    dplyr::filter(!grepl(paste(rep("G", kmer_maximumLengthHompolymers), collapse=""), kmer)) %>%
    dplyr::mutate(As = str_count(kmer, "A")) %>%
    dplyr::mutate(Cs = str_count(kmer, "C")) %>%
    dplyr::mutate(Gs = str_count(kmer, "G")) %>%
    dplyr::mutate(Ts = str_count(kmer, "T")) %>%
    dplyr::filter(As > 0 & Cs > 0 & Gs > 0 & Ts > 0) %>%
    dplyr::mutate(GCcontent = round((Cs+Gs)/kmer_size*100, 1)) %>%
    dplyr::filter(GCcontent >= kmer_GClow & GCcontent <= kmer_GChigh) %>%
    dplyr::mutate_if(is.factor, as.character)
  kmersDF <- rbind(kmersDF, tmp_kmersDF) %>% dplyr::distinct(kmer, .keep_all = TRUE) 
  outputCounter = 1 + nrow(kmersDF)
  seedCounter = seedCounter + 1
}
rm(tmp_kmersDF)
```
