#!/usr/bin/env Rscript

.libPaths(c("./.Rlibs", .libPaths()))

options(timeout = 600)

library(cancereffectsizeR)
library(ces.refset.hg38)
library(Biostrings)
library(stringr)
library(data.table)

maf <- preload_maf(maf="../data/tcga.luad.maf.txt",
                   refset=ces.refset.hg38,
                   keep_extra_columns = T)

cesa <- CESAnalysis(refset="ces.refset.hg38")
cesa <- load_maf(cesa=cesa, maf=maf)

signatures_to_remove <- suggest_cosmic_signatures_to_remove(cancer_type="LUAD",
                                                            treatment_naive=T)

cesa <- trinuc_mutation_rates(
    cesa,
    signature_set=ces.refset.hg38$signatures$COSMIC_v3.2,
    signatures_to_remove=signatures_to_remove)

cesa <- gene_mutation_rates(cesa, covariates = "lung", save_all_dndscv_output=T)

## Not really necessary for the epistasis project, but we might want
## to compare the selection coefficients without the epistasis
## consideration
cesa <- ces_variant(cesa=cesa, run_name = "general")

## Extract computed mutation rates

gene_list = c('KRAS','TP53','LRP1B','BRAF','STK11')

## There are multiple isoforms for KRAS and TP53, we choose the one
## with lowest q value (at least for TP53)
dndscv_results <- cesa$dNdScv_results[[1]][[2]]
## Take lowest q-value record for each gene, and sort by significance.
sig_genes <- dndscv_results[, .SD[which.min(qallsubs_cv)],
                            by = "gene"][order(qallsubs_cv)]

gene_pids <- sig_genes[, .(gene, pid)]

## However for KRAS, we need another one namely PID: "ENSP00000256078.5"
gene_pids[gene_pids$gene == "KRAS"]$pid = "ENSP00000256078.5"

cesa$gene_rates[cesa$gene_rates$pid %in% gene_pids[gene %in% gene_list]$pid]

## Methods to obtain the mutation rate per variant

## this order is hard-coded here just because the
## cancereffectsizeR:::.ces_ref_data object is not named. But it
## can be obtained from
## colnames(cesa@trinucleotide_mutation_weights$trinuc_proportion_matrix)
## for any cesa object
order_of_trinuc_muts <- c("A[C>A]A", "A[C>A]C", "A[C>A]G", "A[C>A]T",
                          "C[C>A]A", "C[C>A]C", "C[C>A]G", "C[C>A]T",
                          "G[C>A]A", "G[C>A]C", "G[C>A]G", "G[C>A]T",
                          "T[C>A]A", "T[C>A]C", "T[C>A]G", "T[C>A]T",
                          "A[C>G]A", "A[C>G]C", "A[C>G]G", "A[C>G]T",
                          "C[C>G]A", "C[C>G]C", "C[C>G]G", "C[C>G]T",
                          "G[C>G]A", "G[C>G]C", "G[C>G]G", "G[C>G]T",
                          "T[C>G]A", "T[C>G]C", "T[C>G]G", "T[C>G]T",
                          "A[C>T]A", "A[C>T]C", "A[C>T]G", "A[C>T]T",
                          "C[C>T]A", "C[C>T]C", "C[C>T]G", "C[C>T]T",
                          "G[C>T]A", "G[C>T]C", "G[C>T]G", "G[C>T]T",
                          "T[C>T]A", "T[C>T]C", "T[C>T]G", "T[C>T]T",
                          "A[T>A]A", "A[T>A]C", "A[T>A]G", "A[T>A]T",
                          "C[T>A]A", "C[T>A]C", "C[T>A]G", "C[T>A]T",
                          "G[T>A]A", "G[T>A]C", "G[T>A]G", "G[T>A]T",
                          "T[T>A]A", "T[T>A]C", "T[T>A]G", "T[T>A]T",
                          "A[T>C]A", "A[T>C]C", "A[T>C]G", "A[T>C]T",
                          "C[T>C]A", "C[T>C]C", "C[T>C]G", "C[T>C]T",
                          "G[T>C]A", "G[T>C]C", "G[T>C]G", "G[T>C]T",
                          "T[T>C]A", "T[T>C]C", "T[T>C]G", "T[T>C]T",
                          "A[T>G]A", "A[T>G]C", "A[T>G]G", "A[T>G]T",
                          "C[T>G]A", "C[T>G]C", "C[T>G]G", "C[T>G]T",
                          "G[T>G]A", "G[T>G]C", "G[T>G]G", "G[T>G]T",
                          "T[T>G]A", "T[T>G]C", "T[T>G]G", "T[T>G]T")


compute_genome_trinucleotide_mut_proportions <- function(cesa){
    variants = c(unlist(cesa$variants[variant_type == 'aac',
                                      constituent_snvs]),
                 cesa$variants[variant_type == 'snv',variant_id])

    temp1 = str_split_fixed(variants, pattern = ':', n=2)
    temp2 = str_split_fixed(temp1[,2], pattern= '_', n=2)
    variants = as.data.table(cbind(temp1[,1],temp2))

    colnames(variants) = c('chr','pos','mut')
    variants$pos = strtoi(variants$pos)

    ## One sample did not had a position (strangely)
    variants = variants[!is.na(variants$pos)]

    variants$trinuc = as.vector(as.character(getSeq(
        BSgenome.Hsapiens.UCSC.hg38,
        names=str_c('chr', variants$chr),
        start=variants$pos-1,
        width=3,
        strand = '+')))

    ## reverse complementing variants as necessary
    middle = str_c(substr(variants$trinuc,2,2),collapse="")
    revcomp_ind = c(str_locate_all(middle,'G|A')[[1]][,1])
    variants[revcomp_ind,
             ':=' (trinuc=as.character(reverseComplement(DNAStringSet(
                       variants[revcomp_ind, trinuc]))),
                   mut=paste0(as.character(reverseComplement(DNAStringSet(
                       variants[revcomp_ind, substr(mut,1,1)]))),
                       '>',
                       as.character(reverseComplement(DNAStringSet(
                           variants[revcomp_ind, substr(mut,3,3)])))))]

    ## generating COSMIC type format for trinuc mutations
    variants$trinuc_mut = paste0(substr(variants$trinuc, 1, 1),
                                 '[', variants$mut, ']',
                                 substr(variants$trinuc, 3, 3))

    ## concatenated table contains the names of trinuc_muts not
    ## represented in the dataset (subtract 1 because table is
    ## auto-iniated with 1s) divides number of observations by number of
    ## variants to get proportions
    trinuc_mut_proportions = (
        c(table(variants$trinuc_mut),
          table(setdiff(order_of_trinuc_muts, variants$trinuc_mut)) - 1)
        / nrow(variants))

  ## TODO: order table before returning value
  return(trinuc_mut_proportions)
}

## it takes a while to compute, so we store it as a global variable
trinucleotide_mut_proportions <- compute_genome_trinucleotide_mut_proportions(cesa)


get_gene_pid <- function(the_gene){
    pid = gene_pids[gene == the_gene]$pid

    return(pid)
}


compute_trinucleotide_contexts <- function(gene){

    pid = get_gene_pid(gene)

    coding_seq = DNAString(paste0(
        as.character(cancereffectsizeR:::.ces_ref_data$
                     ces.refset.hg38$RefCDS[[pid]]$seq_cds1up[1:2]),
        as.character(cancereffectsizeR:::.ces_ref_data$
                     ces.refset.hg38$RefCDS[[pid]]$seq_cds1down)))

    tri_nt_contexts = data.frame(as.list(Biostrings::trinucleotideFrequency(coding_seq)))

    cols_to_remove = c()

    for(tri in colnames(tri_nt_contexts)){
        if(substr(tri,2,2) %in% c('G','A')){
            cols_to_remove = c(cols_to_remove, tri)
            tri_nt_contexts[as.character(
                reverseComplement(DNAString(tri)))] =
                tri_nt_contexts[as.character(
                    reverseComplement(DNAString(tri)))] +
                tri_nt_contexts[[tri]]
        }
    }

    tri_nt_contexts[, cols_to_remove] <- list(NULL)

    return(tri_nt_contexts)
}


unravel_constituent_snv <- function(constituent_snv){
    first_split <- str_split_fixed(constituent_snv, pattern = ':', n=2)

    second_split <- str_split_fixed(first_split[[2]], pattern = '_', n=2)

    unraveled <- as.data.table(cbind(first_split[[1]], second_split))
    colnames(unraveled) <- c('chr','pos','mut')
    unraveled$pos = strtoi(unraveled$pos)
    return(unraveled)
}


get_context <- function(variant){

    variant_info <- unravel_constituent_snv(variant)

    context <-as.character(getSeq(BSgenome.Hsapiens.UCSC.hg38,
                                  names=str_c('chr', variant_info$chr),
                                  start=variant_info$pos-1,
                                  width=3, strand = '-'))

    ## reverse complement if necessary
    if (substr(context, 2, 2) %in% c("G", "A")){
        context <- as.character(reverseComplement(DNAStringSet(context)))
    }

    return(context)
}


get_trinuc_mut <- function(variant){
    ## generate COSMIC type format for trinuc mutations
    mut <- unravel_constituent_snv(variant)$mut
    context <- get_context(variant)

    ## HAVE TO FIX MUT IF THE VARIANT HAS A OR G
    if(substr(mut,1,1) %in% c('G','A')){
        mut <- paste0(as.character(reverseComplement(DNAStringSet(substr(mut,1,1)))),
                      '>',
                      as.character(reverseComplement(DNAStringSet(substr(mut,3,3)))))
    }

    trinuc_mut <- paste0(substr(context,1,1),
                         '[', mut, ']',
                         substr(context,3,3))
    return(trinuc_mut)
}


## get_freq_trinuc_mut_in_gene <- function(trinuc_mut, gene){

##     index_context <- match(trinuc_mut, order_of_trinuc_muts)

##     gene_pid <- get_gene_pid(gene)

##     ## TODO: This line needs to be fixed!
##     frequencies <- cancereffectsizeR:::.ces_ref_data$ces.refset.hg38[["gene_trinuc_comp"]][[gene_pid]]

##     return(frequencies[index_context])
## }


variant_mutation_rate <- function(variant, gene, cesa, trinuc_mut_proportions){
    gene_pid <- get_gene_pid(gene)
    context <- get_context(variant)
    trinuc_mut <- get_trinuc_mut(variant)

    gene_mut_rate <- cesa$gene_rates[pid == gene_pid]$rate

    freq_trinuc_mut <- trinuc_mut_proportions[[trinuc_mut]]

    times_context_in_gene <- compute_trinucleotide_contexts(gene)[[context]]

    return(gene_mut_rate*freq_trinuc_mut/times_context_in_gene)
}



all_variants_in_gene <- function(the_gene, cesa){
    ## unlist because some codons have multiple variants

    variants =  c(
        unlist(cesa$variants[gene == the_gene & variant_type == 'aac',
                             constituent_snvs]),
        cesa$variants[gene == the_gene & variant_type == 'snv',
                      variant_id])
    return(variants)
}

new_all_variants_in_gene <- function(the_gene, cesa){
    ## unlist because some codons have multiple variants

    variants =  unique(cesa$maf[top_gene == the_gene & !is.na(variant_id),
                                variant_id])

    return(variants)
}



## Mutation rates but only for the variants included (all of them pathogenic for KRAS)
sapply(gene_list,
       function(y){
           sum(sapply(new_all_variants_in_gene(y, cesa),
                      function(x){
                          variant_mutation_rate(x,
                                                y,
                                                cesa,
                                                trinucleotide_mut_proportions)}))})

## Results:
##     KRAS         TP53        LRP1B         BRAF        STK11
## 2.030913e-08 2.262881e-07 2.295583e-07 8.843950e-09 4.544432e-08
##     KRAS         TP53        LRP1B         BRAF        STK11
## 1.992647e-08 1.808964e-07 1.559768e-07 6.277220e-09 2.543602e-08
