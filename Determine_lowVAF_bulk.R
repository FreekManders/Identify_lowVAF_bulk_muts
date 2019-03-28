# Determine_lowVAF_bulk.R
# Freek Manders 13-2-2019
# Minor modifications by A. Rosendahl Huber 14-2-2019

library(tidyverse)
library(VariantAnnotation)
library(GenomicRanges)
library(optparse)

####___________________Parse command line arguments____________________________________________####
option_list = list(
    make_option(c("--vcf"), type="character", 
                help="The vcf used as an input", metavar="character"),
    make_option(c("-o", "--out_dir"), default = getwd(), type="character", 
                help="The output path", metavar="character"),
    make_option(c("--bulk"), type="character", 
                help="The bulk sample", metavar="character"),
    make_option(c("--sample_name"), type="character", 
                help="The name of the sample", metavar="character"),
    make_option(c("--gender"), type="character", 
                help="The gender of the sample", metavar="character"),
    make_option(c("--genome"), type="character", 
                help="The genome of the sample. Tested on hg19", default = "hg19", metavar="character") 
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

vcf_fname = opt$vcf
out_dir = opt$out_dir
bulk = opt$bulk
sample_name = opt$sample_name
gender = opt$gender
genome = opt$genome

####_______________Helper functions_____________####
#Get the vaf of a vcf object
get_vaf = function(vcf, sample){
    ad = geno(vcf)$AD[,sample]
    vaf = sapply(ad, function(x) x[[2]] / sum(x))
    return(vaf)
}

#Sorts a vcf object.
sort_vcf = function(vcf){
  if (class(vcf)[1] != "CollapsedVCF"){
    warning(paste0("This function was tested on the Collapsed VCF class. It may not work on the supplied input vcf, which is a ", class(vcf[1]), " object."))
  }
  
  if (length(vcf) == 0){
    warning("vcf is empty. Returning vcf as is.")
    return(vcf)
  }
  gr = granges(vcf)
  gr$id = seq(1, length(gr))
  gr_sorted = sort(gr)
  vcf = vcf[gr_sorted$id,]
  return(vcf)
}

#Helper function which filters one row of a tibble, containing genotype information.
filter_row = function(row, nsamples, min_vaf, min_absent, min_present, min_dp = 20, min_pres_gq, mut_type){
  
  if (mut_type == "snv"){
    min_abs_gq = 10
  } else if (mut_type == "indel"){
    min_abs_gq = 99
    min_pres_gq = 10
  }
  
  absent = grep("0/0", row)
  present = grep("0/1|1/1", row)
  absent_goodgq = row[absent + nsamples] >= min_abs_gq
  present_goodgq = row[present + nsamples] >= min_pres_gq
  absent_gooddp = row[absent + 2*nsamples] >= min_dp
  absent_goodad = row[absent + 3*nsamples] <= 0
  present_gooddp = row[present + 2*nsamples] >= min_dp
  present_goodvaf = row[present + 4*nsamples] >= min_vaf
  
  good_absents = absent_goodgq & absent_gooddp & absent_goodad
  good_presents = present_goodgq & present_gooddp & present_goodvaf
  
  keep_row = sum(good_absents) >= min_absent & sum(good_presents) >= min_present
  return(keep_row)
}

#Helper function which retrieves the genotype information of a vcf object and stores it in a tibble
vcf_gt_inf = function(vcf, gt_clones_cols){
  gt = geno(vcf)$GT %>% as_tibble() %>% rename_all(., function(x) paste0(names(.), "_GT")) %>% dplyr::select(gt_clones_cols)
  gq = geno(vcf)$GQ %>% as_tibble() %>% rename_all(., function(x) paste0(names(.), "_GQ")) %>% replace(., is.na(.), 0) %>% dplyr::select(gt_clones_cols)
  dp = geno(vcf)$DP %>% as_tibble() %>% rename_all(., function(x) paste0(names(.), "_DP")) %>% replace(., is.na(.), 0) %>% dplyr::select(gt_clones_cols)
  ad = geno(vcf)$AD %>% as_tibble() %>% dplyr::select(gt_clones_cols)
  
  ad_alt = apply(ad, 1:2, function(x) x[[1]][2]) %>% as_tibble() %>% rename_all(., function(x) paste0(names(.), "_ADALT"))
  vaf = apply(ad, 1:2, function(x) x[[1]][2] / sum(x[[1]])) %>% as_tibble() %>% rename_all(., function(x) paste0(names(.), "_vaf")) %>% replace(., is.na(.), 0)
  geno_tb = bind_cols(gt, gq, dp, ad_alt, vaf)
  return(geno_tb)
}



#Function that compares clones and identifies sites, that are present and/or absent in the specified number of clones. Can be run for both snvs and indels.
find_tree_sites = function(vcf, sample_name, bulk, min_absent, min_present, trisomy, gender, mut_type){
  
  #Filter on gt
  gt = geno(vcf)$GT
  bulk_column = grep(bulk, colnames(gt))
  gt_clones = gt[, -bulk_column]
  shared_sites_f = apply(gt_clones, 1, function(x) table(x)["0/0"] >= min_absent & (sum("0/1" == x) + sum("1/1" == x)) >= min_present) %>% replace(., is.na(.), FALSE)
  vcf_shared = vcf[shared_sites_f]
  
  if (nrow(vcf_shared) < 1){
    print(paste0("No shared mutations for sample: ", sample_name, ". The function is returning empty"))
    return(0)
  }
  
  samples = samples(header(vcf_shared))
  nsamples_bulk = length(samples)
  gt_cols = seq(1, nsamples_bulk)
  bulk_gt_col = grep(bulk, samples)
  gt_clones_cols = gt_cols[!gt_cols %in% bulk_gt_col]
  nsamples = nsamples_bulk - 1
  
  chroms_vcf = rowRanges(vcf_shared) %>% seqnames() %>% as.vector()
  other_chroms = unique(chroms_vcf)
  
  #Trisomy filter
  if (trisomy != "FALSE"){
    chr_tri = chroms_vcf == trisomy
    vcf_tri = vcf_shared[chr_tri,]
    geno_tb_tri = vcf_gt_inf(vcf_tri, gt_clones_cols)
    nsites_tri = nrow(geno_tb_tri)
    keep_rows_tri = sapply(seq(1, nsites_tri), function(i) filter_row(geno_tb_tri[i,], nsamples, min_vaf = 0.2, min_absent, min_present, min_dp = 20, min_pres_gq = 99, mut_type))
    vcf_tri = vcf_tri[keep_rows_tri]
    other_chroms = other_chroms[!other_chroms == trisomy]
  } else{
    vcf_tri = vcf_shared[0]
  }
  
  #X chromosome filter for men
  if (gender == "M"){
    chr_x = chroms_vcf == "X"
    vcf_x = vcf_shared[chr_x,]
    geno_tb_x = vcf_gt_inf(vcf_x, gt_clones_cols)
    nsites_x = nrow(geno_tb_x)
    keep_rows_x = sapply(seq(nsites_x), function(i) filter_row(geno_tb_x[i,], nsamples, min_vaf = 0.99, min_absent, min_present, min_dp = 10, min_pres_gq = 10, mut_type))
    vcf_x = vcf_x[keep_rows_x]
    other_chroms = other_chroms[!other_chroms == "X"]
  } else{
    vcf_x = vcf_shared[0]
  }
  
  #Filter for the other chromosomes
  chr_rest = chroms_vcf %in% other_chroms
  vcf_rest = vcf_shared[chr_rest,]
  geno_tb_rest = vcf_gt_inf(vcf_rest, gt_clones_cols)
  nsites_rest = nrow(geno_tb_rest)
  keep_rows_and_good_samples_rest = lapply(seq(nsites_rest), function(i) filter_row(geno_tb_rest[i,], nsamples, min_vaf = 0.3, min_absent, min_present, min_dp = 20, min_pres_gq = 99, mut_type))
  keep_rows_rest = lapply(keep_rows_and_good_samples_rest, function(x) x[[1]]) %>% unlist()
  vcf_rest = vcf_rest[keep_rows_rest]
  
  #Combine the output from the different filters
  vcf_shared_filtered = rbind(vcf_rest, vcf_tri, vcf_x)
  vcf_shared_filtered = sort_vcf(vcf_shared_filtered)
  
  return(vcf_shared_filtered)
}

####______________Perform actual filtering___________####

vcf = readVcf(vcf_fname, genome = genome)

#For unique mutations
vcf_uniq = find_tree_sites(vcf, sample_name, bulk, min_absent = 1, min_present = 1, trisomy = "FALSE", gender, mut_type = "snv")

#Check if mutation is present in bulk
vaf = get_vaf(vcf_uniq, bulk)
vaf_f = vaf > 0

#Check if mutation is unique
gt = geno(vcf_uniq)$GT
samples_bulk = samples(header(vcf_uniq))
cols = seq(1, length(samples_bulk))
bulk_col = grep(bulk, samples_bulk)
samples_cols = cols[!cols %in% bulk_col]
clones_called = gt[ ,samples_cols] == "0/1" | gt[ ,samples_cols] == "1/1"
uniq_clone_f = rowSums(clones_called) == 1

vcf_uniq = vcf_uniq[vaf_f & uniq_clone_f,]

nr_shared = nrow(vcf_uniq)
print(paste0(nr_shared, " unique mutations were left after filtering"))
writeVcf(vcf_uniq, paste0(out_dir, "/", sample_name, "_unique_inbulk.vcf"))

#For shared mutations
vcf_shared = find_tree_sites(vcf, sample_name, bulk, min_absent = 1, min_present = 2, trisomy = "FALSE", gender, mut_type = "snv")
vaf = get_vaf(vcf_shared, bulk)
vaf_f = vaf > 0
vcf_shared = vcf_shared[vaf_f,]

nr_shared = nrow(vcf_shared)
print(paste0(nr_shared, " shared mutations were left after filtering"))
writeVcf(vcf_shared, paste0(out_dir, "/", sample_name, "_shared_inbulk.vcf"))









