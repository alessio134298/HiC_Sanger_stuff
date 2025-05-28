library(tidyverse)
library(argparser, quietly=TRUE)

parser <- arg_parser("Keep only Interchromosomal contacts")

parser <- add_argument(parser, "input", help="input")
parser <- add_argument(parser, "output", help="output")

chr_Graft <- c(paste0("hg38_", 1:22), "hg38_X", "hg38_Y", "gg6_1")
chr_Chicken <- c(paste0("gg6_", 1:28, 30:33), "gg6_W", "gg6_Z")
chr_Human <- c(paste0("hg38_", 1:22), "hg38_X", "hg38_Y")

args <- parse_args(parser)

Input <- args$input
output <- args$output

chr_Graft <- c(paste0("hg38_", 1:22), "hg38_X", "hg38_Y", "gg6_1")
chr_Chicken <- c(paste0("gg6_", c(1:28, 30:33)), "gg6_W", "gg6_Z")
chr_Human <- c(paste0("hg38_", 1:22), "hg38_X", "hg38_Y")

if (str_detect(Input, "DT40|SL29")) {
  chr = chr_Chicken
} else if (str_detect(Input, "HEK293T")) {
  chr = chr_Human
} else {chr = chr_Graft}

File <- read_delim(Input, delim = "\t") %>%
    dplyr::filter(chrom1 %in% chr & chrom2 %in% chr)

File_only_interchromosomal <- File %>% 
  dplyr::filter(chrom1 != chrom2)

