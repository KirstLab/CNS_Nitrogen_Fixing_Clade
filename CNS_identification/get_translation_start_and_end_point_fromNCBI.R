require(vroom)
require(tidyverse)
require(reutils)

GFF <- "./GCF_000004515.5_Glycine_max_v2.1_genomic.gff"

gff_file <- vroom(GFF,
                  col_names = F,
                  delim = "\t",
                  comment = "#") %>%
    filter(X3 == "CDS") %>%
    separate(X9, into = c("a", "b"), sep = "GeneID:") %>%
    separate(b, into = c("name", "remove"), sep = ",") %>%
    select(-a, -remove) %>%
    separate(name, into = c("name", "remove"), sep = ";") %>%
    select(-remove) %>%
    select(-X2, -X3, -X6, -X8 )

## plus strand
gff_plus <- gff_file[gff_file$X7 == "+", ]
plus_ltags <- unique(gff_plus$name)

plus_sub_final <- data.frame()
for (i in plus_ltags) {

    plus_sub <- gff_plus[gff_plus$name == i, ]

    chr <- unique(plus_sub[, 1])
    TSS <- min(plus_sub[, 2]) - 1
    TSE <- max(plus_sub[, 3])
    name <- paste0(unique(plus_sub[, 5]))
    strand <- unique(plus_sub[, 4])

    plus_sub_df <- data.frame(chr, TSS, TSE, name, strand)

    plus_sub_final <- rbind(plus_sub_final, plus_sub_df)

}

## minus strand
gff_minus <- gff_file[gff_file$X7 == "-", ]
minus_ltags <- unique(gff_minus$name)

minus_sub_final <- data.frame()
for (i in minus_ltags) {

    minus_sub <- gff_minus[gff_minus$name == i, ]

    chr <- unique(minus_sub[, 1])
    TSS <- min(minus_sub[, 2])
    TSE <- max(minus_sub[, 3])
    name <- paste0(unique(minus_sub[, 5]))
    strand <- unique(minus_sub[, 4])

    minus_sub_df <- data.frame(chr, TSS, TSE, name, strand)

    minus_sub_final <- rbind(minus_sub_final, minus_sub_df)

}

combined <- rbind(plus_sub_final, minus_sub_final)
combined$score <- 0

combined <- combined[, c(1,2,3,4,6,5)]

write.table(combined,
            "Glycine_max_CDS_TSS_to_TSE.bed",
            col.names = F,
            row.names = F,
            quote = F,
            sep = "\t")

## downloading the ids from ncbi to convert to the IDs acceped by plaza.

# Entrez only accepts 500 IDs per query.
gene_ID <- combined$name
max <- 500
x <- seq_along(gene_ID)

gene_ID2 <- split(gene_ID, ceiling(x/max))

gene_ID2_df <- data.frame()
for (i in 1:length(gene_ID2)) {

    # query NCBI and save the data in a list (since each query can return different number of fields)
    sub_set <- esummary(as.character(gene_ID2[[i]]),
                        db = "gene",
                        querykey = c("Name"))

    sub_set_df <- content(sub_set, as = "parsed")

    ids_DF <- data.frame()
    for( c in 1:length(sub_set_df)){

        ids <- data.frame(ID = names(sub_set_df[c]),
                          Name = sub_set_df[[c]]["Name"],
                          Name = sub_set_df[[c]]["OtherAliases"])

        ids_DF <- rbind(ids_DF, ids)

    }

    gene_ID2_df <- rbind(gene_ID2_df, ids_DF)

}

write.table(gene_ID2_df,
            "soybean_genes_names.tsv",
            col.names = T,
            row.names = F,
            quote = F,
            sep = "\t")
