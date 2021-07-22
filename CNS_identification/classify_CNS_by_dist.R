"Usage:
    classify_CNS_by_dist.R (--out1 <O1>) (--out2 <O2>) <input1>
    classify_CNS_by_dist.R -h | --help

Options:
    -h --help   show this screen.
    <input1>    INPUT1
    --out1  OUTPUT1
    --out2  OUTPUT2
" -> doc

suppressMessages(require(vroom))
suppressMessages(require(tidyverse))
suppressMessages(require(docopt))

opts <- docopt(doc)

dist_to_CDS <- suppressMessages(vroom(opts$`<input1>`,
                                      delim = "\t",
                                      col_names = F))

dist_to_CDS$class <- ifelse(dist_to_CDS$X13 >= 2000 & dist_to_CDS$X13 < 10000,
                            "lrang_downstream",
                            ifelse(dist_to_CDS$X13 > 0 & dist_to_CDS$X13 < 2000,
                                   "downstream",
                                   ifelse(dist_to_CDS$X13 == 0, "intronic",
                                          ifelse(dist_to_CDS$X13 < 0 & dist_to_CDS$X13 > -2000,
                                                 "upstream",
                                                 ifelse(dist_to_CDS$X13 <= -2000 & dist_to_CDS$X13 > -10000,
                                                        "lrang_upstream",
                                                        "intergenic")))))

filtered <- dist_to_CDS[!dist_to_CDS$class == "intergenic", ]

unclassified <- dist_to_CDS[dist_to_CDS$class == "intergenic", ]
write.table(unclassified,
            paste(opts$O1),
            col.names = F,
            row.names = F,
            sep = "\t",
            quote = F)

CNS_in_both_strands <- filtered[duplicated(filtered$X4) == TRUE, ] %>%
    select(X4)

CNS_in_both_strands <- unique(CNS_in_both_strands$X4)

dupl_df <- filtered[filtered$X4 %in% CNS_in_both_strands, ]
filtered <- filtered[!filtered$X4 %in% CNS_in_both_strands, ]

selected_from_dup <- data.frame()
for(i in CNS_in_both_strands){

    dupl_df_sub <- dupl_df[dupl_df$X4 == i, ]

    dupl_df_sub <- dupl_df_sub[abs(dupl_df_sub$X13) == min(abs(dupl_df_sub$X13)), ]

    selected_from_dup <- rbind(selected_from_dup, dupl_df_sub)

}

final_CNS_classified <- rbind(selected_from_dup, filtered) %>%
    select(-X5, -X7, -X11)

colnames(final_CNS_classified) <- c("Chr",
                                    "CNS_start",
                                    "CNS_end",
                                    "CNS_name",
                                    "CNS_strand",
                                    "gene_start",
                                    "gene_end",
                                    "gene_name",
                                    "gene_strand",
                                    "distance",
                                    "CNS_classification")

## Check if is there any CNS that have duplicated entrance (i.e. same distance of features in both strands)
CNS_dup <- final_CNS_classified[duplicated(final_CNS_classified$CNS_name), ]
final_CNS_classified_dup <- final_CNS_classified[final_CNS_classified$CNS_name %in% CNS_dup$CNS_name, ]
final_CNS_classified <- final_CNS_classified[!final_CNS_classified$CNS_name %in% CNS_dup$CNS_name, ]

# After investigating the gene models that overlaps, excludes the ones that seems to be in the incorrectly annotated gene
to_keep <- c("MtrunA17_Chr2g0315481", "MtrunA17_Chr2g0298731")

final_CNS_classified_dup <- filter(final_CNS_classified_dup,
                                   gene_name %in% to_keep)

final_CNS_classified <- rbind(final_CNS_classified, final_CNS_classified_dup)

write.table(final_CNS_classified,
            paste(opts$O2),
            col.names = T,
            row.names = F,
            quote = F,
            sep = "\t")

table(final_CNS_classified$CNS_classification)
length(unique(final_CNS_classified$CNS_name))
