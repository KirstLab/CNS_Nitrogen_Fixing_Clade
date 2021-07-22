suppressMessages(require(vroom))
suppressMessages(require(ComplexHeatmap))
suppressMessages(require(tidyverse))
suppressMessages(require(circlize))
suppressMessages(require(docopt))

CNS <- suppressMessages(vroom("class_organized_medicago.txt", delim = "\t"))

final_CNS <- suppressMessages(vroom("nfix_cns/nfix_all_Nfix.all_final.bed",
                                    delim = "\t",
                                    col_names = F))

CNS <-  CNS %>%
    filter(CNS_name %in% final_CNS$X4) %>%
    mutate(CNS_classification = replace(CNS_classification,
                                        CNS_classification == "lrang_downstream",
                                        "distal downstream" ))

nod <- suppressMessages(vroom("nodulation_genes_medicago.csv",
                              col_names = T,
                              delim = ",")) %>%
    rename(gene_name = "IMGAG5.0") %>%
    filter( !is.na(gene_name) ) %>%
    unique()

scaled_rzm <- suppressMessages(vroom("input/ATAC_signal/joined_uniform_globalMeanLogQNorm_RZM.txt",
                                     delim = "\t") %>%
                                   rename(CNS_name = Name))

CNS_ATAC <- merge(CNS, scaled_rzm, by = "CNS_name") %>%
    select("CNS_name",
           "Control",
           "15min",
           "30min",
           "1hr",
           "2hr",
           "4hr",
           "8hr",
           "24hr")

assoc_genes <- suppressMessages(vroom("Mtrun_classified_CNS.bed",
                                      col_names = T)) %>%
    select(CNS_name, gene_name) %>%
    mutate(gene_name = gsub("_", "", gene_name))

CNS_ATAC_expr <- merge(CNS_ATAC,
                       assoc_genes,
                       by = "CNS_name",
                       all.x = T)

rna_seq1 <- suppressMessages( 
    vroom("gene_expression_tpm_qnorm_zeromeaned.tst",
          col_names = T,
          delim = "\t")
)

CNS_ATAC_expr2 <- merge(CNS_ATAC_expr,
                        rna_seq1,
                        by = "gene_name",
                        all.x = T) %>%
    arrange(CNS_name) %>%
    column_to_rownames("CNS_name") %>%
    filter(gene_name %in% nod$gene_name) %>%
    arrange(gene_name)

atac_to_ht <- select(CNS_ATAC_expr2,
                     c("Control.x",
                       "15min.x",
                       "30min.x",
                       "1hr.x",
                       "2hr.x",
                       "4hr.x",
                       "8hr.x",
                       "24hr.x"))

colnames(atac_to_ht) <- c("Control",
                          "15min",
                          "30min",
                          "1hr",
                          "2hr",
                          "4hr",
                          "8hr",
                          "24hr")

rnaseq_to_ht <- select(CNS_ATAC_expr2,
                       c("gene_name",
                         "Control.y",
                         "15min.y",
                         "30min.y",
                         "1hr.y",
                         "2hr.y",
                         "4hr.y",
                         "8hr.y",
                         "24hr.y"))

rnaseq_to_ht <- as.matrix(rnaseq_to_ht[, -1])
rownames(rnaseq_to_ht) <- CNS_ATAC_expr2$gene_name

colnames(rnaseq_to_ht) <- c("Control",
                            "15min",
                            "30min",
                            "1hr",
                            "2hr",
                            "4hr",
                            "8hr",
                            "24hr")

CNS_ATAC_groups <- merge(CNS, scaled_rzm, by = "CNS_name") %>%
    select("CNS_name",
           "CNS_classification") %>%
    rename(Class = CNS_classification)

col_fun = colorRamp2(c(-2.0, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2.0),
                     c("#0800ff", "#0033ff", "#0059ff", "#0084ff", "#FFFFFF", "#ff7300", "#FF6A6A", "#ff3300", "#ff0000"))

col_rna = colorRamp2(c(-3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0),
                     c("#05016e","#0600b0", "#0800ff", "#0033ff", "#0059ff", "#0084ff", "#FFFFFF", "#ff7300", "#FF6A6A", "#ff3300", "#ff0000","#b30000", "#990000"))

atac_lgd = Legend(col_fun = col_fun,
                  title = "ATAC-seq \n signal",
                  legend_height = unit(7, "cm"),
                  at = c(-2.0,
                         -1.5,
                         -1.0,
                         -0.5,
                         0,
                         0.5,
                         1.0,
                         1.5,
                         2.0))

rna_lgd = Legend(col_fun = col_rna,
                 title = "Expression",
                 legend_height = unit(7, "cm"),
                 at = c(-3.0,
                        -2.5,
                        -2.0,
                        -1.5,
                        -1.0,
                        -0.5,
                        0,
                        0.5,
                        1.0,
                        1.5,
                        2.0,
                        2.5,
                        3.0))

##@@@@@@@@@@@@@@@@@@@@##
## Heatmap of all CNS ##
##@@@@@@@@@@@@@@@@@@@@##

atac_to_ht_bar_annot <- atac_to_ht %>%
    mutate(CNS_name = row.names(.)) %>%
    select(CNS_name)

atac_to_ht_bar_annot <- plyr::join(atac_to_ht_bar_annot,
                                   CNS_ATAC_groups, by = "CNS_name")

sig <- suppressMessages(
    vroom("all_CNS_correlation.tst",
          col_names = T,
          delim = "\t") %>%
        mutate(class = ifelse(pvalue < 0.05, 1, 0)) %>%
        select(CNS, class) %>%
        rename(CNS_name = CNS)
)

atac_to_ht_bar_annot <- plyr::join(atac_to_ht_bar_annot,
                                   sig, by = "CNS_name") %>%
    column_to_rownames("CNS_name")

atac_to_ht_bar_annot[is.na(atac_to_ht_bar_annot)] <- 2

CNS.class_col <-c("#F0E442", # intronic
                  "#CC79A7", # upstream
                  "#E69F00", # distal upstream
                  "#56B4E9", # downstream
                  "#b56400", # distal downstream
                  "#009E73") # intergenic

names(CNS.class_col) <- c("intronic",
                          "upstream",
                          "distal upstream",
                          "downstream",
                          "distal downstream",
                          "intergenic")

class_colors <- list(Corr.sig = c("0"="grey60", "1"="black", "2"="white"),
                     CNS.class = CNS.class_col)
row_annot <- rowAnnotation(CNS.class = atac_to_ht_bar_annot$Class,
                           Corr.sig = as.factor(atac_to_ht_bar_annot$class),
                           col = class_colors)

ht <- Heatmap(atac_to_ht,
              column_title = "Time after LCO treatment",
              column_title_side = "bottom",
              cluster_columns = F,
              cluster_rows = F,
              col = col_fun,
              show_row_dend = T,
              show_heatmap_legend = F,
              show_row_names = T,
              row_names_side = "left",
              left_annotation = row_annot)

ht2 <- Heatmap(rnaseq_to_ht,
               column_title = "Time after LCO treatment",
               column_title_side = "bottom",
               cluster_columns = F,
               cluster_rows = F,
               col = col_fun,
               show_row_dend = T,
               show_heatmap_legend = F,
               show_row_names = T)

svg("nod_genes_ATCseq_RNAseq.svg", width = 17, height = 15, pointsize = 12)
draw(ht + ht2,
     annotation_legend_list = list(atac_lgd, rna_lgd))
dev.off()