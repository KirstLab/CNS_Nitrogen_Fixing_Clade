"Usage:
    nfix_complex_heat_map_CNS.R (--out1 <O1>) (--out2 <O2>) (--out3 <O3>) (--out4 <O4>) (--out5 <O5>) (--out6 <O6>) (--out7 <O7>) (--out8 <O8>) (--out9 <O9>)  (--out10 <O10>) (--out11 <O11>) (--out12 <O12>) (--out13 <O13>) (--out14 <O14>) (--out15 <O15>) <input1> <input2> <input3> <input4> <input5>
    nfix_complex_heat_map_CNS.R -h | --help

Options:
    -h --help   show this screen.
    <input1>    INPUT1
    <input2>    INPUT2
    <input3>    INPUT3
    <input4>    INPUT4
    <input5>    INPUT5
    --out1  OUTPUT1
    --out2  OUTPUT2
    --out3  OUTPUT3
    --out4  OUTPUT4
    --out5  OUTPUT5
    --out6  OUTPUT6
    --out7  OUTPUT7
    --out8  OUTPUT8
    --out9  OUTPUT9
    --out10  OUTPUT10
    --out11  OUTPUT11
    --out12  OUTPUT12
    --out13  OUTPUT13
    --out14  OUTPUT14
    --out15  OUTPUT15
" -> doc

suppressMessages(require(vroom))
suppressMessages(require(ComplexHeatmap))
suppressMessages(require(tidyverse))
suppressMessages(require(circlize))
suppressMessages(require(docopt))

opts <- docopt(doc)

filtered_CNS <- suppressMessages(vroom("selected_CNS_same_orthogroup.bed",
                                       delim = "\t",
                                       col_names = F))

CNS <- suppressMessages(vroom(opts$input1, delim = "\t")) %>%
    filter(CNS_name %in% filtered_CNS$X4) %>%
    mutate(CNS_classification = replace(CNS_classification,
                                        CNS_classification == "lrang_downstream",
                                        "distal downstream" ))

nfix_CNS <- suppressMessages(vroom(opts$input5,
                                   col_names = F,
                                   delim = "\t"))

CNS <- CNS %>%
    filter(CNS_name %in% nfix_CNS$X4)

scaled_rzm <- suppressMessages(vroom(opts$input2,
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

assoc_genes <- suppressMessages(vroom(opts$input3,
                                      col_names = T)) %>%
    select(CNS_name, gene_name) %>%
    mutate(gene_name = gsub("_", "", gene_name))

CNS_ATAC_expr <- merge(CNS_ATAC,
                       assoc_genes,
                       by = "CNS_name",
                       all.x = T)

rna_seq1 <- suppressMessages( 
    vroom(opts$input4,
          col_names = T,
          delim = "\t")
)

CNS_ATAC_expr2 <- merge(CNS_ATAC_expr,
                        rna_seq1,
                        by = "gene_name",
                        all.x = T) %>%
    arrange(CNS_name) %>%
    column_to_rownames("CNS_name")


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
                       c("Control.y",
                         "15min.y",
                         "30min.y",
                         "1hr.y",
                         "2hr.y",
                         "4hr.y",
                         "8hr.y",
                         "24hr.y"))

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
    rename(Class = CNS_classification) %>%
    column_to_rownames("CNS_name")

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

ht <- Heatmap(atac_to_ht,
              column_title = "Time after LCO treatment",
              column_title_side = "bottom",
              cluster_columns = F,
              col = col_fun,
              show_row_dend = T,
              show_heatmap_legend = F,
              show_row_names = F)

colors_ht2 = structure(1:6, names = c("intronic",
                                      "upstream",
                                      "distal upstream",
                                      "downstream",
                                      "distal downstream",
                                      "intergenic"))

colors_ht2 <-c("#F0E442", # intronic
               "#CC79A7", # upstream
               "#E69F00", # distal upstream
               "#56B4E9", # downstream
               "#b56400", # distal downstream
               "#009E73") # intergenic

names(colors_ht2) <- c("intronic",
                       "upstream",
                       "distal upstream",
                       "downstream",
                       "distal downstream",
                       "intergenic")

order_ht <- row_order(ht)

CNS_ATAC_groups_ht <- rnaseq_to_ht[order_ht,]
CNS_ATAC_groups_ht <- CNS_ATAC_groups_ht %>%
    mutate(CNS = rownames(.)) %>%
    select(CNS)

CNS_ATAC_groups2 <- CNS_ATAC_groups %>%
    mutate(CNS = rownames(.))

CNS_ATAC_groups_ht <- plyr::join(CNS_ATAC_groups_ht,
                                 CNS_ATAC_groups2,
                                 by = "CNS") %>%
    select("Class")

ht2 <- Heatmap(CNS_ATAC_groups_ht,
               name = "class",
               col = colors_ht2,
               width = 0.4,
               show_heatmap_legend = F,
               column_title = "",
               show_row_names = F)

lgd2 = Legend(labels = c("intronic",
                         "upstream",
                         "distal upstream",
                         "downstream",
                         "distal downstream",
                         "intergenic"),
              legend_gp = gpar(fill = colors_ht2),
              title = "Class",
              legend_height = unit(3, "cm"))

svg(opts$O1, width = 15, height = 30, pointsize = 12)
draw(ht + ht2,
     annotation_legend_list = list(atac_lgd, lgd2))
dev.off()

####@@@@@@@@@@@@@@@@@@@@@@####
### Heatmaps of each class ###
####@@@@@@@@@@@@@@@@@@@@@@####

cns_names <- filter(CNS_ATAC_groups,
                    Class != "intronic") %>%
    rownames()

cns_df <- filter(atac_to_ht,
                 rownames(atac_to_ht) %in% cns_names)

TEST <- CNS_ATAC_expr2[rownames(CNS_ATAC_expr2) %in% cns_names, ]
GENES <- CNS_ATAC_expr[CNS_ATAC_expr$CNS_name %in% cns_names, ]
length(unique(GENES$gene_name))
length(unique(GENES$CNS_name))

## Help function
ATAC_ht <- function(CNS_class = CNS_class) {
    
    cns_names <- filter(CNS_ATAC_groups,
                        Class == CNS_class) %>%
        rownames()
    
    cns_df <- filter(atac_to_ht,
                     rownames(atac_to_ht) %in% cns_names)
    
    Heatmap(cns_df,
            column_title = "Time after LCO treatment",
            column_title_side = "bottom",
            cluster_columns = F,
            col = col_fun,
            show_row_dend = T,
            show_heatmap_legend = F,
            show_row_names = F)
    
}

RNAseq_ht <- function(CNS_class = CNS_class,
                      order_ht = order_ht){
    
    cns_names <- filter(CNS_ATAC_groups,
                        Class == CNS_class) %>%
        rownames()
    
    cns_df <- filter(rnaseq_to_ht,
                     rownames(rnaseq_to_ht) %in% cns_names)
    
    cns_df2 <- cns_df[order_ht, ]
    
    Heatmap(cns_df,
            column_title = "Time after LCO treatment",
            column_title_side = "bottom",
            cluster_columns = F,
            cluster_rows = F,
            col = col_rna,
            show_row_dend = F,
            show_heatmap_legend = F,
            show_row_names = F)
    
}

upstream_cns <- ATAC_ht(CNS_class = "upstream")

upstream_cns_rna <- RNAseq_ht(CNS_class = "upstream",
                              order_ht = row_order(upstream_cns))

svg(opts$O2, width = 15, height = 30, pointsize = 12)
draw(upstream_cns + upstream_cns_rna,
     annotation_legend_list = list(atac_lgd, rna_lgd))
dev.off()

intronic_cns <- ATAC_ht(CNS_class = "intronic")

intronic_cns_rna <- RNAseq_ht(CNS_class = "intronic",
                              order_ht = row_order(intronic_cns))

svg(opts$O3, width = 15, height = 30, pointsize = 12)
draw(intronic_cns + intronic_cns_rna,
     annotation_legend_list = list(atac_lgd, rna_lgd))
dev.off()

downstream_cns <- ATAC_ht(CNS_class = "downstream")

downstream_cns_rna <- RNAseq_ht(CNS_class = "downstream",
                                order_ht = row_order(downstream_cns))

svg(opts$O4, width = 15, height = 30, pointsize = 12)
draw(downstream_cns + downstream_cns_rna,
     annotation_legend_list = list(atac_lgd, rna_lgd))
dev.off()

distal_upstream_cns <- ATAC_ht(CNS_class = "distal upstream")

distal_upstream_cns_rna <- RNAseq_ht(CNS_class = "distal upstream",
                                     order_ht = row_order(distal_upstream_cns))

svg(opts$O5, width = 15, height = 30, pointsize = 12)
draw(distal_upstream_cns + distal_upstream_cns_rna,
     annotation_legend_list = list(atac_lgd, rna_lgd))
dev.off()

distal_downstream_cns <- ATAC_ht(CNS_class = "distal downstream")

distal_downstream_cns_rna <- RNAseq_ht(CNS_class = "distal downstream",
                                       order_ht = row_order(distal_downstream_cns))

svg("nfix_all_CNS_heatmap_distal_downstream.svg", width = 15, height = 30, pointsize = 12)
draw(distal_downstream_cns + distal_downstream_cns_rna,
     annotation_legend_list = list(atac_lgd, rna_lgd))
dev.off()

####@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@####
### Heatmaps of each class after correlation filter ###
####@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@####

upst_or_down_CNS <- CNS_ATAC_groups %>%
    filter(Class %in% c("upstream", "downstream", "distal upstream", "distal downstream" )) %>%
    rownames()

orig_corr2 <- vroom("all_CNS_correlation.tst",
                    col_names = T,
                    delim = "\t") %>%
    filter(CNS %in% upst_or_down_CNS) %>%
    column_to_rownames("CNS")

orig_corr2 <- orig_corr2[!is.na(orig_corr2$pvalue), ]

orig_corr2$color_code <- ifelse(orig_corr2$pvalue < 0.05 & abs(orig_corr2$orig_corr) >= 0.5,
                                "sig",
                                "not_sig")

orig_corr2$direction <- ifelse(orig_corr2$color_code == "sig", 
                               ifelse(orig_corr2$orig_corr < 0, "Negative", "Positive"), "NA")

print("The set of all CNS (up or down a gene) is classified as:")
table(orig_corr2$color_code)
table(orig_corr2$direction)

corr_plot_sig <- function (dataset = dataset, CNS_class = "all") {
    
    if (CNS_class == "all") {
        
        print(table(dataset$color_code))
        print(table(dataset$direction))
        
        ggplot(dataset, aes(x=orig_corr, fill= color_code)) +
            geom_histogram(binwidth = 0.025, color = "black", alpha = 0.7) +
            ylab("") +
            xlab("") +
            xlim(-1, 1) +
            theme_bw() +
            theme(axis.text = element_text(size = 22),
                  axis.title = element_text(size = 22),
                  legend.position = "") +
            scale_fill_manual( values = c("#0021A5", "#FA4616"))
        
    } else {
        
        cns_names <- filter(CNS_ATAC_groups,
                            Class == CNS_class) %>%
            rownames()
        
        cns_df <- filter(dataset,
                         rownames(dataset) %in% cns_names)
        
        print(table(cns_df$color_code))
        print(table(cns_df$direction))
        
        ggplot(cns_df, aes(x=orig_corr, fill= color_code)) +
            geom_histogram(binwidth = 0.025, color = "black", alpha = 0.7) +
            ylab("") +
            xlab("") +
            xlim(-1, 1) +
            theme_bw() +
            theme(axis.text = element_text(size = 22),
                  axis.title = element_text(size = 22),
                  legend.position = "") +
            scale_fill_manual( values = c("#0021A5", "#FA4616"))
        
    }
    
}


svg(opts$O6, width = 5, height = 6, pointsize = 12)
corr_plot_sig(orig_corr2)
dev.off()

tokeep_CNS_pos <- orig_corr2 %>%
    filter(color_code == "sig") %>%
    filter(orig_corr >= 0.5) %>%
    rownames()

tokeep_CNS_neg <- orig_corr2 %>%
    filter(color_code == "sig") %>%
    filter(orig_corr <= 0.5) %>%
    rownames()

rnaseq_to_ht_cor_pos <- rnaseq_to_ht[rownames(rnaseq_to_ht) %in% tokeep_CNS_pos, ]
rnaseq_to_ht_cor_neg <- rnaseq_to_ht[rownames(rnaseq_to_ht) %in% tokeep_CNS_neg, ]

atac_to_ht_cor_pos <- atac_to_ht[rownames(atac_to_ht) %in% tokeep_CNS_pos, ]
atac_to_ht_cor_neg <- atac_to_ht[rownames(atac_to_ht) %in% tokeep_CNS_neg, ]

ATAC_ht_corr <- function(CNS_class = CNS_class, dataset = dataset) {
    
    cns_names <- filter(CNS_ATAC_groups,
                        Class == CNS_class) %>%
        rownames()
    
    cns_df <- filter(dataset,
                     rownames(dataset) %in% cns_names)
    
    Heatmap(cns_df,
            column_title = "Time after LCO treatment",
            column_title_side = "bottom",
            cluster_columns = F,
            col = col_fun,
            show_row_dend = T,
            show_heatmap_legend = F,
            show_row_names = F)
    
}

RNAseq_ht_corr <- function(CNS_class = CNS_class,
                           order_ht = order_ht,
                           dataset = dataset){
    
    cns_names <- filter(CNS_ATAC_groups,
                        Class == CNS_class) %>%
        rownames()
    
    cns_df <- filter(dataset,
                     rownames(dataset) %in% cns_names)
    
    cns_df2 <- cns_df[order_ht, ]
    
    Heatmap(cns_df,
            column_title = "Time after LCO treatment",
            column_title_side = "bottom",
            cluster_columns = F,
            cluster_rows = F,
            col = col_rna,
            show_row_dend = F,
            show_heatmap_legend = F,
            show_row_names = F)
    
}

## Ploting the significance of the correlations in each subset ##
svg(opts$O7, width = 6, height = 5, pointsize = 12)
corr_plot_sig(orig_corr2, CNS_class = "upstream")
dev.off()

svg(opts$O8, width = 6, height = 5, pointsize = 12)
corr_plot_sig(orig_corr2, CNS_class = "downstream")
dev.off()

svg(opts$O9, width = 6, height = 5, pointsize = 12)
corr_plot_sig(orig_corr2, CNS_class = "distal upstream")
dev.off()

svg("nfix_correlation_significance_distal_downstream_all.svg", width = 6, height = 5, pointsize = 12)
corr_plot_sig(orig_corr2, CNS_class = "distal downstream")
dev.off()

### Ploting the heatmaps of the positive correlation

upstream_cns_pos <- ATAC_ht_corr(CNS_class = "upstream",
                                 dataset = atac_to_ht_cor_pos)

upstream_cns_rna_pos <- RNAseq_ht_corr(CNS_class = "upstream",
                                       order_ht = row_order(upstream_cns_pos),
                                       dataset = rnaseq_to_ht_cor_pos)

svg(opts$O10, width = 6, height = 12, pointsize = 12)
draw(upstream_cns_pos + upstream_cns_rna_pos,
     annotation_legend_list = list(atac_lgd, rna_lgd))
dev.off()

downstream_cns_pos <- ATAC_ht_corr(CNS_class = "downstream",
                                   dataset = atac_to_ht_cor_pos)

downstream_cns_rna_pos <- RNAseq_ht_corr(CNS_class = "downstream",
                                         order_ht = row_order(downstream_cns_pos),
                                         dataset = rnaseq_to_ht_cor_pos)

svg(opts$O11, width = 6, height = 12, pointsize = 12)
draw(downstream_cns_pos + downstream_cns_rna_pos,
     annotation_legend_list = list(atac_lgd, rna_lgd))
dev.off()

distal_upstream_cns_pos <- ATAC_ht_corr(CNS_class = "distal upstream",
                                        dataset = atac_to_ht_cor_pos)

distal_upstream_cns_rna_pos <- RNAseq_ht_corr(CNS_class = "distal upstream",
                                              order_ht = row_order(distal_upstream_cns_pos),
                                              dataset = rnaseq_to_ht_cor_pos)

svg(opts$O12, width = 6, height = 12, pointsize = 12)
draw(distal_upstream_cns_pos + distal_upstream_cns_rna_pos,
     annotation_legend_list = list(atac_lgd, rna_lgd))
dev.off()

### Ploting the heatmaps of the neagtive correlation

upstream_cns_neg <- ATAC_ht_corr(CNS_class = "upstream",
                                 dataset = atac_to_ht_cor_neg)

upstream_cns_rna_neg <- RNAseq_ht_corr(CNS_class = "upstream",
                                       order_ht = row_order(upstream_cns_neg),
                                       dataset = rnaseq_to_ht_cor_neg)

svg(opts$O13, width = 6, height = 12, pointsize = 12)
draw(upstream_cns_neg + upstream_cns_rna_neg,
     annotation_legend_list = list(atac_lgd, rna_lgd))
dev.off()

downstream_cns_neg <- ATAC_ht_corr(CNS_class = "downstream",
                                   dataset = atac_to_ht_cor_neg)

downstream_cns_rna_neg <- RNAseq_ht_corr(CNS_class = "downstream",
                                         order_ht = row_order(downstream_cns_neg),
                                         dataset = rnaseq_to_ht_cor_neg)

svg(opts$O14, width = 6, height = 12, pointsize = 12)
draw(downstream_cns_neg + downstream_cns_rna_neg,
     annotation_legend_list = list(atac_lgd, rna_lgd))
dev.off()

distal_upstream_cns_neg <- ATAC_ht_corr(CNS_class = "distal upstream",
                                        dataset = atac_to_ht_cor_neg)

distal_upstream_cns_rna_neg <- RNAseq_ht_corr(CNS_class = "distal upstream",
                                              order_ht = row_order(distal_upstream_cns_neg),
                                              dataset = rnaseq_to_ht_cor_neg)

svg(opts$O15, width = 6, height = 12, pointsize = 12)
draw(distal_upstream_cns_neg + distal_upstream_cns_rna_neg,
     annotation_legend_list = list(atac_lgd, rna_lgd))
dev.off()