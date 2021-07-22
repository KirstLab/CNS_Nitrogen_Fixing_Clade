suppressMessages(require(vroom))
suppressMessages(require(tidyverse))
suppressMessages(require(docopt))

medicago_class <- suppressMessages(vroom("Mtrun_classified_CNS.bed",
                                         col_names = T)) %>%
    filter(CNS_classification != "intergenic") %>%
    select(CNS_name, CNS_classification) %>%
    mutate_if(is.character,
              str_replace_all,
              pattern = "lrang_upstream",
              replacement = "distal upstream") %>%
    mutate_if(is.character,
              str_replace_all,
              pattern = "lrang_downstream",
              replacement = "distal downstream")

medicago_inter <- suppressMessages(vroom("Mtrun_unclassified_CNS.bed",
                                         col_names = F)) %>%
    filter(!X4 %in% medicago_class$CNS_name) %>%
    select(X4) %>%
    unique()

medicago_inter <- as.data.frame(medicago_inter) %>%
    mutate(CNS_classification = "intergenic") %>%
    rename(CNS_name = X4)

medicago <- rbind(medicago_class, medicago_inter)

Soybean_class <- suppressMessages(vroom("Soybean_classified_CNS.bed",
                                        col_names = T)) %>%
    filter(CNS_classification != "intergenic") %>%
    select(CNS_name, CNS_classification) %>%
    mutate_if(is.character,
              str_replace_all,
              pattern = "lrang_upstream",
              replacement = "distal upstream") %>%
    mutate_if(is.character,
              str_replace_all,
              pattern = "lrang_downstream",
              replacement = "distal downstream")

Soybean_inter <- suppressMessages(vroom("Soybean_unclassified_CNS.bed",
                                        col_names = F)) %>%
    filter(!X4 %in% Soybean_class$CNS_name) %>%
    select(X4) %>%
    unique()

Soybean_inter <- as.data.frame(Soybean_inter) %>%
    mutate(CNS_classification = "intergenic") %>%
    rename(CNS_name = X4)

Soybean <- rbind(Soybean_class, Soybean_inter)

write.table(medicago,
            "class_organized_medicago.txt",
            col.names = T,
            row.names = F,
            quote = F,
            sep = "\t")

write.table(Soybean,
            "class_organized_Soybean.txt",
            col.names = T,
            row.names = F,
            quote = F,
            sep = "\t")

comb_med_soy <- merge(medicago, Soybean, by = "CNS_name") %>%
    rename(Medicago = CNS_classification.x,
           Soybean = CNS_classification.y) %>%
    mutate(CNS_name = "CNS")

comb_med_soy$same_class <- ifelse(comb_med_soy$Medicago == comb_med_soy$Soybean,
                         "Yes",
                         "No")

write.table(comb_med_soy,
            "ALL_CNS_classification_MED_vs_SOY.tst",
            col.names = T,
            row.names = T,
            quote = F,
            sep = "\t")

### Sankey plot ###

require(flipPlots)
require(grid)
library(webshot)
webshot::install_phantomjs()

comb_med_soy <- comb_med_soy[,-4]

SankeyDiagram(comb_med_soy,
              link.color = "Source",
              weights = rep(3, nrow(comb_med_soy)),
              label.show.varname = F,
              font.size = 24,
              label.show.percentages = T,
              variables.share.values = T,
              label.show.counts = F,
              sinks.right =  F,
              node.position.automatic = F)

#webshot("./sankey.html" , "output.pdf", delay = 0.2)
