"Usage:
    orthologous_soy_med.R (--out1 <O1>) (--out2 <O2>) (--out3 <O3>) <input1> <input2> <input3> <input4> <input5> <input6>
    orthologous_soy_med.R -h | --help

Options:
    -h --help   show this screen.
    <input1>    CNS classification of Medicago
    <input2>    CNS classification of Soybean
    <input3>    Soybean gene names
    <input4>    Medicago map of the names in the v5 and the names in the v4
    <input5>    List of othologous of Medicago and Soybean as listed by Plaza
    <input6>    Set of CNS
    --out1      File.
    --out2      File.
    --out3      File.

" -> doc

suppressMessages(require(vroom))
suppressMessages(require(tidyverse))
suppressMessages(require(docopt))
suppressMessages(require(splitstackshape) )

# retrieve the command-line arguments
opts <- docopt(doc)

medicago <- suppressMessages(vroom(opts$`<input1>`))
soybean <- suppressMessages(vroom(opts$`<input2>`))

CNS <- merge(medicago, soybean, by = "CNS_name")

CNS$same_class <- ifelse(CNS$CNS_classification.x == CNS$CNS_classification.y,
                         "Yes",
                         "No")

############################################################################
# Checking the orthologous only of the CNS with the correct classification #
############################################################################

CNS_correct <- filter(CNS, same_class == "Yes") %>%
    select(CNS_name, gene_name.x, gene_name.y, CNS_classification.y) %>%
    rename(gene_name_med = gene_name.x,
           ID = gene_name.y)

print("Where are the CNS that were classification in both genomes?")
print(table(CNS_correct$CNS_classification.y))

## Adds common name of soybean genes
names_soy <- suppressMessages(vroom(opts$`<input3>`))
names_soy <- cSplit(names_soy, 'OtherAliases', ',')

names_soy <- select(names_soy,
                    "ID", "Name", "OtherAliases_1") %>%
    rename(OtherAliases = OtherAliases_1)

CNS_correct2 <- merge(CNS_correct, names_soy, by = "ID")

## adds v4 names for medicago.
med_v4 <- suppressMessages(vroom(opts$`<input4>`))[, c(1,7)] %>%
    rename(gene_name_med = "#objectA",
           names_v4 = "objectB")

CNS_correct3 <- merge(CNS_correct2, med_v4, by = "gene_name_med", all.x = T )

# Using the list of orthologous "Orthologous gene families"
plaza <- suppressMessages(vroom(opts$`<input5>`,
                                col_names = F,
                                delim = "\t",
                                comment = "#"))

plaza_med <- plaza[plaza$X2 == "mtr", ]
plaza_med_soy <- plaza_med[plaza_med$X4 == "gma", c(1,3)] %>%
    rename(names_v4 = X1,
           soybean = X3) %>%
    group_by(names_v4) %>%
    summarise(ortho_list=paste(soybean, collapse=","))

CNS_correct4 <- merge(CNS_correct3,
                      plaza_med_soy,
                      by = "names_v4", all.x = T)
INCOMPLET <- CNS_correct4[!complete.cases(CNS_correct4), ]

## Some of the CNS doesn't have the IDs for both genomes
CNS_correct4 <- CNS_correct4[complete.cases(CNS_correct4), ]

write.table(CNS_correct4,
            opts$O1,
            col.names = F,
            row.names = F,
            quote = F,
            sep = "\t")

ortho_test <- data.frame()
for( c in 1:nrow(CNS_correct4)) {

    soy <- as.character(CNS_correct4$OtherAliases[c])
    ortho <- unique( strsplit(CNS_correct4$ortho_list, split = ",")[[c]] )

    Result <- as.data.frame(soy %in% ortho)

    ortho_test <- rbind(ortho_test, Result)
}

CNS_correct4$True_ortho <- ortho_test[, 1]

## Get the names of the CNS with duplication and select in a new data.frame
dup_in_orthogroup <- CNS_correct4[duplicated(CNS_correct4$CNS_name), ]
dup_in_orthogroup <- unique(dup_in_orthogroup$CNS_name)

dupl_df <- CNS_correct4[CNS_correct4$CNS_name %in% dup_in_orthogroup, ]
filtered <- CNS_correct4[!CNS_correct4$CNS_name %in% dup_in_orthogroup, ]

filtered_dups <- data.frame()
for(i in dup_in_orthogroup) {

    dup_sub <- dupl_df[dupl_df$CNS_name == i, ]

    dup_sub2 <- dup_sub[dup_sub$True_ortho == "TRUE", ]

    if( nrow(dup_sub2) > 0) {

        if(nrow(filtered_dups) == 0) {

            filtered_dups <- dup_sub2[1, ]

        } else {

            filtered_dups <- rbind(filtered_dups, dup_sub2[1, ])

        }

    } else {

        if(nrow(filtered_dups) == 0) {

            filtered_dups <- dup_sub[1, ]

        } else {

            filtered_dups <- rbind(filtered_dups, dup_sub[1, ])

        }
    }

}

final_ortho <- rbind(filtered_dups, filtered)
length(unique(final_ortho$CNS_name))

table(final_ortho$True_ortho)

write.table(final_ortho,
            opts$O2,
            col.names = F,
            row.names = F,
            quote = F,
            sep = "\t")

final_ortho_positive <- final_ortho[final_ortho$True_ortho == TRUE, ]

all_selected_CNS <- vroom(opts$`<input6>`,
                          col_names = F) %>%
    filter(X4 %in% final_ortho_positive$CNS_name)

write.table(all_selected_CNS,
            opts$O3,
            col.names = F,
            row.names = F,
            quote = F,
            sep = "\t")
