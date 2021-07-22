"Usage:
    get_translation_start_and_end_point.R (--out1 <O1>) <input1>
    get_translation_start_and_end_point.R -h | --help

Options:
    -h --help   show this screen.
    <input1>    INPUT1
    --out1  OUTPUT1
" -> doc

suppressMessages(require(vroom))
suppressMessages(require(tidyverse))
suppressMessages(require(docopt))

opts <- docopt(doc)

gff_file <- vroom(opts$`<input1>`,
                  col_names = F) %>%
    filter(X3 == "CDS") %>%
    separate(X9, into = c("a", "b"), sep = "locus_tag \"") %>%
    separate(b, into = c("name", "remove"), sep = "\";") %>%
    select(-a, -remove, -X2, -X3, -X6, -X8 )

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
    
    plus_sub_df <- data.frame(chr, TSS, TSE, name, 0, strand)
    
    plus_sub_final <- rbind(plus_sub_final, plus_sub_df)
    
}

gff_minus <- gff_file[gff_file$X7 == "-", ]
minus_ltags <- unique(gff_minus$name)

minus_sub_final <- data.frame()
for (i in minus_ltags) {
    
    minus_sub <- gff_minus[gff_minus$name == i, ]
    
    chr <- unique(minus_sub[, 1])
    TSS <- max(minus_sub[, 3]) - 1
    TSE <- min(minus_sub[, 2])
    name <- paste0(unique(minus_sub[, 5]))
    strand <- unique(minus_sub[, 4])
    
    minus_sub_df <- data.frame(chr, TSE, TSS, name, 0, strand)

    minus_sub_final <- rbind(minus_sub_final, minus_sub_df)
    
}

names(plus_sub_final)[2:3] <- c("start", "end")
names(minus_sub_final)[2:3] <- c("start", "end")

combined <- rbind(plus_sub_final, minus_sub_final)

write.table(combined,
            paste(opts$O1),
            col.names = F,
            row.names = F,
            quote = F,
            sep = "\t")