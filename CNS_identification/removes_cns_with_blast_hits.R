"Usage:
    removes_cns_with_blast_hits.R [--length <lgt>] [--evalue <ev>] (--out1 <O1>) <input1> <input2>
    removes_cns_with_blast_hits.R -h | --help

Options:
    -h --help   show this screen.
    <input1>    INPUT1  bed file containing the CNS information.
    <input2>    INPUT2  file containing the blastx hits.
    --out1  bed file containing the filtered CNS.
    --evalue=EVALUE threshold to consider a hit significant enought to exclude the CNS [default: 0.01]
    --length=LENGTH minimum size of the alignment to consider a hit [default: 0].
" -> doc

suppressMessages(require(docopt))
suppressMessages(require(vroom))
suppressMessages(require(tidyverse))

opts <- docopt(doc)
evalue <- as.numeric(opts$evalue)
alignment_length <- as.numeric(opts$length)

cns_file <- suppressMessages(vroom(opts$input1, col_names = F))
blast_out <- suppressMessages(vroom(opts$input2, col_names = F))

colnames(blast_out) <- c("query id",
                         "subject id",
                         "alignment length",
                         "query length",
                         "subject length",
                         "q. start",
                         "q. end",
                         "s. start",
                         "s. end",
                         "evalue")

blast_out_fit <- blast_out[blast_out$evalue <= evalue, ]
blast_out_fit <- blast_out_fit[blast_out_fit$`alignment length` >= alignment_length, ]

blast_out_fit <- separate(blast_out_fit, 
                          "query id",
                          c("cns_name", "cns_coord"),
                          sep = "::")

cns_file_final <- cns_file[!cns_file$X4 %in% blast_out_fit$cns_name, ]

write.table(cns_file_final,
            opts$O1,
            col.names = F,
            row.names = F,
            quote = F,
            sep = "\t")