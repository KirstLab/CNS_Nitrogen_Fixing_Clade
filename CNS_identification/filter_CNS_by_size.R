"Usage:
    filter_CNS_by_size.R (--size <sz>) (--out1 <O1>) <input1>
    filter_CNS_by_size.R -h | --help

Options:
    -h --help   show this screen.
    <input1>    INPUT1  file containing the genomes information on NCBI
    --out1  List of all genomes and ftp links to download
    --size  Size in nucleotides
" -> doc

suppressMessages(require(docopt))
suppressMessages(require(vroom))

opts <- docopt(doc)

bed_file <- suppressMessages(vroom(opts$`<input1>`, 
                                   col_names = F))

size <- opts$sz
bed_file_filt <- bed_file[bed_file$X3 - bed_file$X2 >= size, ]

write.table(bed_file_filt,
            opts$O1,
            col.names = F,
            row.names = F,
            quote = F,
            sep = "\t")