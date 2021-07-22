"Usage:
    convert_cns_coords_to_bed.R (--out1 <O1>) <input1>
    convert_cns_coords_to_bed.R -h | --help

Options:
    -h --help   show this screen.
    <input1>    INPUT1
    --out1  OUTPUT1
" -> doc

suppressMessages(require(vroom))
suppressMessages(require(tidyverse))
suppressMessages(require(docopt))

opts <- docopt(doc)

coord <- vroom(opts$`<input1>`,
               col_names = F,
               delim = " ") %>%
    separate(X1, into = c("Chr", "b"), sep = ":") %>%
    separate(b, into = c("start", "end"), sep = "-") %>%
    mutate(start = as.numeric(start) - 1) %>%
    mutate(score = 0)
    
coord$strand <- "+"

coord_minus <- coord
coord_minus$strand <- "-"

coord <- rbind(coord, coord_minus)

write.table(coord, 
            paste(opts$O1),
            col.names = F,
            row.names = F,
            sep = "\t",
            quote = F)