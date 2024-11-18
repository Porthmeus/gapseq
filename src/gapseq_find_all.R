# Porthmeus
# 25.10.24

# this script takes a output from a diamond blast and reformats it into the expected form of the pathway and reactions tables of gapseq

# libraries
library(data.table)
library(digest)

# arguments and variables
args <- commandArgs(trailingOnly = TRUE)
diamond_out <- args[1]
taxonomy <- args[2]
taxRange <- args[3]
seqSrc <- args[4]
srcDir <- args[5]
pwyKey <- "Pathways|Enzyme-Test|seed|kegg" # corresponds to -p all


dmnd_colnames <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

#### debug ####
diamond_out <- "./ecoli.fna_prodigal.fasta_blastresult.tsv"
taxonomy <- "Bacteria"
taxRange <- "all"
srcDir <-"./src"
seqSrc <- 2



# code
# load pathway databases
pwys <- c("meta","kegg","seed","custom")
pwfls <- file.path(srcDir, "..","dat", paste0(pwys,"_pwy.tbl"))
names(pwfls)<- pwys
pathways <- lapply(names(pwfls),function(x) cbind(source = x, fread(pwfls[x])))
pathways[["fill"]] <- TRUE
pathways <- do.call(rbind, pathways)

# remove duplicated pathways
dup_ids <- pathways[duplicated(id),id]
pathways <-pathways[!(id %in% dup_ids) | source == "custom",]

# helper function to split correctly
splitBy <- function(x, split){
    # similar to strsplit but applicable only to a single string. Will add a "" in the end of the output vector if the string ends with a a split character
    x_split <- gsub(paste0("(",split,")"), "\\1__",x)
    x_split <- strsplit(x_split, split = split)[[1]]
    x_split <- gsub("__","",x_split)
    return(x_split)
}

# add empty string for reaEc when there is a reacId but no ec number
pathways[reaEc == "" & reaId != "",reaEc := "__"]

# get a table with the pathways split by individual reactions
reaEc <- sapply(pathways[,reaEc], function(x) splitBy(x, split = "[,;]"))
reaId <- sapply(pathways[,reaId], function(x) splitBy(x, split = ","))
reaName <- sapply(pathways[,reaName], function(x) splitBy(x, split = ";"))
lengthEc <- sapply(reaEc, length)
lengthId<- sapply(reaId, length)
lengthName<- sapply(reaName, length)

# check if all vectors are of the same length, if not report an error
if(any(lengthEc != lengthId|lengthEc != lengthName|lengthId != lengthName)){
    stop("Check database input as reactions, ec numbers and reaction names have different length")
}

# get the right length for the pathway id
pathwayId <-sapply(1:length(pathways[,id]), function(x) rep(pathways[x,id],lengthEc[x]))

rxns <- data.table(id = unlist(pathwayId),
                   reaEc = unlist(reaEc),
                   reaId = unlist(reaId),
                   reaName = unlist(reaName))

#sel <- 35760:35764
#rxns <- rxns[sel,]
#rxns[,.(reaEc = sapply(reaEc, splitBy, split = "/")), by = id]
# now split multiple ec numbers for a single reaction which have been concatenated by a "/"
reaEc <- sapply(rxns[,reaEc], function (x) splitBy(x, split = "/")) # split all ec entries
reaEc[sapply(reaEc, rlang::is_empty)] <- "" # add an empty string if reaEc is empty
lengthEc <- sapply(reaEc, length) # count the number of entries (splits)
sel <- which(lengthEc != 1) # select those which are not 1
lengthEc <- lengthEc[sel] # reduce to the selection
reaEc <- reaEc[sel] # reduce to the selection

rxns2 <- rxns[sel,]
pathwayId <-sapply(1:length(rxns2[,id]), function(x) rep(rxns2[x,id],lengthEc[x]))
reaId <-sapply(1:length(rxns2[,reaId]), function(x) rep(rxns2[x,reaId],lengthEc[x]))
reaName <-sapply(1:length(rxns2[,reaName]), function(x) rep(rxns2[x,reaName],lengthEc[x]))
rxns2 <- data.table(id = unlist(pathwayId),
                   reaEc = unlist(reaEc),
                   reaId = unlist(reaId),
                   reaName = unlist(reaName))
# combine the extended table with the old table
rxns <- rbind(rxns[-sel,], rxns2)

# add the hash for the reaName
rxns[,reaMD5 := sapply(reaName, digest, serialize = FALSE)]

#fwrite(pathways, file = "test.csv")
#rxnNamd <- "exo-chitinase (non-reducing end)"
#a <- digest(rxnNamd, serialize = F) # this is it! column in pathways is reaName
#b <- system2("echo", paste("-n",rxnNamd,"| md5sum | cut -f1 -d' ' "), stdout =TRUE)
#print(paste(a,b))
#a==b

# load diamond result
dmnd <- fread(diamond_out)
colnames(dmnd) <-  dmnd_colnames

# split the fasta header again
split_fasta_header <- function(header, element = 1){
    header_res <- strsplit(header, split ="\\|")[[1]][1]
    header_res <- strsplit(header_res, split = "_")[[1]][element]
    return(header_res)
}
dmnd[,c("seq_source", "enzyme") := list(split_fasta_header(sseqid, 1), split_fasta_header(sseqid,2)), by = sseqid]
dmnd[,sseqid := gsub("^\\w+_.*\\|(.*)$","\\1", sseqid)]

# apply the seqSrc filter
if(seqSrc == 1){
    dmnd <- dmnd[seq_source == "rev",]
} else if(seqSrc == 2){
    enz_rev <- unique(dmnd[seq_source == "rev",enzyme])
    enz_unrev <- unique(dmnd[seq_source == "unrev", enzyme])
    enz_unrev <-setdiff(enz_unrev, enz_rev)
    dmmd <- dmnd[seq_source == "rev"| enzyme %in% enz_unrev,] 
} 







print(args)

print("Not implemented yet")
