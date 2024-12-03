# Porthmeus
# 25.10.24

# this script takes a output from a diamond blast and reformats it into the expected form of the pathway and reactions tables of gapseq
print(args)

# libraries
library(data.table)
library(digest)

# arguments and variables
args <- commandArgs(trailingOnly = TRUE)
diamond_out <- args[1] # path to the diamond output file
taxonomy <- args[2] # included taxonomy
taxRange <- args[3] # included taxonomyRange TODO needs to be implemented
seqSrc <- args[4] # which sequences should be included user, rev, unrev?
srcDir <- args[5] # directory of the script to find the data
noSuperpathways <- args[6] # boolean, should super pathways be included
identcutoff <- args[7] # cutoff for the identity in the blast to be considered 'good'
identcutoff_exception <- args[8] # some enzymes require higher identcutoff, because of high similarity to other enzymes, how high the cutoff must be is defined here
database <- args[9] # vmh or seed?
subunit_cutoff <- args[10] # how much of subunit needs to be found to be considered present in the sequences? default 0.5
bitcutoff <- args[11] # cutoff on the bitscore for good blasts
covcutoff <- args[12] # coverage cutoff - how much needs to be covered of the enzyme in the blast
completenessCutoff <- args[13] # cutoff for the completeness of the pathway to be considered present, if key reactions are found (default 0.66)
completenessCutoffNoHints <- args[14] # cutoff for the completeness of the pathway to be considered present if no key reactions are found (default 0.8)
pwyKey <- "Pathways|Enzyme-Test|seed|kegg" # corresponds to -p all


dmnd_colnames <- c("qseqid","pident","evalue","bitscore","scovhsp","sseqid","sstart","send")

#### debug ####
#diamond_out <- "./ecoli.fna_prodigal.fasta_blastresult.tsv"
#taxonomy <- "Bacteria"
#taxRange <- "all"
#srcDir <-"./src"
#noSuperpathways <- "true"
#seqSrc <- 2
#identcutoff <- 0
#identcutoff_exception <- 70
#database <- "seed"
#subunit_cutoff <- 0.5
#bitcutoff <- 200 # bitscore cutoff
#covcutoff <- 75 # 
#completenessCutoff <- 66
#completenessCutoffNoHints <- 80


# code
# load pathway databases
pwys <- c("meta","kegg","seed","custom")
pwfls <- file.path(srcDir, "..","dat", paste0(pwys,"_pwy.tbl"))
names(pwfls)<- pwys
pathways <- lapply(names(pwfls),function(x){ cbind(source = x, fread(pwfls[x], fill = F))})
pathways[["fill"]] <- TRUE
pathways <- do.call(rbind, pathways)

# remove duplicated pathways
dup_ids <- pathways[duplicated(id),id]
pathways <-pathways[!(id %in% dup_ids) | source == "custom",]

# remove pathways without reaction
pathways <- pathways[!is.na(reaId),]

# remove super-pathways if set
if(noSuperpathways == "true"){
    pathways <- pathways[!grepl("Super-Pathways", hierarchy),]
}

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

#fwrite(pathways, file = "test.csv")
#rxnNamd <- "exo-chitinase (non-reducing end)"
#a <- digest(rxnNamd, serialize = F) # this is it! column in pathways is reaName
#b <- system2("echo", paste("-n",rxnNamd,"| md5sum | cut -f1 -d' ' "), stdout =TRUE)
#print(paste(a,b))
#a==b
# add the hash for the reaName
rxns[,reaMD5 := sapply(reaName, digest, serialize = FALSE)]

# add information about spontaneous and key reactions in pathways
# key reactions
keyReactions <- sapply(pathways[,keyRea], splitBy,split = ",")
lengthKeyRea <- sapply(keyReactions, length)
pathwayId <-sapply(1:length(pathways[,id]), function(x) rep(pathways[x,id],lengthKeyRea[x]))
keyReactions <- data.table(id = unlist(pathwayId),
                           keyRea = unlist(keyReactions))
keyReactions <-keyReactions[!(is.na(keyRea)),]

if(!all(keyReactions[,keyRea] %in% rxns[,reaId])){
    sel <- !(keyReactions[,keyRea] %in% rxns[,reaId])
    out <- capture.output(print(keyReactions[sel,]))
    out <- paste(out, collapse = "\n")
    warning(paste0("The following key reactions are not part of the current pathway reactions:\n\n",out ))
}

# add the keyReaction indicator
rxns[paste(id,reaId, sep = "__") %in% keyReactions[,paste(id,keyRea, sep= "__")], keyRea := TRUE]
rxns[is.na(keyRea), keyRea := FALSE]


# spontaneous reactions
spontRea <- sapply(pathways[,spont], splitBy, split = ",")
lengthSpontRea <- sapply(spontRea, length)
pathwayId <- sapply(1:length(pathways[,id]), function(x) rep(pathways[x,id],lengthSpontRea[x]))
spontReactions <- data.table(id = unlist(pathwayId),
                           spontRea = unlist(spontRea))
spontReactions <- spontReactions[!is.na(spontRea),]

if(!all(spontReactions[,spontRea] %in% rxns[,reaId])){
    sel <- !(spontReactions[,spontRea] %in% rxns[,reaId])
    out <- capture.output(print(spontReactions[sel,]))
    out <- paste(out, collapse = "\n")
    warning(paste0("The following spontaneous reactions are not part of the current pathway reactions:\n\n",out ))
}

rxns[paste(id,reaId, sep = "__") %in% spontReactions[,paste(id,spontRea, sep= "__")], spontRea := TRUE]
rxns[is.na(spontRea), spontRea := FALSE]

# add reactions which need an adjustment for the bitscore
rxns[,identcutoff := identcutoff]
if(identcutoff_exception > identcutoff){
    exceptions.tbl <- fread(file.path(srcDir, "..","dat","exception.tbl"), header = TRUE)
    rxns[reaName %in% exceptions.tbl[[1]] |
         reaId %in% exceptions.tbl[[1]]   | 
         reaEc %in% exceptions.tbl[[1]],
         identcutoff:= identcutoff_exception]
}

# add the alternative EC numbers
altec <- scan(file.path(srcDir, "..","dat","altec.csv"), what = "text")
altec <- lapply(altec, splitBy, split = ",")
altec_num <- lapply(1:length(altec), function(x) rep(x, length(altec[[x]])))
altec <- data.table(cluster = paste0("alt",unlist(altec_num)),
                    ec = unlist(altec))

# get also alternatives from brendaec, remove all brendaecs where alread alternatives in the altec exists
brendaec <- fread(file.path(srcDir, "..","dat","brenda_ec_edited.csv"))
brendaec[,ec2 := gsub(".* (\\d*\\.\\d*\\.\\d*\\.\\d*).*","\\1",History)]
brendaec <- brendaec[grepl("\\d*\\.\\d*\\.\\d*\\.\\d*",ec2)]
brendaec <- unlist(lapply(1:nrow(brendaec), function(x) c(brendaec[x,`EC Number`],brendaec[x,ec2])))
brendaec <- data.table(cluster = paste0("brenda",sort(rep(1:(length(brendaec)/2),2))),
                       ec = brendaec)

sel <- brendaec[(ec %in% altec[,ec]),cluster]
brendaec <- brendaec[!(cluster %in% sel),]
# combine both tables
altec <- rbind(altec,brendaec)
altec <- merge(altec,altec, by = "cluster",allow.cartesian=TRUE)
# remove alternative ecs with unknown acceptors
altec <-altec[!(grepl("(\\.99\\.\\d*$)",ec.x)|grepl("(\\.99\\.\\d*$)",ec.y)),]

# add alternative ecs to the rxns
rxns2 <- merge(rxns, altec[,.(ec.x,ec_alt = ec.y)], by.x = "reaEc",by.y = "ec.x", all.x = TRUE, allow.cartesian=TRUE)
rxns2[!is.na(ec_alt) & ec_alt != reaEc,reaEc := ec_alt]
rxns2[,ec_alt := NULL]
rxns <- rxns2



# after this preparation, we can create the results tables for reactions and pathways
##### debug #####
# rxn.tbl <- fread("ecore-min-Reactions.tbl")
# pwy.tbl <- fread("ecore-min-Pathways.tbl")

# load subunit table, which describe the subunits formation of a reaction
subunits <- fread(file.path(srcDir, "..","dat","seq",taxonomy,"Subunits.csv"))
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

# add the subunit information
dmnd <- merge(dmnd,
              subunits[,.(directory,
                    file = gsub("\\.fasta$","", file),
                    seqID,
                    subunit,
                    max_subunit)],
              by.x = c("seq_source","enzyme","sseqid"),
              by.y = c("directory","file","seqID"), all.x=TRUE)


# split into three tables, one for identifier derived by EC, one for md5sum and one for rxnID
dmnd.ec <- dmnd[enzyme %in% rxns[,unique(reaEc),],]
dmnd.md5 <- dmnd[enzyme %in% rxns[,unique(reaMD5),],]
dmnd.rxn <- dmnd[enzyme %in% rxns[,unique(reaId),],]

# merge the table to the rxns
rxn.dmnd.ec <- merge(dmnd.ec, rxns, by.x = "enzyme", by.y = "reaEc", allow.cartesian=TRUE)
rxn.dmnd.ec[,reaEc := enzyme]
rxn.dmnd.md5 <- merge(dmnd.md5, rxns, by.x = "enzyme", by.y = "reaMD5", allow.cartesian=TRUE)
rxn.dmnd.md5[,reaMD5 := enzyme]
rxn.dmnd.rxn <- merge(dmnd.rxn, rxns, by.x = "enzyme", by.y = "reaId", allow.cartesian=TRUE)
rxn.dmnd.rxn[,reaId := enzyme]
# merge the tables again
rxn.dmnd <- rbind(rxn.dmnd.ec, rxn.dmnd.md5, rxn.dmnd.rxn)

# check for good blasts
rxn.dmnd[,status := "bad_blast"]
rxn.dmnd[pident > identcutoff & bitscore > bitcutoff & scovhsp > covcutoff, status := "good_blast"]




print(args)

print("Not implemented yet")
