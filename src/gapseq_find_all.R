# Porthmeus
# 25.10.24


# this script takes a output from a diamond blast and reformats it into the expected form of the pathway and reactions tables of gapseq
#print(args)

# libraries
library(data.table)
library(digest)
library(stringr)

# arguments and variables
args <- commandArgs(trailingOnly = TRUE)
diamond_out <- args[1] # path to the diamond output file
taxonomy <- args[2] # included taxonomy
taxrange <- args[3] # included taxonomyRange, this is used to pre-filter pathways which should be searched
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
vagueCutoff <- args[15] # cutoff maximal allowed of vague reactions in a pathway to be predicted (vague reactions do not have a reference sequence in the database)
strictCandidates <- args[16] # consider vague reactions as part of the pathway (strictCandidates == false) or not (strictCandidates == true)
outfile.rxns <- args[17] # output file for the reactions
outfile.pathways <- args[18] # output file for the patways
pwyKey <- "Pathways|Enzyme-Test|seed|kegg" # corresponds to -p all




dmnd_colnames <- c("qseqid","pident","evalue","bitscore","scovhsp","sseqid","sstart","send","qstart","qend")

#### debug ####
#diamond_out <- "./ecore.faa_blastresult.tsv"
#taxonomy <- "Bacteria"
#taxrange <- "all"
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
#vagueCutoff <- 0.3
#strictCandidates <- "false"
#outfile.rxns<- "ecore.faa-dbg-Reactions.tbl"
#outfile.pathways<- "ecore.faa-dbg-Pathway.tbl"


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

# remove pathways outside the taxrange
if(taxrange != "all"){
    # get the actual taxrange
    taxonomies <- fread(file.path(srcDir, "..","dat","taxonomy.tbl"))
    taxonomies.character <- apply(taxonomies, 1, paste, collapse = "")
    sel <- grep(taxrange, taxonomies.character, ignore.case = TRUE)
    taxranges <- paste0("\\|TAX-",taxonomies[sel, tax], "\\|")
    sel <- unique(unlist(lapply(taxranges, grep, x = pathways[,taxrange])))
    pathways <- pathways[sel,]
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

# add information for vague reactions 
if(seqSrc == 1){ # check for unrev only if option is set
    checkdirs <- c("rev","user","rxn")
} else if(seqSrc == 2){
    checkdirs <- c("rev","user","unrev","rxn")
} 
seqDir <- file.path(srcDir, "..","dat","seq",taxonomy)
files <- lapply( checkdirs, function(x) list.files(file.path(seqDir, x)))
non_vague_rxn <- gsub(".fasta$","",unlist(unique(files)))
rxns[, vague := !(reaId %in% non_vague_rxn|reaMD5 %in% non_vague_rxn | reaEc %in% non_vague_rxn)]


# get a table for metacycid/keggid to target database
# define patterns to look for in the tables
KO_pattern<-"R[0-9][0-9][0-9][0-9][0-9]"
EC_pattern<-"[0-9]+\\.[0-9]+\\.[0-9]+\\.[0-9]+"

expandToTable <- function(ids, matches){
    # small function to save some space
    # takes a vector of ids and a list of matches (e.g. kegg ids) of the same length as ids
    # ids will be expanded to the number of matches in each list element and combined to a data table
    lens <- lapply(matches, length)
    ids.exp <- lapply(1:length(ids), function(x) rep(ids[x],lens[[x]]))
    return(unique(data.table(id = unlist(ids.exp),
                      match = unlist(matches))))
}

# get the different tables for metacyc
# 1) Kegg
meta2kegg <- fread(file.path(srcDir,"..","dat","meta_rea.tbl"))
meta2kegg <- unique(meta2kegg[,.(id = gsub("\\|","", id), kegg)])
keggs <- lapply(meta2kegg[,kegg], splitBy, split =",")
meta2kegg <- expandToTable(meta2kegg[,id],keggs)
# 2) EC
meta2ec <- unique(rxns[grep(EC_pattern,reaEc), .(id = reaId, ec = reaEc)])



if(database == "seed"){
    # load the database where the data is stored
    targetDB <- fread(file.path(srcDir, "..","dat","seed_reactions_corrected.tsv"))
    targetDB <- targetDB[status == "OK",] # filter out reactions which are not ok (whatever that means)
    targetString <- apply(targetDB, 1, paste, collapse = "_;_") # create a string per row, as we do not know in which string the kegg ids are stored (maybe even in flow text)
    seed.keggs <- str_extract_all(targetString, KO_pattern) # extract all kegg identifier
    seed2kegg1 <- expandToTable(targetDB[[1]],seed.keggs)

    # do the same for ec numbers
    seed.ecs <- str_extract_all(targetString, EC_pattern)
    seed2ec1 <- expandToTable(targetDB[[1]],seed.ecs)

    # there is a second table for seed to do this, here we know that the information is stored in 'other'
    targetDB <- fread(file.path(srcDir, "..","dat","mnxref_seed-other.tsv"))
    seed2kegg2 <- targetDB[grep(KO_pattern, other),.(id = seed, kegg = other)]
    seed2ec2 <- targetDB[grep(EC_pattern, other), .(id = seed, ec = gsub("-RXN","", other))]
    # match meta to seed via metanetx
    meta2mnx2target <- targetDB[grep("-RXN$|^RXN-", other), .(metaId = other, tarId = seed)]
    # in this table we can also search for mnx references
    
    # for ec-numbers there is even a third table to look for, so lets do this
    targetDB <- fread(file.path(srcDir, "..","dat","seed_Enzyme_Class_Reactions_Aliases_unique_edited.tsv"))
    seed2ec3 <- targetDB[grep(EC_pattern, `External ID`),.(id = `MS ID`, ec = `External ID`)]
    length(unique(seed2ec3[,ec]))
    seed2ec3 <- seed2ec3[,.(id = unique(unlist(strsplit(id, split = "\\|")))), by = ec][,.(id,ec)]


    # combine both tables and remove duplicates
    target2kegg <- unique(rbind(seed2kegg1[,.(id,kegg = match)],seed2kegg2))
    target2ec <- unique(rbind(seed2ec1[,.(id, ec = match)], seed2ec2, seed2ec3))

    # and another table which matches seed2seed
    seed.rxns <- unique(splitBy(split = ",",pathways[source == "seed" & reaId != "", paste(reaId, collapse =",")]))
    target2seed <- data.table(metaId = seed.rxns,
                              tarId = seed.rxns)
    

} else if(database == "vmh"){
    # similar approach as above, search the string of rows instead of a certain column
    targetDB <- fread(file.path(srcDir,"..","dat","vmh_reactions.tsv"))
    targetString <- apply(targetDB, 1, paste, collapse = "_;_")
    vmh.keggs <- str_extract_all(targetString, KO_pattern)
    target2kegg <- expandToTable(targetDB[[1]], vmh.keggs)
    colnames(target2kegg)[2] <- "kegg"
    vmh.ecs <- str_extract_all(targetString, EC_pattern)
    target2ec <- expandToTable(targetDB[[1]], vmh.ecs)
    colnames(target2ec)[2] <- "ec"

    # match meta to vmh via metanetx
    targetDB <- fread(file.path(srcDir, "..","dat","mnxref_bigg-other.tsv"))
    meta2mnx2target <- targetDB[grep("-RXN$|^RXN-", other), .(metaId = other, tarId = bigg)]

    # match seed to vmh 
    target2seed <- data.table(metaId = c(),
                              tarId= c()) # not implemented yet
    warning("No SEED to VMH conversion implemented yet!")
    
    
    # match meta to vmh via bigg
    #targetDB <- fread(file.path(srcDir, "..","dat", "bigg_reactions.tbl"))
    #targetString <- apply(targetDB, 1, paste, collapse = "_;_")
    #vmh.meta <- str_extract_all(targetString, "^RXN-|-RXN$")
    #target2bigg2meta <- expandToTable(targetDB[[1]], vmh.meta)[,.(metaId = match, tarId = id)]
    # returns an empty list

}
# merge the mnx id to the target id
meta2kegg2target <- merge(meta2kegg[,.(metaId = id, kegg=match)], target2kegg[,.(tarId = id, kegg)], by = "kegg")
meta2ec2target <- merge(meta2ec[,.(metaId = id, ec)], target2ec[,.(tarId = id, ec)], by = "ec", allow.cartesian = TRUE)

meta2target <- unique(rbind(meta2ec2target[,.(metaId, tarId)], # metacyc to target via ec number
                     meta2kegg2target[,.(metaId, tarId)], # metacyc to target via kegg number
                     target2kegg[,.(metaId= kegg, tarId = id)], # kegg to target
                     target2seed[,.(metaId, tarId)], # seed to target
                     meta2mnx2target[,.(metaId, tarId)])) # metacyc to target via metanetx

# after this preparation, we can create the results tables for reactions and pathways
##### debug #####
# rxn.tbl <- fread("test/ecore-min-Reactions.tbl")
# pwy.tbl <- fread("test/ecore-min-Pathways.tbl")
#################

# load subunit table, which describe the subunits formation of a reaction
subunits <- fread(file.path(srcDir, "..","dat","seq",taxonomy,"Subunits.csv"))
#subunits[max_subunit == 0, max_subunit := 1]
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
dmnd[!is.na(subunit),complex := paste0("Subunit ", subunit)] # add subunit numbers
dmnd[is.na(subunit) & !is.na(max_subunit) ,complex := "Subunit undefined"] # add undefined subunits

## keep for each enzyme only the best blast
#dmnd <- dmnd[order(bitscore, decreasing = TRUE),]
#dmnd.max <- dmnd[,.(max_bitscore = max(bitscore)),by = .(enzyme, subunit)]
#dmnd <- merge(dmnd, dmnd.max, by =c("enzyme","subunit"))
#dmnd <- dmnd[bitscore == max_bitscore,]

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

# find subunits
subunits_coverage <- rxn.dmnd[!(is.na(max_subunit)),
                              .(subunits_coverage = length(unique(subunit))/max(max_subunit),
                                good.blast.status = any(status == "good_blast")),
                              by = enzyme]  # note subunits coverage, as well as if at least one good blast is in the list
covered_enzymes <- subunits_coverage[subunits_coverage>= subunits_coverage & good.blast.status == TRUE, enzyme] # every blast hit which has either no subunits, or with sufficient subunits coverage will get a complex.status = 1

rxn.dmnd[is.na(max_subunit), complex.status := 1]
rxn.dmnd[enzyme %in% covered_enzymes, complex.status := 1]

# add the database hits
dmnd2target <- meta2target[metaId %in% rxn.dmnd[,unique(reaId)]]
dmnd2target <- dmnd2target[,.(dbhit = paste(tarId, collapse = " ")), by =.(reaId = metaId)]
rxn.dmnd <- merge(rxn.dmnd, dmnd2target, by = "reaId", all.x =TRUE)

# reduce the findings to the relevant fields
reaction.tbl.dmnd <- rxn.dmnd[,.(rxn = reaId,
                                 name = reaName,
                                 ec = reaEc,
                                 bihit = NA, # this is depricated
                                 qseqid = sseqid, # the blast was done revers to initial implementation, thus to keep compatibility return
                                 pident = pident,
                                 evalue = evalue,
                                 bitscore = bitscore,
                                 qcovs = scovhsp, # again reverse s and q, scovhsp is the equivalent to qcovs in diamond
                                 stitle = qseqid, # did not keep the complete title but just id, maybe something to add at some point
                                 sstart = qstart,
                                 send = qend,
                                 pathway = id,
                                 status = status,
                                 pathway.status = as.character(NA), # not yet determined
                                 dbhit = dbhit,
                                 exception = as.integer(identcutoff==identcutoff_exception),
                                 complex.status = complex.status)]


# now get the pathways which are covered 
if(strictCandidates=="true"){
    pathway.rxns <- rxns[reaId %in% reaction.tbl.dmnd[,rxn] | spontRea == TRUE| vague == TRUE,]
} else {
    pathway.rxns <- rxns[reaId %in% reaction.tbl.dmnd[,rxn] | spontRea == TRUE,]
}
# test which pathways are present
# calculate the number of total found, good blasts, bad blasts and vague reactions
totalNrInPwy <- rxns[id %in% pathway.rxns[,id],.(reaNr=length(unique(reaId))), by=id]
foundNrInPwy <- pathway.rxns[,.(foundNr = length(unique(reaId))),by= id]
badNrInPwy <- pathway.rxns[reaId %in% reaction.tbl.dmnd[status == "bad_blast",rxn], .(badNr = length(unique(reaId))),by= id]
vagueNrInPwy <-pathway.rxns[vague ==TRUE, .(vagueNr = length(unique(reaId))), by = id]

totalFoundPwy <- merge(totalNrInPwy,foundNrInPwy, by = "id")
totalFoundPwy <- merge(totalFoundPwy,badNrInPwy, by = "id", all.x = TRUE)
totalFoundPwy <- merge(totalFoundPwy,vagueNrInPwy, by = "id", all.x = TRUE)

totalFoundPwy[is.na(badNr),badNr := 0]
totalFoundPwy[is.na(vagueNr),vagueNr := 0]
totalFoundPwy[,c("Completeness","Good","Bad","Vague") := list((foundNr)/reaNr,
                                                            (foundNr-badNr)/reaNr,
                                                            badNr/reaNr,
                                                            vagueNr/reaNr)]
# check if a key reaction was found
totalFoundPwy[, keyFound := FALSE]
totalFoundPwy[id %in% pathway.rxns[keyRea == TRUE & reaId %in% reaction.tbl.dmnd[,rxn],id], keyFound := TRUE]

# make the checks - the order matters here, so be careful when changed
totalFoundPwy[,c("Prediction","pathway.status") := list(FALSE,as.character(NA))]
# test for treshold pass for pathways with key reactions
totalFoundPwy[keyFound == TRUE & Completeness*100 >= completenessCutoff & Vague <= vagueCutoff,
              c("Prediction", "pathway.status") := list(TRUE, "keyenzyme")]
# test for completeness above threshold
totalFoundPwy[Completeness*100 >= completenessCutoffNoHints & Vague <= vagueCutoff,
              c("Prediction", "pathway.status") := list(TRUE, "threshold")] 
# test for full completeness
totalFoundPwy[Completeness == 1 & Vague <= vagueCutoff,
              c("Prediction", "pathway.status") := list(TRUE, "full")]

# create the correct output
# get the pathway name
pwy.id2name <- unique(pathways[,.(ID=id, name)])
setkey(pwy.id2name,"ID")
# get the number of key reactions
pwy.id2foundKeys <- pathway.rxns[keyRea == TRUE & reaId %in% reaction.tbl.dmnd[,rxn],.(KeyReactions = .N), by =.(ID = id)]
setkey(pwy.id2foundKeys, "ID")
pwy.id2keys <- rxns[keyRea ==TRUE, .(KeyReactions = .N), by = .(ID = id)]
setkey(pwy.id2keys, "ID")

pathways.out <- totalFoundPwy[,.(ID = id,
                                 name = pwy.id2name[id,name],
                                 Prediction,
                                 Completeness = Completeness*100,
                                 VagueReactions = vagueNr,
                                 KeyReactions = pwy.id2keys[id,KeyReactions],
                                 KeyReactionsFound = pwy.id2foundKeys[id,KeyReactions]
                                 )]
# remove NAs and replace with 0
pathways.out[is.na(VagueReactions), VagueReactions :=0]
pathways.out[is.na(KeyReactions), KeyReactions :=0]
pathways.out[is.na(KeyReactionsFound), KeyReactionsFound :=0]

# add the reactions found in each pathway
pathways.out <- unique(merge(pathways.out,
      pathway.rxns[spontRea == FALSE & vague == FALSE, .(ReactionsFound = paste(unique(reaId), collapse =" ")), by = id],
      by.x = "ID", by.y = "id"))


# now add the reactions without blast support but predicted by the pathway completeness to the reactions table
# first add the pathway status to the list of the of blast reactions
setkey(totalFoundPwy, "id")
reaction.tbl.dmnd[,pathway.status := totalFoundPwy[pathway, pathway.status]]

# now add all the reactions which are not in the original blast result
predicted.rxns <- rxns[id %in% pathways.out[Prediction == TRUE,unique(ID)]]
predicted.rxns <- predicted.rxns[!(paste(reaId,id, sep = "_;_") %in% reaction.tbl.dmnd[, paste(rxn, pathway, sep = "_;_")]),]
predicted.rxns <- merge(predicted.rxns, totalFoundPwy[,.(id, pathway.status)], by = "id", all.x = TRUE)
# add status of the reaction - the order matters here
predicted.rxns[,status := "no_blast"]
predicted.rxns[vague == TRUE, status := "no_seq_data"]
predicted.rxns[spontRea == TRUE, status := "spontaneous"]
# add the target hits to the reactions
predicted.rxns <- merge(predicted.rxns, meta2target[,.(dbhit = paste(tarId, collapse = " ")), by = metaId], by.x = "reaId", by.y = "metaId")

# create the actual table
reaction.tbl.pred <- predicted.rxns[,.(rxn = reaId,
                                 name = reaName,
                                 ec = reaEc,
                                 bihit = NA, # this is depricated
                                 qseqid = "",
                                 pident = "",
                                 evalue = "",
                                 bitscore = "",
                                 qcovs = "", 
                                 stitle = "",
                                 sstart = "",
                                 send = "",
                                 pathway = id,
                                 status,
                                 pathway.status, 
                                 dbhit = dbhit,
                                 exception = as.integer(identcutoff==identcutoff_exception),
                                 complex.status = as.character(NA))]
reaction.tbl.out <- unique(rbind(reaction.tbl.dmnd,reaction.tbl.pred))

# write the files to disc
fwrite(reaction.tbl.out, file = outfile.rxns, sep ="\t")
fwrite(pathways.out, file = outfile.pathways, sep = "\t")
