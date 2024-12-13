#!/usr/bin/Rscript
# Porthmeus
# 07.11.24

library(data.table)
library(parallel)

# first argument: input file
# second argument: output file

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  stop("First argument should be an input directory and the second one the output file.")
}


tax_dic <- args[1]
# tax_dic <- "Bacteria"
#tax_dic <- file.path("..","dat","seq",tax_dic)

dirs <- list.files(tax_dic)
fls <- lapply(dirs, function(x) list.files(file.path(tax_dic,x), pattern = "*.fasta"))
names(fls) <- dirs

subunit_pattern <- 'subunit|chain|polypeptide|component'

isRomanWord <- function(word,minRomans = 0.6){
    # simple function to determine if the word is a roman number
    # @word is the word to check
    # minRomans is the fraction of roman numbers in the word to classify as roman
    # returns TRUE/FALSE
    romans <- ("[I|V|X]")
    romans <- c("I","V","X")
    wordsplit <- strsplit(word, "")[[1]]
    frac_romans <- sum(toupper(wordsplit) %in% romans)/length(wordsplit)
    return(frac_romans > minRomans)
}


isMostlyRoman <- function(headers){
    
    headers <- gsub(paste0(".(",subunit_pattern,").") ," \\1 ", headers, ignore.case=TRUE) # remove possible hyphens or similar to combine subunits with numbers (e.g. subunit-gamma -> subunit gamma)
    headers <- strsplit(headers, split = " +") # split the header by spaces
    numberIndex <- lapply(headers, function(x) grep( pattern =subunit_pattern,x, ignore.case=TRUE) +c(-1,1)) # get the words index in each header of the file which potentially denote the number of the subunit
    roman_header <- sapply(1:length(headers),function(x) any(c(
                                            isRomanWord(headers[[x]][numberIndex[[x]][1]]),
                                            isRomanWord(headers[[x]][numberIndex[[x]][2]]))))
    return(sum(roman_header)/length(roman_header) > 0.5)
}
    
romanToArabic <- function(words){
    arabic <- 1:15
    names(arabic) <- c("I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII", "XIII", "XIV", "XV")
    arabic <- rev(arabic)
    for(i in 1:length(words)){
        matches <- names(arabic) == toupper(words[i])
        if(any(matches)){
            words[i] <- arabic[matches]
        }
    }
    return(words)
}

toNumeral <- function(words){
    # create patterns to search for
    greek <- 1:17
    names(greek) <- c("alpha","beta","gamma","delta","epsilon","zeta","eta","theta","iota","kappa","lambda","my","ny","omikron","pi","rho","sigma")
    largeSmall <- 1:2
    names(largeSmall) <- c("large","small")
    latin <- 1:26
    names(latin) <- letters
    patterns <- list(greek, latin, largeSmall)
#
    # look only for exact matches
    for(pat in patterns){
        for(i in 1:length(words)){
            matches <- names(pat) == tolower(words[i])
            if(any(matches)){
                words[i] <- pat[matches]
            }
        }
    }
    return(words)
}
    
extractPrePostWord <- function(header){
    header <- gsub(paste0(".(",subunit_pattern,").") ," \\1 ", header, ignore.case=TRUE) # remove possible hyphens or similar to combine subunits with numbers (e.g. subunit-gamma -> subunit gamma)
    header <- strsplit(header, split = " +")[[1]] # split the header by spaces
    indx<-grep( pattern =subunit_pattern, header, ignore.case=TRUE)+c(-1,1)
    return(header[indx])
}

getNumber <- function(words){
    if(!is.na(words[2])){
        return(words[2])
    } else if(!is.na(words[1])){
        return(words[1])
    } else {
        return(NA)
    }
}
        

# go through the files and try to parallelize
threads <-  detectCores()
subunits_tab <- data.table()
for(dr in dirs){
    drfls <- fls[[dr]]
    cl <- makeForkCluster(threads)
    sub_tabs <- parLapply(cl = cl, drfls, function(fl){
   # for(fl in fls[[dr]]){
        if(!is.null(fl)){
            grp.fl <- file.path(tax_dic, dr, fl)
            arguments <- paste("-E 'subunit|chain|polypeptide|component'", gsub('([;|"])',"\\\\\\1",grp.fl)) # remove some naming artifacts
            header <- suppressWarnings(system2("grep",arguments, stdout =TRUE))
            if(!rlang::is_empty(header)){
                prepost <- lapply(header, extractPrePostWord)
                if(isMostlyRoman(header)){
                    prepost <- lapply(prepost, romanToArabic)
                }
            prepost <- lapply(prepost, toNumeral)           
            prepost <- suppressWarnings(lapply(prepost,as.numeric))
            subunits <- sapply(prepost, getNumber)
           # 
            # save the subunit information in a table
            seqID <- sapply(header, function(x) gsub(">", "", strsplit(x,split = " ")[[1]][1]),USE.NAMES = FALSE)
            # get the number of subunits (count the different subunits instead of taking the max) -NA
            max.subunits <- length(unique(subunits))# - any(is.na(subunits)) # -any(...) was the idea to not penalize blast due to undefined subunits, yet original code does not check that either
            sub_tab <- data.table(directory = dr,
                                  file = fl,
                                  seqID = seqID,
                                  header = header,
                                  subunit = subunits,
                                  max_subunit = max.subunits)
            return(sub_tab)
            }
        }
    })
    stopCluster(cl)
    subunits_tab <- rbind(subunits_tab, do.call(rbind, sub_tabs))
}

# add a 1 where only NA subunits where detected
#subunits_tab[max_subunit==-Inf, max_subunit := 1]
fwrite(subunits_tab, args[2])

