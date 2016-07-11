dependancies <- c("dplyr", "tidyr", "magrittr", "sonicLength", "devtools")
null <- sapply(dependancies, function(x){
  suppressPackageStartupMessages(
    try(library(x, character.only = TRUE), silent = TRUE))
})

dependancies_present <- sapply(dependancies, function(package){
  package <- paste0("package:", package)
  logic <- package %in% search()
})

if(FALSE %in% dependancies_present){
  Unloaded_Packages <- data.frame(package=as.character(dependancies),
                                  loaded=dependancies_present)
  write.table(Unloaded_Packages,
              file = "Unloaded_Packages.tsv",
              quote = FALSE,
              row.names = FALSE)
  stop("Load required packages. Check Unloaded_Packages.tsv for missing
       dependancies.")
}else{
  remove(dependancies, dependancies_present)
  message("Required packages loaded.")
}

suppressMessages(
  source_url(
    "https://raw.githubusercontent.com/cnobles/cloneTracker/master/cloneTracker.SOURCE_ME.R"
))

args <- commandArgs(trailingOnly = TRUE)

uniqData <- args[ grep("--uniq", args)+1 ]
multiData <- args[ grep("--multi", args)+1 ]

message(paste0("uniqPath: ", uniqData))
message(paste0("multihitPath: ", multiData))

load(uniqData)
load(multiData)

sites.uniq <- db_to_granges(sites.uniq, keep.additional.columns=TRUE)
sites.multi <- db_to_granges(sites.multi, keep.additional.columns=TRUE)

std.uniq <- normalize_intsite_positions(sites.uniq, gap = 5L)
cond.uniq <- condense_intsites(std.uniq, return.abundance = TRUE, method = "estAbund", replicates = "sampleName")

std.multi <- normalize_intsite_positions(sites.multi, gap = 5L)
norm.multi <- normalize_multihit_clutsers(std.multi)
multi.spec <- split(norm.multi, norm.multi$specimen)

multi.sonic.estimate <- do.call(rbind, lapply(multi.spec, function(multi.spec){
  multi.spec.rep <- split(multi.spec, multi.spec$sampleName)
  sonic.table <- do.call(rbind, lapply(1:length(multi.spec.rep), function(i){
    multi.rep <- multi.spec.rep[[i]]
    multi.clus <- split(multi.rep, multi.rep$stdmultihitID)
    data.frame(
      "id" = Rle(
        values = sapply(multi.clus, function(x) unique(x$stdmultihitID)),
        lengths = sapply(multi.clus, function(x) length(unique(width(x))))),
      "length" = as.numeric(unlist(sapply(multi.clus, function(x) unique(width(x))))),
      "replicates" = i
    )
  }))
  if(length(unique(sonic.table$replicates)) == 1){
    sonic.estimate <- with(sonic.table, estAbund(id, length))
  }else if(length(unique(sonic.table$replicates)) > 1){
    sonic.estimate <- with(sonic.table, estAbund(id, length, replicates))
  }
  data.frame(
    "stdmultihitID" = names(sonic.estimate[[1]]),
    "estAbund" = round(sonic.estimate$theta),
    row.names = NULL
  )
}))

n.uniq.sites <- length(cond.uniq)
message(paste0("Number of unique sites: ", n.uniq.sites))
abund.uniq <- sum(cond.uniq$estAbund)
message(paste0("Number of unique sonicLengths: ", abund.uniq))
n.multihits <- length(unique(norm.multi$stdmultihitID))
message(paste0("Number of multihit clusters: ", n.multihits))
abund.multi <- sum(multi.sonic.estimate$estAbund)
message(paste0("Number of multihit sonicLengths: ", abund.multi))
abundances <- c(cond.uniq$estAbund, multi.sonic.estimate$estAbund)
message(paste0("Ranges of sonicAbundances: ", range(abundances)))
ave.abund <- sum(abund.uniq, abund.multi) / sum(n.uniq.sites, n.multihits)
message(paste0("Average integration site abundance: ", ave.abund))
save.image("abundStats.RData")
q()
