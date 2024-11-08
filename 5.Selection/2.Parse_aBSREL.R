# Parses multiple json files and returns a table of all omega values and percent of sites under selection for each rate class
# must be used in directory with all the json files
# cd /scratch/csm6hg/codeml/mega_OG
# mkdir jsonOutput
# cp OG*/*json jsonOutput

# module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1; R

# Libraries
library(data.table)
library(tidyverse)
library(jsonlite)
library(foreach)

# Working directory 
setwd("/scratch/csm6hg/codeml/mega_OG/jsonOutput/")

# remember to replace Inf values with "Inf" using
# sed -i 's/inf/"inf"/'

# Files
files <- list.files(full.names = TRUE, pattern = '*.json$')
print(files)
string <- "p2n.ABSREL.json"

# Go through each file and extract information - Only significant < 0.05 adjusted p-value
foreach(file = 1:length(files), .errorhandling = "remove") %do% {
  #file=2
  
  # Which file we are working on
  print(files[file])
  
  # Extracting results
  global_results <- data.frame(name = character(0), branch = character(0), 
                               omega1 = character(0), percent1 = character(0), omega2 = character(0),
                               percent2 = character(0), omega3 = character(0), percent3 = character(0),
                               stringsAsFactors = FALSE)
  
  name <- strsplit(paste(files[file]), split = string)
  simplejson <- fromJSON(files[file])
  index <- 1
  branchnames <- list()
  branchnames <- simplejson$`branch attributes`$`0`
  
  # Forloop for each branch - extract values
  for(i in branchnames){
    branch <- names(branchnames[index])
    local_results <- data.frame(name = character(0), branch = character(0), omega1 = character(0), 
                                percent1 = character(0), omega2 = character(0),
                                percent2 = character(0), omega3 = character(0), 
                                percent3 = character(0), pval = character(0), stringsAsFactors = FALSE)
    local_results[1, 1] <- name
    local_results[1, 2] <- branch
    local_results[1, 9] <- i$`Corrected P-value`
    if(i$`Rate classes` == 1){
      local_results[1, 3] <- i$`Rate Distributions`[1,1]
      local_results[1, 4] <- i$`Rate Distributions`[1,2]
      local_results[1, 5] <- "NA"
      local_results[1, 6] <- "NA"
      local_results[1, 7] <- "NA"
      local_results[1, 8] <- "NA"
    } else if (i$`Rate classes` == 2){
      local_results[1, 3] <- i$`Rate Distributions`[1,1]
      local_results[1, 4] <- i$`Rate Distributions`[1,2]
      local_results[1, 5] <- i$`Rate Distributions`[2,1]
      local_results[1, 6] <- i$`Rate Distributions`[2,2]
      local_results[1, 7] <- "NA"
      local_results[1, 8] <- "NA"
    } else {
      local_results[1, 3] <- i$`Rate Distributions`[1,1]
      local_results[1, 4] <- i$`Rate Distributions`[1,2]
      local_results[1, 5] <- i$`Rate Distributions`[2,1]
      local_results[1, 6] <- i$`Rate Distributions`[2,2]
      local_results[1, 7] <- i$`Rate Distributions`[3,1]
      local_results[1, 8] <- i$`Rate Distributions`[3,2]
    }
    global_results <- dplyr::bind_rows(global_results, local_results)
    rm(local_results)
    index <- index + 1
    
  }
  
  if (nrow(global_results) > 0) {
    write.table(global_results, file = "parsed_absrel_json.txt", sep = "\t", 
                row.names = FALSE, col.names = FALSE, append = TRUE, quote = FALSE)
  }
}
