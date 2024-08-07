## Code to produce UK weights file based on UK landcover (what size neightbourhood)
  # Load in test data from Oli's original markdown document
    s <- read.csv(file = "../Background files/frescaloR-master/clusterTestDat.csv", stringsAsFactors = F, header = T)
  # Keep only columns we need/want
    s <- s[,c(4,2,3)]
  # Now save as package dataset (ovewriting if necessary)
    usethis::use_data(s, overwrite = TRUE)


## Code to produce test dataset unicorns
  # Load in weights file from Oli's original markdown document
    d <- read.delim(file = "../Background files/frescaloR-master/GB_LC_Wts.txt", header = F, sep = "")
  # keep only columns we need/want
    d <- d[,c("V1","V2","V3")]
  # Rename the columns
    names(d) <- c("location1", "location2", "w") # perhaps focal and neighb as names rather than location1 and location2?
  # Now save as package dataset (ovewriting if necessary)
    usethis::use_data(d, overwrite = TRUE)

