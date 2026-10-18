###### -- MGnifyR workflow ----------------------------------------------------

library(MGnifyR)
options(timeout = 600)
# use a cache in R's temp directory
mg <- MgnifyClient(useCache = TRUE)

# danish waste water collection study from copenhagen
study_id <- "MGYS00001135"

analysis_ids <- searchAnalysis(mg, "studies", study_id)

# if you haven't changed R's default timeout option, this may complain
analysis_meta_data <- getMetadata(x = mg,
                                  accession = analysis_ids)

# this is slightly, idiomatic without specifying `FUN`, tapply will just return integer
# ids for the shared factor combinations, in this case, latitute and longitude
# it can contain NAs when NAs are present in the index
sample_locations <- tapply(X = analysis_meta_data,
                           INDEX = list(analysis_meta_data$`sample_geographic location (longitude)`,
                                        analysis_meta_data$`sample_geographic location (latitude)`))

# any(is.na(analysis_meta_data$`sample_geographic location (latitude)`))
# any(is.na(analysis_meta_data$`sample_geographic location (longitude)`))

# we want to select the sampling location with the fewest samples, for simplicity
# table the sample location index, find the lowest value, get the name
# that in this case is just an integer
w1 <- as.integer(names(which.min(table(sample_locations))))

# which() manages NAs for us
analysis_subset <- analysis_meta_data[which(sample_locations == w1), ]

# getting a single tree summarized experiment
tse <- getResult(x = mg,
                 accession = analysis_subset[["analysis_accession"]][1],
                 get.func = FALSE)

# getting a multiassayexperiment containing the tses for several accessions
mse <- getResult(x = mg,
                 accession = analysis_subset[["analysis_accession"]])

save(tse,
     mse,
     file = "~/Repos/MetagenomicsWithBioconductor/inst/extdata/workshop_data.RData",
     compress = "xz")
