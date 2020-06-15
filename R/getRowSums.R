
getRowSums <- function(rds.path) {
  table <- readRDS(rds.path)
  sums <- as.numeric(rowSums(table))
  return(sums)
}

