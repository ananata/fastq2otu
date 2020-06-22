#' Denoise sequence data
#' 
#' @param derep Output(s) of derepSeqs function
#' @param err Output(s) of learnErrors function
#' @return  dada object
#' @export
dadaSeqs <- function(derep, err, object) {
  # Execute the function
  dada_obj <- dada::dada(derep, err=err, BAND_SIZE = object@dadaBandSize, OMEGA_A = object@dadaOmegaA)

  # Return derep object
  return(dada_obj)
}

