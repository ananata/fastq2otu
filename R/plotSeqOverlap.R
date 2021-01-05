#' Plot the quality profile of paired-end data
#' The figure created can be used to visually determine the best merging and trimming parameters.
#'
#' @param tabF (Required). Forward table obtained from the getMeanQualityProfile() function. 
#' 
#' @param tabR (Required). Reverse table obtained from the getMeanQualityProfile() function. 
#'
#' @param primerF (Required). Forward primer sequence provided in the 5' to 3' orientation (must NOT feature any degenerate bases or invalid
#' characters).
#'
#' @param primerR (Required). Reverse primer sequence provided in the 5' to 3' orientation (must NOT feature any degenerate bases or invalid
#' characters).
#'
#' @param projectLabel (Required). Default is "Paired-end Study". 
#' Character label (i.e. project accession number) that can be used to describe data being analyzed.
#' 
#' @param matchScore (Default is 1). Match alignment score.
#'
#' @param mismatchScore (Default is -3). Mismatch alignment score.
#'
#' @param gapOpeningScore (Default is 5). Gap opening score.
#'
#' @param gapExtensionScore (Default is 2). Gap extension score.
#'
#' @param minScore (Default is 10). Minimum alignment score to accept.
#'
#' @import Biostrings
#' @import ggplot2
#' @return ggplot figure
#'
#' @export
#'
plotSeqOverlap <- function(tabF, tabR, primerF, primerR, projectLabel = "Paired-end Study", matchScore = 1, mismatchScore = -3,
                           gapOpeningScore = 5, gapExtensionScore = 2, minScore = 4) {
  
  ## E. coli reference sequence (NCBI Accession: J01859.1)
  ## TODO: Export to individual object
  reference <- "AAATTGAAGAGTTTGATCATGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGT
  AACAGGAAGAAGCTTGCTCTTTGCTGACGAGTGGCGGACGGGTGAGTAATGTCTGGGAAACTGCCTGATG
  GAGGGGGATAACTACTGGAAACGGTAGCTAATACCGCATAACGTCGCAAGACCAAAGAGGGGGACCTTCG
  GGCCTCTTGCCATCGGATGTGCCCAGATGGGATTAGCTAGTAGGTGGGGTAACGGCTCACCTAGGCGACG
  ATCCCTAGCTGGTCTGAGAGGATGACCAGCCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGG
  CAGCAGTGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTATGAAGAAGGCCTT
  CGGGTTGTAAAGTACTTTCAGCGGGGAGGAAGGGAGTAAAGTTAATACCTTTGCTCATTGACGTTACCCG
  CAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTAATCGGAAT
  TACTGGGCGTAAAGCGCACGCAGGCGGTTTGTTAAGTCAGATGTGAAATCCCCGGGCTCAACCTGGGAAC
  TGCATCTGATACTGGCAAGCTTGAGTCTCGTAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGCGT
  AGAGATCTGGAGGAATACCGGTGGCGAAGGCGGCCCCCTGGACGAAGACTGACGCTCAGGTGCGAAAGCG
  TGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGTCGACTTGGAGGTTGTGCCC
  TTGAGGCGTGGCTTCCGGAGCTAACGCGTTAAGTCGACCGCCTGGGGAGTACGGCCGCAAGGTTAAAACT
  CAAATGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGATGCAACGCGAAGAACCT
  TACCTGGTCTTGACATCCACGGAAGTTTTCAGAGATGAGAATGTGCCTTCGGGAACCGTGAGACAGGTGC
  TGCATGGCTGTCGTCAGCTCGTGTTGTGAAATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTATCCT
  TTGTTGCCAGCGGTCCGGCCGGGAACTCAAAGGAGACTGCCAGTGATAAACTGGAGGAAGGTGGGGATGA
  CGTCAAGTCATCATGGCCCTTACGACCAGGGCTACACACGTGCTACAATGGCGCATACAAAGAGAAGCGA
  CCTCGCGAGAGCAAGCGGACCTCATAAAGTGCGTCGTAGTCCGGATTGGAGTCTGCAACTCGACTCCATG
  AAGTCGGAATCGCTAGTAATCGTGGATCAGAATGCCACGGTGAATACGTTCCCGGGCCTTGTACACACCG
  CCCGTCACACCATGGGAGTGGGTTGCAAAAGAAGTAGGTAGCTTAACCTTCGGGAGGGCGCTTACCACTT
  TGTGATTCATGACTGGGGTGAAGTCGTAACAAGGTAACCGTAGGGGAACCTGCGGTTGGATCACCTCCTT
  A"
  
  # Remove all unwanted whitespace
  clean_ref <- gsub("[[:space:]]", "", reference)
  clean_primerF <- gsub("[[:space:]]", "", primerF)
  clean_primerR <- gsub("[[:space:]]", "", primerR)
  
  ecoli_seq <- DNAString(clean_ref)
  primerF_seq <- DNAString(clean_primerF)
  primerR_seq <- reverseComplement(DNAString(clean_primerR)) # Take the reverse complement of the sequence
  
  # First use a fixed substitution matrix
  mat <- nucleotideSubstitutionMatrix(match = matchScore, mismatch = mismatchScore, baseOnly = TRUE)
  
  # Align the forward primer to the ecoli rRNA gene
  alignF <- pairwiseAlignment(ecoli_seq, primerF_seq, type = "local", substitutionMatrix = mat,
                              gapOpening = gapOpeningScore, gapExtension = gapExtensionScore)
  if (alignF@score < minScore) {
    stop("Primer could not be detected")
  }
  
  # Get forward sequence starting position
  align_startF <- as.numeric(as.character(alignF@pattern@range@start))
  seqStartF <- align_startF + length(primerF)
  
  # Align the reverse primer to the ecoli rRNA gene
  alignR <- pairwiseAlignment(ecoli_seq, primerR_seq, type = "local", substitutionMatrix = mat,
                              gapOpening = gapOpeningScore, gapExtension = gapExtensionScore)
  if (alignR@score < minScore) {
    stop("Primer could not be detected")
  }
  
  # Get reverse sequence starting position 
  align_startR <- as.numeric(as.character(alignR@pattern@range@start))
  seqStartR <- align_startR + length(primerR)

  # Prepare graph
  COLORS <- c("Forward" = "#00AFBB", "Reverse" = "#E7B800")
  overlap <- ggplot(data=tabF, aes(x=Position + seqStartF, y=Mean)) +
                geom_line(color = COLORS["Forward"], size = 1) +
                geom_area(fill = COLORS["Forward"], alpha = 0.4) +
                geom_line(data = tabR, aes(x = rev(Position + seqStartR), y = Score), color = COLORS["Reverse"], size = 1) +
                geom_area(data = tabR, aes(x = rev(Position + seqStartR), y = Score), fill = COLORS["Reverse"], alpha = 0.4) +
                labs(
                  title = projectLabel,
                  subtitle = "Aggregated Quality Distribution of Unmerged Reads",
                  color = "Legend"
                ) +
                scale_color_manual(values = COLORS) + 
                xlab("Position (bp)") + 
                ylab("Quality Score (Phred)") +
                theme_bw(base_size = 25) + coord_fixed(ratio = 3) +
                theme(plot.title = element_text(face="bold")) +
                scale_x_continuous(minor_breaks = seq(
                    (signif(max(rev(tabR$Position + seqStartR)), 1)/500) * 10, 
                    signif(max(rev(tabR$Position + seqStartR)), 1), 
                    (signif(max(rev(tabR$Position + seqStartR)), 1)/500) * 10),
                    limits = c(0, signif(max(rev(tabR$Position + seqStartR) + 100), 1))
                    )
  # Save plot 
  ggsave(paste0(projectLabel, "_unmerged_paired_reads.pdf"),
         plot = overlap, 
         device = "pdf")
  message("Created: ", paste0(projectLabel, "_unmerged_paired_reads.pdf"))
  message("Forward Start Position: ", seqStartF)
  message("Reverse Start Position: ", seqStartR)

  # Return plot
  return(overlap)
}
