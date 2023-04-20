##
library(tidyverse)

### Explanation of variables
# uR - transport velocity of mRNA granule, micrometers / second;
# DR - diffusion coefficient for mRNA, micrometers^2 / second;
# kR - degradation rate of mRNA;
# betaR - transcription rate of mRNA;
# lifTR - half-life of a transcript, hours;
# The function is invariant with respect to transcription rate
input = list(name = "yada",
             tr_hl = 4.94,
             nbins = 100,
             ddist = 300,
             DR = 3.0 * 10^(-3),
             uR = 5.8 * 10^(-3),
             betaR = 0.0109)
# ------
## Function to compute mRNA density over distance given input parameters
rdd <- function(input) {
  
  # -----
  # Compute intermediate variables
  LifTR = input$tr_hl * 60 * 60
  kR = log(2) / LifTR
  lambdaR = (sqrt(input$uR^2 + 4 * kR * input$DR) - input$uR) / (2 * input$DR)
  
  # -----
  # Compute distance grid
  ddgrid <- seq(1, input$ddist, length = input$nbins)
  
  # -----
  # Estimate mRNA probability density
  rdgrid <- purrr::map(ddgrid, 
                       ~ ((input$betaR * lambdaR) / kR) * exp(-lambdaR * .x)) %>%
    flatten_dbl()
  
  # -----
  # Standartize values
  rdgrid_st <- purrr::map(rdgrid[2:length(rdgrid)], 
                          ~ .x / rdgrid[[1]]) %>% 
    flatten_dbl() %>% 
    c(1, .)
  
  # -----
  # Prepare tibble output
  tib_out <- tibble(group = input$name, 
                    rna_hl = input$tr_hl,
                    dend_dist = ddgrid,
                    rna_pdens = rdgrid_st)
  return(tib_out)
}
