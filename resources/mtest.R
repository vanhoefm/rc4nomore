# Copyright (c) 2015 Mathy Vanhoef
# 
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
# 
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
# CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
library("Rmpfr")

#' Perform goodness-of-fit or independence test using the M-test of Fuchs and Kenett.
#'
#' @examples
#' # Two-way contingency table of two one bit-values. On its own each bit is uniformly
#' # distributed, but there is a clear correlation between both bit values.
#' pair <- matrix(data=c(30,70,70,30), nrow=2, ncol=2)
#' util.mtest(pair)
#' 
#' bytevalues <- c(rmultinom(1, 2^20, rep(2^-8, 2^8))) # simulate uniform distributed byte
#' util.mtest(bytevalues) # calculate p-value with as null hypothesis uniform distribution
#' 
#' @param n  A vector of counts of successes, or a matrix representing a two-way
#'           contingency table.
#' @param p0 [optional] Expected probabilities under the null hypothesis. If not given,
#'                      and n is a matrix, these are set to the marginal probabilities
#'                      of the two-way contingency table. Otherwise set to uniform.
#' @param N  [optional] Number of trails performed, set to sum(n) if not given.
#'
#' @references Fuchs, Camil and Kenett, Ron. A test for Detecting Outlying Cells in
#'             the Multinomial Distribution and Two-Way Contingency Tables.
#' @author Mathy Vanhoef
util.mtest <- function(n, p0 = NULL, N = NULL)
{
  if (is.null(N)) N <- sum(n)
  k <- length(n)
  
  # Process (default) arguments
  if (class(n) == "matrix") {
    # For two-way contingency tables we calculate the expected (marginal) probabilities p0.
    # See section 4: "for two-sided test ... using the methods of section 3 with k = IJ."
    pi <- apply(n, c(1), sum) / N
    pj <- apply(n, c(2), sum) / N
    p0 <- sapply(pj, function(pr) pr * pi)
  } else if (is.null(p0)) {
    # Assume uniform distribution
    p0 <- rep(1/k, k)
  }
    
  # formula (2.1): calculate adjusted residuals Z_i
  z <- (n - N*p0) / sqrt(N * p0 * (1 - p0))
  maxz <- max(abs(z))
  
  # Calculate (an upper bound of) the p-value for a two-sided test. Based on formula 3.5,
  # the definitions in Lemma 1, and the remark for two-sided tests in the last paragraph
  # of section 3. We need high precision `mpfr` otherwise `pnorm` returns zero too quickly.
  #pval <- as.numeric(2 * k * (1 - pnorm(mpfr(maxz, 2048))))
  
  # Improved upper bound for the two-sided test by Sidak (1968). See formula (3.8)
  # at the end of section 3. Doesn't offer that much advantage though.
  pval <- as.numeric(1 - (2 * pnorm(mpfr(maxz, 2048))-  1)^k)
  
  if (pval > 1) pval <- 1
  return (pval)
} 
