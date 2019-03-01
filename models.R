# Copyright:    (C) 2017-2019 Sachs Undergraduate Research Apprentice Program 
#               (URAP) group, University of California at Berkeley.
#               This program and its accompanying materials are distributed 
#               under the terms of the GNU General Public License v3. As 
#               detailed in that license no warranty, explicit or implied, 
#               comes with this suite of R scripts.
# Filename:     virtual_data.R 
# Purpose:      Concerns radiogenic human tumorigenesis. Generates virtual data
#               based on a modified version of the Shuryak-Brenner radiation 
#               model.
# Contact:      Rainer K. Sachs 
# Website:      https://github.com/rainersachs/Shuryak_Brenner
# Attribution:  This R script was developed at UC Berkeley. Written by HG pod
#               and Rainer K. Sachs 2019.

# this and next few lines are in R notation 
# L = LET. 
# lambda = AR interpreted as a parameter in itself replacing the parameter A 
T_0 <- 1 # d (days); L= c(1, 75*1:3) # t_0 is a very short time., during which NTE build up and saturate. 
# next a standard pattern of how parameters depend on LET 
iL1 <- 1 * exp(-2); iL2 <- 1 * exp(-2); iL <- iL1 * L * exp(-iL2 * L) 
lambda0 <- exp(-4); lambda <- lambda0 * iL #Lubin miner data (lubdat) suggests lambda < 365 days 
M <- 0.4 #estimate from lubdat, assumed the same for all LET 
q0 <- 3 * exp(-4); 
q <- q0 / iL # q0 estimated from lubdat; divided by iL because q is in denominator 
rate0 <- 2 * exp(-3); rate <- rate0 / c(1, iL[2:4] / length(L)) # approximates R, the LET dependent dose rate above LEO 
K <- 4 * exp(-3) # from lubdat


# proposed parameter ordering: r > t > t0 > q > ep > k > b > a > m

# r = radiation dose rate
# q = dose rate at which 50% of all NTE-susceptible cells become activated
# ep = equilibrium probability
# t = time
# t0 = t0 is a very short time, during which NTE build up and saturate. 
# k = rate of exponential decay of NTE signals
# b = background tumors
# a = TE tumor induction
# m = maximum nte contribution

#===================== Shuryak-Brenner Equation (SBE) =========================#
equilibrium_probability <- function(r, q) {
  return(1 / (1 + (q / r)))
}

nte_term <- function(t, ep, k) { 
  return(ep * t + 1 / (k * ep ^ (-1))) 
}

te_term <- function(r, t) {
  return(r * t)
}

yield <- function(r, t, q, k, b, a, m) { # Mean yield of tumors per individual.
  pe <- equilibrium_probability(r, q)
  return(b * (1 + a * te_term(r, t) + m * nte_term(t, pe, k)))
}

#================= Modified Shuryak-Brenner Equation (MSBE) ===================#
modified_yield <- function(r, t, t0, q, k, b, a, m) {
  pe <- equilibrium_probability(r, q)
  return(b * (1 + a * te_term(r, t) + m * nte_term(t, pe, k) 
              * (1 - exp(- t / t0)))) # Hazard term
}

