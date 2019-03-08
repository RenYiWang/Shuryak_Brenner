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

# 1. Need virtual data: start with oversimplified parametrized model; 
#    add noise and some systematic bias; generate virtual data. Then try to 
#    calibrate a model and see what the parameters look like.
# 2. Simulate HG protracted dosing first; eventually use CA pod data to guide simulated data
# 3. Assume K in MSBE is LET independent
# 4. Assume most LET dependent quantities in MSBE have standard LET dependence 

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


# proposed parameter ordering: r > l > t > t0 > q > ep > k > b > a > m > adjustable parameters

# r = radiation dose rate
# l = LET
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

nte_term <- function(t, pe, k) { 
  return(pe * t + 1 / (k * pe ^ (-1))) 
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
  return(b * (1 + a * te_term(r, t) + m * nte_term(t, pe, k)) 
              * (1 - exp(- t / t0))) # Hazard term
}

# We will need simulated data because no constant flux NSRL data is as yet 
# available for the very low fluxes where NTE are expected to dominate.

# The general plan will be the following. Choose a version of the Shuryak-Brenner 
# equation that uses not only dose-rate R and time but also LET as a categorial 
# predictor variable with 5 values. 

# Modified NTE term with categorical LET dependence
mod_nte_term <- function(l, t, pe, k, alpha) { # Includes l and alpha.
  return(alpha * l * nte_term(t, pe, k))
}

# Modified yield equation with categorical LET dependence (MLSBE)
LET_yield <- function(r, l, t, t0, q, k, b, a, m, alpha) { # l is for LET, alpha is a parameter.
  pe <- equilibrium_probability(r, q)
  return(b * (1 + a * te_term(r, t) + m * mod_nte_term(l, t, pe, k, alpha)) 
              * (1 - exp(- t / t0)))
}

# Choose reasonable numerics so that the equation can produce many simulated 
# data points Add some random noise and a systematic bias to the simulated 
# data points. Now devise a similar but not identical version of the 
# Shuryak-Brenner equation and don't fill in all the numerics -- leave in some 
# unknown adjustable parameters. Calibrate these from the simulated noise-added
# bias-added data points. That gives an HZE effect equation with flux and time
# as predictor variables. 

## EGH: assume chosen LET values are 25, 70, 193, 250, 464
## Plotting set up
library(colorspace)
palette <- c("orange", "red", "purple", "blue", "green") # Color palette for LET ions
alt_palette <- rev(heat_hcl(5, h = c(0, -100), c = c(40, 80), l = c(75, 40), power = 3)) 
LET <- c(25, 60, 193, 250, 464) # Ne, Si, Fe 6, Fe 3, Nb

## Figure 0: Noisy plots of MLSBE with arbitrary parameter settings
plot(NULL, xlim = c(0, 100), ylim = c(0, 10000),
     ylab = "Tumor yield count", xlab = "Dose rate (cGy/day)")
for (i in 1:length(LET)) {
  points(1:100, LET_yield(r=0.01 * 1:100, l=LET[i], # dose rate up to 1 cGy/day,
                          # LET taken from values set earlier.
                          t=365, t0=14, q=0.5, # Over a year, saturation point 
                          # is two weeks, 0.5 cGy/day is point in which 
                          # 50% of susceptible cells are activated.
                          k=2, b=2, a=1, m=10, alpha = 0.004) # k, b, a, m, and
                          # alpha chosen arbitrarily.
         + rnorm(100, mean=0, sd=400), # Assumes Gaussian, second argument 
         # can bias the noise.
         col = palette[i])
}


# Repeat four times for the 
# remaining LET values. 


# Choose the LET-independent SLI effect equation effect = alpha * flux * time 
# where alpha is the only adjustable parameter. Generate
# noisy SLI data data= alpha*[flux +C* flux^2]*time by choosing a numerical 
# value for flux and then adjusting C so the quadratic component is, say, 1/10 
# the linear component and then adding noise.

# sli_data <- alpha * (flux + c * flux ^ 2) * time

# Now we have six effect equations, 5 HZE labelled by LET and one SLI. 5 have 
# adjustable parameters and differ only by having different LETs. The other has 
# adjustable parameter alpha. We have six corresponding data sets. Now consider 
# mixtures, e.g. of two HZE ions or of more than two ions. Combine these using 
# incremental effect additivity, which should now be easy since the damage is in 
# fact arriving incrementally.

# Draw some figures, write some paragraphs, publish, devise some jokes to use 
# when we team-tag our talk accepting our joint award of the Gray medal at the 
# 2021 International Radiation Research Society meeting. Celebrate at the 
# I-house with lots of Lox.