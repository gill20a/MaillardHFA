# MaillardHFA
## Code associated with weibull model fits and figures for Maillard HFA decomposition study
This repository contains all code to replicate model fits, analyses, and figures associated with *Maillard et al. XXXX. Title*.

This project contains two complementary studies exploring home/away effects on the *Quercus petraea* leaf litter decomposition process. In Experiment 1, leaf litter was reciprocally transplanted between two sites (factorial cross) and harvested at approximately 3-month intervals over a two year period to evaluate the rate of substrate mass loss. In Experiment 2, leaf litter from both sites was decomposed at a single site and harvested at approximately 6-month intervals, and additional micorbial community and decomposition processes were evaluated at each harest point.

Litter bags were deployed in a transect and removed randomly at each harvest point (experiment one: 3 litter bags per harvest; experiment two: 13 litter bags per harvest). We fit weibull decomposition models to the fraction litter mass remaining at each litter harvest point to describe substrate decomposition dynamics throughout each litter decay period. The weibull model characterizes the litter decay process as a continuous Weibull distribution of residence times described by scale (β),  and shape (α) parameters: 

t^2^~2~ 

*X*=e^(〖-(t/β)〗^α )			

We then used model parameters to calculate the time to 10, 25 and 50% mass loss (t^~1/10~ , t^1/4, t^½.  
t^(1-p)=β(lna(1/p) )^(1/α)				

, where p is the proportion of litter mass remaining. We also calculated the litter mean residence time (MRT; (28)).
MRT= βΓ(1+1/α)										Eqn. 6.
, where Γ is the gamma function. The ratio of the litter MRT to t1/2 provides information about the shape of the substrate loss curve. When MRT/( t_(1/2) ) is high, labile substrate pools decompose quickly, while substrate loss at later stages of decay is slow (28). 


