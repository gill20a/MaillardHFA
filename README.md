# MaillardHFA
## Code associated with weibull model fits and figures for Maillard HFA decomposition study
This repository contains all code to replicate model fits, analyses, and figures associated with *Maillard et al. XXXX. Title*.

This project contains two complementary studies exploring home/away effects on the *Quercus petraea* leaf litter decomposition process. In Experiment 1, leaf litter was reciprocally transplanted between two sites (factorial cross) and harvested at approximately 3-month intervals over a two year period to evaluate the rate of substrate mass loss. In Experiment 2, leaf litter from both sites was decomposed at a single site and harvested at approximately 6-month intervals, and additional micorbial community and decomposition processes were evaluated at each harest point.

Litter bags were deployed in a transect and removed randomly at each harvest timepoint (**Experiment one**: 3 litter bags per harvest; **Experiment two**: 13 litter bags per harvest). We used weibull decomposition models (Cornwell *et al.* 2014) to describe  substrate decomposition dynamics throughout each litter decay period. The weibull model characterizes the litter decay process as a continuous Weibull distribution of residence times described by scale (β),  and shape (α) parameters: 

*X* = e<sup>-(t/β)<sup>α</sup>	</sup>		

We then used model parameters to calculate the time to 10, 25 and 50% mass loss (*t<sub>1/10</sub>* , *t<sub>1/4</sub>* ,*t<sub>1/2</sub>*) using the function:

*t<sub>(1-p)</sub>*=β(ln 1/p)<sup>1/α</sup> , where *p* is the proportion of litter mass remaining. We also calculated the litter mean residence time (*MRT*):

*MRT*= βΓ(1+1/α), where Γ is the gamma function. The ratio of the litter *MRT* to *t<sub>1/2</sub>* provides information about the shape of the substrate loss curve. When *MRT/t<sub>1/2</sub>* is high, labile substrate pools decompose quickly, while substrate loss at later stages of decay is slow. *MRT/t<sub>1/2</sub>* is typically correlated with the α shape parameter, where lower α values indicate larger differences between the early and late stage litter decomposition rates.

In addition to the weibull model approach, we also fit standard first-order exponential decomposition models to the litter bag harvest series:

*X*=e<sup>-k<sub>s</sub>t</sup>, where *k<sub>s</sub>* defines a constant rate of substrate mass loss per unit time throughout the decomposition period. We  compared model performance using individual model AIC and RMSE. 

## **Model Fitting**
*Experiment 1:* Only three litter bags were harvested at each time point during experiment one, limiting our ability to estimate variability in model fit using the bootstrapping appraoch described below. Therefore, we fit both single and weibull models to all harvest points for each litter decomposition treatment (Home/Away). 
 - **Data**: MaillardHFA_Exp1_masslossdata.csv
 - **MetaData**: MaillardHFA_Exp1_masslossdata_metadata.csv
 - **Code**: MaillardHFA_Exp1_modelfitting.R
 
 *Experiment 2:* 

**References**
Cornwell WK, Weedon JT (2014) Decomposition trajectories of diverse litter types: a model selection analysis. Methods in Ecology and Evolution 5(2):173–182.


