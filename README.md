# R package: Envelopes and Reduced-rank Regression

A project developing a R-package which enables us to build Reduced-rank envelopes regression models

## Contributors

* __Pedro Orozco del Pino__ (PhD student in Biostatistics at University of Michigan, Ann Arbor)
* __Sunghun Lee__ (Masters student in Applied Statistics at University of Michigan, Ann Arbor)
* __Yixi Li__ (Masters student in Survery Methodology, Statistical Science at University of Michigan, Ann Arbor)

### Reference

This project is based on a research paper: 
R. Dennis Cook, Liliana Forzani, Xin Zhang; Envelopes and reduced-rank regression, Biometrika, Volume 102, Issue 2, 1 June 2015, Pages 439–456, https://doi.org/10.1093/biomet/asv001


### Description

In multivariate regression analysis there is the case in which the coefficients matrix is not full rank. In such cases the reduced-rank methodology can get a better estimate than OLS. Another improvement that can be done in multivariate analysis is when the estimation of the beta matrix can be done with less information without compromissing accuracy. For this, we look a subspace in which the projection the response is independent of the covariates. The envelope regression methodology estimates this subspace and also outperforms OLS. In this project we build a program that combines reduced-rank and envelope regression methodologies. Each of the methodologies don't compite with each other since they solve different challenges. It is because of this that the combine methodology outperforms OLS, reduced-rank and envelopes.

### Data simulation

We compared the performance (in terms of the fit of model) of our new Reduced-rank envelope package and existing R packages (for OLS, Reduced-rank regression and Envelope regression). As you can see from the below plot, the Reduced-rank envelope regression model with our new R package always showed the best performance in all the three cases. In addition, it is also impressive that Reduced-rank envelope was stable between different sample sizes from 160 to 2000. In particular, it is noteworthy that in a small sample size, its performance was superior to other methods.

![ScreenShot](https://github.com/shnlee-ds/Rpackages-Envelopes_and_ReducedRankReg/blob/master/simulation.png)


### Selection of parameter d and u

When you have a data, a function _rrenv.choose_du(X, Y, Beta=NULL)_ in our package recommends the combination of parameter (d, u) for a reduced-rank envelope regression modeling based on the error calculation. When you input the matrices of X, Y and the true Beta matrix(the left case of below figure), the function calculates  ![](https://latex.codecogs.com/gif.latex?%5Cleft%20%5C%7C%20%5Chat%7B%5Cbeta%7D-%5Cbeta%20%5Cright%20%5C%7C_F), while it calculates ![](https://latex.codecogs.com/gif.latex?AIC%20%3D%202N_%7BRE%7D%20-%202L_%7Bu%2Cd%7D%20%3D%20%28p&plus;r-d%29d&plus;r%28r&plus;1%29/2%20-%202%5Chat%7BL%7D_%7Bu%2Cd%7D) when you only input the matrices of X and Y without the true Beta matrix (the right case of below figure).

The function _rrenv.choose_du(X, Y, Beta=NULL)_ returns not only the heatmap but also the best combination of (d,u) based on the calculated errors.

![ScreenShot](https://github.com/shnlee-ds/Rpackages-Envelopes_and_ReducedRankReg/blob/master/choose_du.png) 
