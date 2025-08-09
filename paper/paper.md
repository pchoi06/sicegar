---
title: 'Modifications to **sicegar**: Analysis of Single-Cell Viral Growth Curves'
tags:
  - R
  - time course data
authors:
  - name: Samuel Butler
    equal-contrib: true 
    affiliation: '1'
  - name: Phineus Choi
    equal-contrib: true 
    affiliation: '1'
  - name: Thomas Matheis
    equal-contrib: true 
    affiliation: '1'
  - name: Mira Terdiman
    affiliation: '1'
  - name: Johanna Hardin
    orcid: 0000-0001-6251-1955
    affiliation: '1'
affiliations:
 - name: Pomona College, United States
   index: 1
date: 13 August 2025
bibliography: paper.bib
---

# Summary

The R package **sicegar** aims to quantify time intensity data by using sigmoidal and double sigmoidal curves.
It fits straight lines, sigmoidal, and double sigmoidal curves on to time vs intensity data.
Each of the fits are used to make a decision on which model best describes the data. 
The method was originally developed in the context of single-cell viral growth analysis (for details, see @caglar2018), and the package name stands for "SIngle CEll Growth Analysis in R". 
Beyond **sicegar**'s ability to categorize fits, it also provides parameter estimations for each curve which can also provide important information to researchers. 

In particular, the sigmoidal function is given as follows (with parameters estimated using the Levenberg-Marquardt algorithm [@levenberg1944;@marquardt1963]).
$h_0$ represents the lower asymptote (as $x$ approaches negative infinity).
$t_1$ is the onset time, the midpoint between $h_0$ and $h_1$.
$a$ determines the magnitude of the slope of the sigmoidal curve.
$h_1$ is the upper asymptote (as $x$ approaches positive infinity).

\begin{equation}
I(x) = h_0 + \frac{h_1-h_0}{1 + e^{-a(x - t_1)}}
\end{equation}

A similar, but slightly more complicated, formula for the double sigmoidal function is also used for parameter estimation in **sicegar**.
In the original implementation of **sicegar**, the parameter $h_0$ is set to zero for both the sigmoidal and double-sigmoidal models.

@Adams uses **sicegar** to investigate the onset time of RNA expression in genes in *E. coli* undergoing stress. 
They observed limitations in **sicegar**'s fits, which motivated the improvements proposed in this paper. 
Our primary update to the **sicegar** package is the estimation of an additional parameter, $h_0$, which is the lower asymptote of both the sigmoidal and double-sigmoidal curves, which had previously been set equal to 0.
Based on simulated data, the free estimation of $h_0$ provides both a better fit (lower SSE) and more accurate parameter estimations than when $h_0$ is forced to be equal to zero. 
Other smaller adjustments to the package include improvements to how one of the parameters is un-normalized and adjusting the threshold for an error catch. 
For backward compatibility, the updated package is designed so that the $h_0 = 0$ is the default value.
As seen in **Figure 1**, the argument `use_h0` in the overarching function `fitAndCategorize` allows the user to decide whether to allow the package to estimate $h_0$.


# Statement of Need

@caglar2018 discuss **sicegar**'s ability to correctly identify sigmoidal and double-sigmoidal curves on simulated data.
They report, "Overall, we can conclude that our algorithm results in reliable fits, that it fails gradually with increasing noise levels, and that it is conservative in assessing whether it has correctly identified a sigmoidal or double-sigmoidal curve or not."
The focus of their results is on whether **sicegar** is able to correctly identify sigmoidal and double-sigmoidal curves on simulated data, not on whether it is able to accurately report parameter estimations.
Though they acknowledge **sicegar**'s ability to estimate parameters associated with sigmoidal and double-sigmoidal curves, they do not report on the accuracy of the parameter estimates.


In dozens of research projects that use **sicegar** for modeling time-intensity data, researchers are interested in extracting specific parameter estimates, like midpoints and slopes, as the parameter values represent biologically meaningful information.
@Adams extracted midpoint values ($t_1$) to investigate onset time of RNA expression in genes in *E. coli* undergoing stress.
@wittemeier used **sicegar** to estimate molar carbon assimilation. 
They used **sicegar**'s estimation of maximum slope in sigmoidal curves to understand maximum assimilation rate, and midpoints ($t_1$) in sigmoidal curves to extract the point at which maximum assimilation is reached. 
@rajarathinam also used **sicegar** to analyze carbon assimilation and to extract estimations of maximum slope.
Our addition of $h_0$ to the set of estimated parameters greatly improves the package's ability to provide accurate parameter estimations.
Through simulations, with varying levels of noise, generating parameters, and both sigmoidal and double sigmoidal curves, we were able to accurately estimate all parameter estimates, even when the lower asymptote is not zero.
Though the categorization of the model as sigmoidal or double-sigmoidal is hugely important, it is not the only important aspect of the **sicegar** modeling.
Our updated implementation, which includes the estimation of the lower asymptote, is prevailingly important, and thus our adjustments fit the needs of current research. 

# Features 

### Core Functions

`fitAndCategorize`, the overarching function in **sicegar**, takes time-intensity data as an argument and runs the data through a series of nested functions.
The structure of this process is outlined in **Figure 1**.
First, the data are normalized in the function `normalizeData` (1).
Then they are passed through `multipleFitFunction` (2) which uses the Levenberg-Marquardt algorithm [@levenberg1944;@marquardt1963] to fit both sigmoidal and double-sigmoidal curves to the data, (3), (4).
Additional parameters are then added to the output vector in `parameterCalculations` to prepare the model to be plotted.
The user decides whether to allow the function to estimate $h_0$ using the argument `use_h0 = TRUE` in `fitAndCategorize`. 
If `use_h0` is set to `FALSE` (the default), the algorithm will run **sicegar** as it was originally written, with $h_0$ fixed at zero.
If they choose to allow the function to estimate $h_0$, it will follow the same function flow, except that each function will account for the estimation of $h_0$.

### Novel Contributions

Each function in **sicegar** was rewritten to include the parameter $h_0$ in addition to a some small technical changes to the functions. 
Our new version is outlined in the right-hand branch of **Figure 1**.

**Figure 1** Structure of the `fitAndCategorize` function.

![Structure of the `fitAndCategorize` function](JOSS_Graphic.png)

# Example

I'm open to feedback here for sure.
Totally fine if you don't like it and take it out!
And the code I've included isn't sufficient, because maybe we need pictures like in the vignette.
But... I'm thinking we might be able to run a mini-version of the example at the bottom of the new vignette?
https://hardin47.github.io/sicegar/articles/h0_functions.html

(We will have to run the code in R / Rmd, save the image(s) if needed, and then copy over the code into this md file)

I'm including some code here. 
What do you think? 
Also, I think Tommy is going to change the example, so probably all this will need to change.

```r
time <- seq(1, 24, 0.5)
noise_parameter <- 0.2
intensity_noise <- runif(n = length(time), min = 0, max = 1) * noise_parameter
intensity <- doubleSigmoidalFitFormula_h0(time,
                                       finalAsymptoteIntensityRatio = .3,
                                       maximum = 10,
                                       slope1Param = 1,
                                       midPoint1Param = 7,
                                       slope2Param = 1,
                                       midPointDistanceParam = 8,
                                       h0 = 2)
intensity <- intensity + intensity_noise
dataInput <- data.frame(time, intensity)
```

Recall that the original model parameters (which generated the data) are given as `finalAsymptoteIntensityRatio = 0.3`, `maximum = 10`, `slope1Param = 1`, `midPoint1Param = 7`, `slope2Param = 1`, `midPointDistanceParam = 8`, `h0 = 2`.

```r
fitObj_zero <- fitAndCategorize(dataInput,
                           threshold_minimum_for_intensity_maximum = 0.3,
                           threshold_intensity_range = 0.1,
                           threshold_t0_max_int = 0.05,
                           use_h0 = FALSE)   # Default
```

```
$finalAsymptoteIntensityRatio_Estimate
[1] 0.1264636

$maximum_Estimate
[1] 10.13

$slope1Param_Estimate
[1] 1.007482

$midPoint1Param_Estimate
[1] 7.006687

$slope2Param_Estimate
[1] 1.021521

$midPointDistanceParam_Estimate
[1] 7.976228
```

```r
fitObj_free <- fitAndCategorize(dataInput,
                           threshold_minimum_for_intensity_maximum = 0.3,
                           threshold_intensity_range = 0.1,
                           threshold_t0_max_int = 0.05,
                           use_h0 = TRUE)
```

```
$finalAsymptoteIntensityRatio_Estimate
[1] 0.3080899

$maximum_Estimate
[1] 10.13

$slope1Param_Estimate
[1] 1.007482

$midPoint1Param_Estimate
[1] 7.006687

$slope2Param_Estimate
[1] 1.021521

$midPointDistanceParam_Estimate
[1] 7.976228

$h0_Estimate
[1] 2.106237
```



# Availability

The **sicegar** package is available on CRAN (https://CRAN.R-project.org/package=sicegar) and GitHub (https://github.com/hardin47/sicegar).
Documentation, including vignettes and examples, is provided to facilitate adoption.

# Acknowledgements

The authors gratefully acknowledge Dan Stoebel for bringing the application to our attention and Federica Domecq Lacroze for sharing her explorations of the **sicegar** package.





# References









