---
title: 'Modifications to Sicegar: Analysis of Single-Cell Viral Growth Curves'
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

Sicegar aims to quantify time intensity data by using sigmoidal and double sigmoidal curves.
It fits straight lines, sigmoidal, and double sigmoidal curves on to time vs intensity data.
Then all the fits are used to make a decision on which model best describes the data. 
his method was first developed in the context of single-cell viral growth analysis (for details, see @caglar2018), and the package name stands for "SIngle CEll Growth Analysis in R". 
Beyond Sicegar's ability to categorize fits, it also provides parameter estimations for each curve which can also provide important information to researchers. 
@Adams used Sicegar to investigate the onset time of RNA expression in genes in $\textit{E. Coli}$ undergoing stress. 
They observed limitations in Sicegar's fits, which motivated the improvements proposed in this paper. 
The primary change is the estimation of an additional parameter, $h_0$, which is the lower asymptote of both the sigmoidal and double sigmoidal curves. 
It had previously been set equal to 0. Based on thousands of iterations on simulated data, we've concluded that the free estimation of $h_0$ both provides better fits (lower SSE), and more accurate parameter estimations than when $h_0$ is forced to be equal to 0. 
We've also made some smaller adjustments to the package with respect to how one of the parameters is un-normalized, and adjusting the threshold for an error catch, however these are more technical and less important than the primary change which is the free estimation of $h_0$. 
The package is designed so that previous work done using Sicegar can be reproduced. 
As seen in Figure blah blah blah, the argument `use_h0` in the overarching function `fitAndCategorize` allows the user to decide whether to allow the package to estimate $h_0$.


# Statement of Need

the peer j article has been cited 46 times - Dimensions Badge, which is a free tool from dimensions ai that provides visualizations of citation data for scholarly publications states with regards to the peerj article: Compared to other publications in the same field, this publication is extremely highly cited and has received approximately 5.29 times more citations than average.

@wittemeier used Sicegar to estimate molar carbon assimilation. They used Sicegar's estimation of maximum slope in sigmoidal curves to understand maximum assimilation rate, and $t_1$ in sigmoidal curves to extract the point at which maximum assimilation is reached. (this is important because they extracted parameter estimations and thats what we improved)

@rajarathinam also used Sicegar to analyze carbon assimilation, and to extract estimations of maximum slope.

# Features 
(or maybe "Useage" instead of "Features"????)

### Core Functions

### Novel Contributions

![Structure of the original and revised `fitAndCategorize` function](JOSS Graphic.png)

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

The voice package is available on CRAN (https://CRAN.R-project.org/package=sicegar) and GitHub (https://github.com/hardin47/sicegar).
Documentation, including vignettes and examples, is provided to facilitate adoption.

# Acknowledgements

The authors gratefully acknowledge Dan Stoebel for bringing the application to our attention and Federica Domecq Lacroze for sharing her explorations of the **sicegar** package.





# References









