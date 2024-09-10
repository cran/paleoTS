# paleoTS 0.6-2

## Bug fix
* corrected a bug, reported by N. Hohmann, in test.var.het() so that it gracefully handles a situation 
where some sample variances are zero.

## Minor changes
* Added a warning message to fit3models() and fit4models() that cautions against interpreting model support
for very short sequences.


# paleoTS 0.6.1

## Major changes
* added capacity to fit many models using state-space models (SSMs) and the Kalman filter
* model-fitting functions now return the convergence code from optim()
* added additional information to paleoTSfit class, which also has a new print method to show an informative summary of these objects
* improved default setting for all optimization functions that improves convergence behavior
* function akaike.wts() is now documented and exported


## Minor changes and bug fixes
* the plot method for paleoTS objects now allows the user more control, including specifying axis labels
* bug fix for fitMult() that previously wouldn't use the hess argument even when specified
* standardized how bounds are handled in L-BFGS-B optimization for non-negative parameters
* bug fix for the function sub.paleoTS(), which returned an incorrect start.age in some situations
