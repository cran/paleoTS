---
title: "NEWS"
author: "Gene Hunt"
date: "3/12/2019"
output: html_document
---



## Changes for v. 0.5-2

## Bug fixes
* sim.OU() now returns a vector for $mm, rather than a matrix
* fitMult() now correctly passes the hess argument to the fitting functions
* corrected error in sub.paleoTS() which caused that function to fail if one omitted the first point of a series
* added a check to test.var.het() so that it gracefully handles when all vv = 0


## Other changes
* changed plotting options so that the bg for points is now "white" and model fits are shown in grey
* added argument to compareModels(sort = FALSE); if TRUE the results table is sorted from best to worst model. If FALSE, no sorting is done.

