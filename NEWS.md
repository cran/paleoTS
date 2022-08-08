# paleoTS 0.5.3

## Bug fixes
* corrected error that caused the fitting of mode shift models to fail when a start age was specified
for the paleoTS object to be analyzed.
* corrected error in fitMult() that always skipped the calculation of standard errors
* changed fit3models(), fit4models(), and fit9models() so that they did not fail with an obscure error when all sample sizes are equal to 1.


## Other changes
* no other changes
