# BBS_Seasonal_Shift

In response to some interesting reviewer comments, we set out to see if population trends for Allen's Hummingbird might be biased due to a climate-change induced, shifting seasonal peak in bird activity. 

We added a seasonal GAM smooth that was allowed to vary by decade to the base GAMYE trend model in bbsBayes.
The seasonal GAM parameters were estimated as random effects using a first-difference time-series structure that allowed the shape of the seasonal pattern to drift across decades. 
If the changing seasonal patterns were causing a reduction in counts, then adjusted trends should be less negative than the original estimates.

In a nutshell, they weren't.

