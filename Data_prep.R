### modeling the seasonal pattern in BBS observations through time


# packages ----------------------------------------------------------------

library(bbsBayes)
library(tidyverse)
library(mgcv)
library(cmdstanr)




# prep the BBS data -------------------------------------------------------

st_dat = stratify(by = "bbs_usgs")

species = "Allen's Hummingbird"


spsData = prepare_jags_data(st_dat,species_to_run = species,
                            model = "slope")


#require(lubridate)
spsData$doy = lubridate::yday(as.Date(paste(spsData$r_year,
                            spsData$month,
                            spsData$day,
                            sep = "-")))

spsData$decadeF = cut(spsData$r_year,breaks = c(1966,
                                                1979.5,
                                                1989.5,
                                               1999.5,
                                               2009.5,
                                               2019.5),
                     labels = c("70s",
                                "80s",
                                "90s",
                                "00s",
                                "10s"),
                     ordered_result = TRUE)

spsData$decade = as.integer(spsData$decadeF)

tmp = data.frame(year = spsData$year,
                 route = factor(spsData$route),
                 obser = factor(spsData$obser),
                 count = spsData$count,
                 lcount = log(spsData$count+1,base = 10),
                 doy = spsData$doy,
                 decadeF = spsData$decadeF,
                 yr_d = as.integer(str_sub(spsData$r_year,4,4)),
                 week = cut(spsData$doy,breaks = seq(134,195,by = 7)))



bp = ggplot(data = tmp,aes(y = lcount,x = doy,colour = yr_d))+
  #geom_boxplot(varwidth = TRUE)+
  #geom_violin()+
  geom_point(alpha = 0.3,position = position_jitter(width = 0.25))+
  #scale_y_log10()+
  scale_colour_viridis_c(begin = 0.1,end = 0.7)+
  geom_smooth()+
  facet_wrap(~decadeF,nrow = 2)
print(bp)





