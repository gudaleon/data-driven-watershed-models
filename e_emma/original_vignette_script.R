###########################################################################################
#
# Demonstration vignette script for the ensemble end-member mixing analysis routine eemma.R




# version 1  build 2023.04.14
# Author: James Kirchner, ETH Zurich
#
# Copyright (C) 2023 ETH Zurich and James Kirchner
# Public use of this script is permitted under GNU General Public License 3 (GPL3); for details see <https://www.gnu.org/licenses/>
# 
#
# READ THIS CAREFULLY:
# ETH Zurich and James Kirchner make ABSOLUTELY NO WARRANTIES OF ANY KIND, including NO WARRANTIES, expressed or implied, that this software is
#    free of errors or is suitable for any particular purpose.  Users are solely responsible for determining the suitability and
#    reliability of this software for their own purposes.





setwd("~/EEMMA paper")   # set this to wherever eemma.R and your data are

rm(list=ls())            # clean up the environment

source("eemma.R")        # source the eemma script



# read in the isotope data
dat <- read.table(file="vignette_data.txt", header=TRUE, sep="\t")


# if we need to see what the data look like...:
str(dat)


# for convenience, unbind the data frame
list2env(dat, envir=.GlobalEnv)


# now let's try a first analysis: xylem1 as a mixture of precipitation, soil water, and groundwater
mm <- eemma(y=xylem1, X=cbind(precip, soil, gw))
mm$coeffs

# or we can do this in one step
eemma(y=xylem1, X=cbind(precip, soil, gw))$coeffs


# what if we only had measurements of precipitation and soil water, but not groundwater?
eemma(y=xylem1, X=cbind(precip, soil), complete=FALSE)$coeffs


# OK, so if we have a missing end-member, can we estimate its tracer time series?
ts <- eemma(y=xylem1, X=cbind(precip, soil), complete=FALSE)$timeseries
plot(timestep, ts$unk, col='blue', pch=19, ylim=c(min(ts$unk-ts$unk_std_err), max(ts$unk+ts$unk_std_err)))
arrows(timestep, ts$unk-ts$unk_std_err, timestep, ts$unk+ts$unk_std_err, length=0.05, angle=90, code=3, col="blue")
lines(timestep, gw, col='gray', lwd=2)


# and what happens if we specify complete=FALSE, but the model is actually complete?
eemma(y=xylem1, X=cbind(precip, soil, gw), complete=FALSE)$coeffs


# and what happens if we just fit to precipitation and soil water, without even testing for missing end-members?
eemma(y=xylem1, X=cbind(precip, soil))$coeffs


# now let's look for memory effects, knowing we shouldn't find them...
eemma(y=xylem1, X=cbind(precip, soil, gw), memory=TRUE)$coeffs


# and then look for memory effects in the other xylem time series, where we certainly should find them:
eemma(y=xylem2, X=cbind(precip, soil, gw), memory=TRUE)$coeffs


# note that if we ignore these memory effects, we get distorted mixing coefficients:
eemma(y=xylem2, X=cbind(precip, soil, gw), memory=FALSE)$coeffs


# now we can look for sources of streamwater
eemma(y=stream, X=cbind(precip, soil, gw))$coeffs

# and test for missing end-members (we won't find any)
eemma(y=stream, X=cbind(precip, soil, gw), complete=FALSE)$coeffs


# and test for memory effects (again, we won't find any)
eemma(y=stream, X=cbind(precip, soil, gw), memory=TRUE)$coeffs


# what if we only had measurements of precipitation?  Would we be able to tell that something was missing?
eemma(y=stream, X=cbind(precip), complete=FALSE)$coeffs


# OK, so if we have a missing soil + groundwater end-member, can we estimate its tracer time series?
ts <- eemma(y=stream, X=cbind(precip), complete=FALSE)$timeseries
plot(timestep, ts$unk, col='orange', pch=19, ylim=c(min(ts$unk-ts$unk_std_err), max(ts$unk+ts$unk_std_err)))
arrows(timestep, ts$unk-ts$unk_std_err, timestep, ts$unk+ts$unk_std_err, length=0.05, angle=90, code=3, col="orange")
lines(timestep, 0.25*soil+0.75*gw, col='gray', lwd=2)


# OK, so if we have a missing soil + groundwater end-member, can we estimate its tracer time series? (second plot)
ts <- eemma(y=stream, X=cbind(precip), complete=FALSE)$timeseries
plot(timestep, precip, col='lightblue', pch=19, ylim=c(min(precip), max(precip)))
points(timestep, stream, col='blue', pch=19)
lines(timestep, stream, col='blue', pch=19, lwd=2)
points(timestep, ts$unk, col='orange', pch=19)
arrows(timestep, ts$unk-ts$unk_std_err, timestep, ts$unk+ts$unk_std_err, length=0.05, angle=90, code=3, col="orange")
lines(timestep, 0.25*soil+0.75*gw, col='gray', lwd=2)




