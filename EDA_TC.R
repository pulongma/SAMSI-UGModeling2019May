
rm(list=ls())

setwd("~/Misc/SAMSI/MUMS/UGModeling")


### load library to read .nc files
library(ncdf4)

### load all storms available in the IBTrACS record
#storms <- nc_open("./data/IBTrACS.ALL.v04r00.nc")
storms = nc_open("./data/IBTrACS.since1980.v04r00.nc")
#print(storms)
#nc_close(storms)

### Storm name
name = ncvar_get(storms, "name")

nstorms = length(name)

### Year based on season
season = ncvar_get(storms, "season")
count = as.numeric(table(season))
year = as.numeric(names(table(season)))
parset = par(las=1)
plot(year, count, type="h", col="blue")
par(parset)

### Type of storm
nature = ncvar_get(storms, "nature")
barplot(round(table(nature)[2:7]/360))

### Storm at each basin
basin = ncvar_get(storms, "basin")
barplot(round(table(basin)[2:8]/360))


### Storm center 
Lat = ncvar_get(storms, "lat")
Lon = ncvar_get(storms, "lon")

### Distance to land: km
dist2land = ncvar_get(storms, "dist2land")

### Maximum sustained wind speed
mws = ncvar_get(storms, "wmo_wind") * 0.514444  # kt to m/s

summary(c(mws))

hist(mws, 100, col="lightblue", border="gray",
     main="Maximum Sustained Wind Speed", xlab="m/s")

### Saffir-Simpson hurricane wind scale (in m/s)
# https://en.wikipedia.org/wiki/Saffir–Simpson_scale

tstep = dim(mws)[1]

SSHS = matrix(NA, tstep, nstorms)
SSHS[mws<33] = 0
SSHS[mws>=33 & mws<43] = 1
SSHS[mws>=43 & mws<50] = 2
SSHS[mws>=50 & mws<58] = 3
SSHS[mws>=58 & mws<70] = 4
SSHS[mws>=70] = 5


### Global map for MWS
library(mapproj)
library(fields)
library(scales)
library("viridis")

mws.df = data.frame(long=c(Lon), lat=c(Lat), value=c(mws))
mws.df$long[mws.df$long>180] = NA
mws.df = mws.df[complete.cases(mws.df), ]


proj.type = "mollweide"
proj.orient = c(90, 0, 0)

png("./figures/global_map_mws.png", width=800, height=400, 
    units="px", pointsize=14)
parset = par(mar=c(0,0,0,5)+1)
map("world", mar=rep(0,4), 
    proj=proj.type, orient=proj.orient,
    col=alpha("gray50", 0.8), fill=T, lwd=0.05)
map.grid(col="gray", labels=F, lty=2)
image.plot(legend.only=T, zlim=c(0,95), col=viridis(6),
           legend.line=0.5)
xyproj = mapproject(x=mws.df$long, y=mws.df$lat, proj=proj.type, orient=proj.orient)
val = color.scale(c(mws.df$value), col=alpha(viridis(6), 1))
points(xyproj$x, xyproj$y, col=val, pch=19, cex=0.25)
title("Maximum sustained wind speed")
par(parset)
dev.off()





### Global map for SSHS
library(mapproj)
library(fields)
library(scales)
library("viridis")

sshs.df = data.frame(long=c(Lon), lat=c(Lat), value=c(SSHS))
sshs.df$long[sshs.df$long>180] = NA
sshs.df = sshs.df[complete.cases(sshs.df), ]

proj.type = "mollweide"
proj.orient = c(90, 0, 0)

png("./figures/global_map_sshs.png", width=800, height=400, 
    units="px", pointsize=14)
parset = par(mar=c(0,0,0,5)+1)
map("world", mar=rep(0,4), 
    proj=proj.type, orient=proj.orient,
    col=alpha("gray50", 0.95), fill=T, lwd=0.05)
map.grid(col="gray50", labels=F, lty=2)
image.plot(legend.only=T, zlim=c(0,5), col=viridis(6))
xyproj = mapproject(x=sshs.df$long, y=sshs.df$lat, proj=proj.type, orient=proj.orient)
val = color.scale(sshs.df$value+1, col=alpha(viridis(6), 1))
points(xyproj$x, xyproj$y, col=val, pch=15, cex=0.25)
title("Saffir-Simpson hurricane scale")
par(parset)
dev.off()





### North Atlantic

### mws 
window = (mws.df$long>-100 & mws.df$long<=-60 & 
            mws.df$lat>=15 & mws.df$lat<=50)
mws.df.NA = mws.df[window, ]

png("./figures/global_map_mws_NA.png", width=800, height=600, 
    units="px", pointsize=16)
parset = par(mar=c(0,0,0,5)+2)
map("world", mar=rep(0,4)+1, xlim=c(-95,-60),
    ylim=c(15,50), col="white", fill=T, lwd=0.05)

image.plot(legend.only=T, zlim=c(0,95), col=viridis(6))
val = color.scale(mws.df.NA$value, col=alpha(viridis(6), 1))
points(mws.df.NA$long, mws.df.NA$lat, col=val, pch=16, cex=0.6)

title("Maximum sustained wind speed")
map("world", mar=rep(0,4)+1, xlim=c(-95,-60),
    ylim=c(15,50), lwd=1.25, add=T)
map("state", mar=rep(0,4)+1, xlim=c(-95,-60),
    ylim=c(15,50), lwd=1.25, add=T)
par(parset)
dev.off()

### SSHS 
window = (sshs.df$long>-100 & sshs.df$long<=-60 & 
            sshs.df$lat>=15 & sshs.df$lat<=50)
sshs.df.NA = sshs.df[window, ]

png("./figures/global_map_sshs_NA.png", width=800, height=600, 
    units="px", pointsize=16)
parset = par(mar=c(0,0,0,5)+2)
map("world", mar=rep(0,4)+1, xlim=c(-95,-60),
    ylim=c(15,50), col="white", fill=T, lwd=0.05)

image.plot(legend.only=T, zlim=c(0,5), col=viridis(6))
val = color.scale(sshs.df.NA$value+1, col=alpha(viridis(6), 1))
points(sshs.df.NA$long, sshs.df.NA$lat, col=val, pch=16, cex=0.6)

title("Saffir-Simpson hurricane scale")
map("world", mar=rep(0,4)+1, xlim=c(-95,-60),
    ylim=c(15,50), lwd=1.25, add=T)
map("state", mar=rep(0,4)+1, xlim=c(-95,-60),
    ylim=c(15,50), lwd=1.25, add=T)
par(parset)
dev.off()

