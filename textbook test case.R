# https://pages.cms.hu-berlin.de/EOL/gcg_quantitative-methods/Lab14_Kriging.html

library("gstat")   # geostatistics
library("mapview") # map plot
library("sf")      # spatial vector data
library("stars")   # spatial-temporal data
library("terra")   # raster data handling 
library("ggplot2") # plotting
mapviewOptions(fgb = FALSE)

meuse <- read_sf('data/meuse.gpkg')
head(meuse)

mapview(meuse['zinc'])


# make eab_dat a spatial object
# https://cengel.github.io/R-spatial/intro.html

eab_spat <- st_as_sf(eab_dat, coords = c("lon", "lat"))

mapview(eab_spat['Y2022'])

# alright were in business
# https://pages.cms.hu-berlin.de/EOL/gcg_quantitative-methods/Lab14_Kriging.html

# make a box around the eab data
bbox <- st_bbox(eab_spat)


# make a box around meuse data

bbox <- st_bbox(meuse)

cell_size <- 40

x <- seq(bbox$xmin, bbox$xmax, by=cell_size)
y <- seq(bbox$ymin, bbox$ymax, by=cell_size)


meuse_grid <- expand.grid(x=x, y=y)
plot(meuse_grid$x, meuse_grid$y, pch=19, cex=0.1)


meuse_grid$tmp <- 1
meuse_grid <- st_as_stars(meuse_grid, crs=st_crs(meuse))
st_crs(meuse_grid) <- st_crs(meuse) # re-assign crs to be safe


zn.idw <- gstat::idw(zinc ~ 1, locations=meuse, newdata=meuse_grid, idp = 2)

zn.idw

plot(rast(zn.idw["var1.pred"]))
plot(meuse["zinc"], col=1, cex=0.5, add=T, type="p")

# lagged scatterplots
gstat::hscat(log(zinc) ~ 1, meuse, (0:9) * 100)


# variogram cloud

meuse.vc <- variogram(log(zinc) ~ 1, meuse, cloud = TRUE)
plot(meuse.vc, ylab=bquote(gamma), xlab=c("h (separation distance in m)"))


# sample variogram
meuse.v <- gstat::variogram(log(zinc) ~ 1, meuse)
plot(meuse.v, ylab=bquote(gamma), xlab=c("h (separation distance in m)"))


# theoretical variogram
myVariogramModel <- vgm(psill=1, "Sph", range=100, nugget=0.5)

plot(myVariogramModel, cutoff=150)



# fit a variogram
meuse.v <- variogram(log(zinc) ~ 1, meuse)

meuse.sph <- vgm(psill=0.6, "Sph", range=800, nugget=0.1)
plot(meuse.v, meuse.sph, cutoff=1000)


######## https://pages.cms.hu-berlin.de/EOL/gcg_quantitative-methods/Lab14_Kriging.html
# Kriging utilizes the theoretical variogram to interpolate values at any location based on distant-variance relationship. 
# We’ll perform Ordinary Kriging at the meuse grid locations. 
# Recall, ordinary kriging has a constant intercept, denoted in the formula with ~ 1:.


# create sample variogram
meuse.v <- gstat::variogram(log(Totals) ~ 1, eab_spat)

# fit variogram model
meuse.vfit <- gstat::fit.variogram(meuse.v, vgm(1, "Sph", 800, 1))

# ordinary kriging
lz.ok <- gstat::krige(log(zinc) ~ 1, meuse, meuse_grid, meuse.vfit)

plot(rast(lz.ok['var1.pred']))
plot(meuse["zinc"], col="blue", cex=0.5, type="p", add=T)
