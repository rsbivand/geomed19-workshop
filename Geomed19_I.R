

## ---- echo = TRUE--------------------------------------------------------
library(sf)

## ---- echo = TRUE, cache=TRUE--------------------------------------------
olinda <- st_read("data/olinda.gpkg")
st_crs(olinda) <- 22525

## ---- echo = TRUE--------------------------------------------------------
library(mapview)
mapview(olinda)


## ---- echo = TRUE--------------------------------------------------------
all(st_is(olinda, "XY"))
str(st_coordinates(st_geometry(olinda)[[1]]))


## ---- echo = TRUE, mysize=TRUE, size='\\tiny'----------------------------
library(sp)
str(slot(as(st_geometry(olinda), "Spatial"), "polygons")[[1]])


## ---- echo = TRUE, mysize=TRUE, size='\\tiny'----------------------------
strwrap(st_as_text(st_geometry(olinda)[[1]]))


## ---- echo=TRUE----------------------------------------------------------
sf_extSoftVersion()


## ---- echo=TRUE----------------------------------------------------------
p1 = st_point(c(3,5))
class(p1)
p2 = st_point(c(4,6))
p3 = st_point(c(4,4))
pts = st_sfc(p1, p2, p3)
class(pts)
sf = st_sf(a = c(3,2.5,4), b = c(1,2,4), geom = pts)
class(sf)
sf


## ---- echo=TRUE----------------------------------------------------------
olinda_sirgas2000 <- st_transform(olinda, 31985)
st_crs(olinda_sirgas2000)


## ---- echo=TRUE----------------------------------------------------------
library(sf)
nc <- st_read(system.file("shapes/sids.shp", package="spData")[1], quiet=TRUE)
st_crs(nc) <- "+proj=longlat +datum=NAD27"
row.names(nc) <- as.character(nc$FIPSNO)
head(nc)


## ---- echo=TRUE, warning=FALSE-------------------------------------------
library(spdep)
gal_file <- system.file("weights/ncCR85.gal", package="spData")[1]
ncCR85 <- read.gal(gal_file, region.id=nc$FIPSNO)
ncCR85


## ---- echo=TRUE, warning=TRUE, out.width='90%', fig.align='center', width=7, height=4----
plot(st_geometry(nc), border="grey")
plot(ncCR85, st_centroid(st_geometry(nc), of_largest_polygon), add=TRUE, col="blue")


## ---- echo=TRUE----------------------------------------------------------
set.seed(1)
nc$rand <- rnorm(nrow(nc))
lw <- nb2listw(ncCR85, style="B")
moran.test(nc$rand, listw=lw, alternative="two.sided")


## ---- echo=TRUE----------------------------------------------------------
nc$LM <- as.numeric(interaction(nc$L_id, nc$M_id))
alpha <- 1
beta <- 0.5
sigma <- 2
nc$trend <- alpha + beta*nc$LM + sigma*nc$rand
moran.test(nc$trend, listw=lw, alternative="two.sided")


## ---- echo=TRUE----------------------------------------------------------
lm.morantest(lm(trend ~ LM, nc), listw=lw, alternative="two.sided")


## ---- echo=TRUE----------------------------------------------------------
aggLM <- aggregate(nc[,"LM"], list(nc$LM), head, n=1)
(aggnb <- poly2nb(aggLM))


## ---- echo=TRUE----------------------------------------------------------
plot(st_geometry(aggLM))


## ---- echo=TRUE----------------------------------------------------------
set.seed(1)
LMrand <- rnorm(nrow(aggLM))


## ---- echo=TRUE----------------------------------------------------------
moran.test(LMrand, nb2listw(aggnb, style="B"))


## ---- echo=TRUE----------------------------------------------------------
nc$LMrand <- LMrand[match(nc$LM, aggLM$LM)]
plot(nc[,"LMrand"])


## ---- echo=TRUE----------------------------------------------------------
moran.test(nc$LMrand, listw=lw, alternative="two.sided")


## ---- echo = TRUE--------------------------------------------------------
library(elevatr)


## ---- echo = TRUE, cache=TRUE--------------------------------------------
elevation <- get_elev_raster(as(olinda_sirgas2000, "Spatial"), z = 14)
elevation[elevation$layer < 1] <- NA


## ---- echo = TRUE--------------------------------------------------------
mapview(elevation, col=terrain.colors)


## ---- echo = TRUE--------------------------------------------------------
elevation


## ---- echo = TRUE--------------------------------------------------------
sp::gridparameters(as(elevation, "SpatialGrid"))


## ---- echo = TRUE--------------------------------------------------------
library(stars)
e1 <- st_as_stars(elevation)
e1


## ---- echo = TRUE--------------------------------------------------------
fn <- system.file("tif/L7_ETMs.tif", package = "stars")
system.time(L7 <- read_stars(fn))
L7


## ---- echo = TRUE--------------------------------------------------------
plot(L7)


## ---- echo=TRUE----------------------------------------------------------
ndvi <- function(x) (x[4] - x[3])/(x[4] + x[3])
(s2.ndvi <- st_apply(L7, c("x", "y"), ndvi))
system.time(plot(s2.ndvi)) 


## ---- echo = TRUE--------------------------------------------------------
system.time(L7 <- read_stars(fn, proxy=TRUE))
L7


## ---- echo=TRUE----------------------------------------------------------
(s2.ndvi = st_apply(L7, c("x", "y"), ndvi))
system.time(plot(s2.ndvi)) 


## ---- echo=TRUE----------------------------------------------------------
library(stars)


## ---- echo=TRUE----------------------------------------------------------
tif <- system.file("tif/L7_ETMs.tif", package = "stars")
x <- read_stars(tif)

## ---- echo=TRUE----------------------------------------------------------
str(x)


## ---- echo=TRUE----------------------------------------------------------
plot(x, rgb = c(3, 2, 1))


## ---- echo=TRUE----------------------------------------------------------
plot(x, rgb = c(4, 3, 2))


## ---- echo=TRUE----------------------------------------------------------
(x6 <- split(x, "band"))


## ---- echo=TRUE----------------------------------------------------------
plot(x6)

## ---- echo=TRUE----------------------------------------------------------
str(x6)


## ---- echo=TRUE----------------------------------------------------------
x6$mean <- (x6[[1]] + x6[[2]] + x6[[3]] + x6[[4]] + x6[[5]] + x6[[6]])/6


## ---- echo=TRUE----------------------------------------------------------
xm <- st_apply(x, c("x", "y"), mean)
all.equal(xm[[1]], x6$mean)

## ---- echo=TRUE----------------------------------------------------------
str(xm)


## ---- echo=TRUE----------------------------------------------------------
library(classInt)
args(classIntervals)


## ---- echo=TRUE----------------------------------------------------------
(cI <- classIntervals(olinda_sirgas2000$DEPRIV, n=7, style="fisher"))


## ---- echo=TRUE----------------------------------------------------------
library(RColorBrewer)
pal <- RColorBrewer::brewer.pal((length(cI$brks)-1), "Reds")
plot(cI, pal)


## ---- echo=TRUE----------------------------------------------------------
display.brewer.all()


## ---- echo=TRUE----------------------------------------------------------
plot(olinda_sirgas2000[,"DEPRIV"], breaks=cI$brks, pal=pal)


## ---- echo=TRUE----------------------------------------------------------
plot(olinda_sirgas2000[,"DEPRIV"], nbreaks=7, breaks="fisher", pal=pal)


## ---- echo=TRUE----------------------------------------------------------
mapview(olinda_sirgas2000, zcol="DEPRIV", col.regions=pal, at=cI$brks)


## ---- echo=TRUE----------------------------------------------------------
library(tmap)
tmap_mode("plot")
o <- tm_shape(olinda_sirgas2000) + tm_fill("DEPRIV", style="fisher", n=7, palette="Reds")
class(o)


## ---- echo=TRUE----------------------------------------------------------
o


## ---- echo=TRUE----------------------------------------------------------
o + tm_borders(alpha=0.5, lwd=0.5)


## ---- echo=TRUE----------------------------------------------------------
tmap_mode("view")


## ---- echo=TRUE----------------------------------------------------------
o + tm_borders(alpha=0.5, lwd=0.5)


## ---- echo=TRUE----------------------------------------------------------
tmap_mode("plot")


## ---- echo=TRUE, eval=FALSE----------------------------------------------
## tmaptools::palette_explorer()


## ---- echo=TRUE----------------------------------------------------------
library(cartography)
display.carto.all()


## ---- echo=TRUE----------------------------------------------------------
choroLayer(olinda_sirgas2000, var="DEPRIV", method="fisher-jenks", nclass=7, col=pal, legend.values.rnd=3)


## ---- echo=TRUE----------------------------------------------------------
library(ggplot2)


## ---- echo=TRUE----------------------------------------------------------
g <- ggplot(olinda_sirgas2000) + geom_sf(aes(fill=DEPRIV))
g


## ---- echo=TRUE----------------------------------------------------------
g + theme_void()


## ---- echo=TRUE----------------------------------------------------------
g + theme_void() + scale_fill_distiller(palette="Reds", direction=1)


## ---- echo=TRUE----------------------------------------------------------
library(colorspace)
hcl_palettes("sequential (single-hue)", n = 7, plot = TRUE)


## ---- echo=TRUE, eval=FALSE----------------------------------------------
## pal <- hclwizard()
## pal(6)


## ---- echo=TRUE----------------------------------------------------------
wheel <- function(col, radius = 1, ...)
  pie(rep(1, length(col)), col = col, radius = radius, ...) 
opar <- par(mfrow=c(1,2))
wheel(rainbow_hcl(12))
wheel(rainbow(12))
par(opar)

