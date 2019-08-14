
## ---- echo=TRUE----------------------------------------------------------
library(sf)
packageVersion("sf")


## ---- echo=TRUE----------------------------------------------------------
sf_extSoftVersion()


## ---- echo=TRUE----------------------------------------------------------
st_crs(22525)


## ---- echo=TRUE----------------------------------------------------------
cat(system("projinfo EPSG:22525", intern=TRUE), sep="\n")


## ---- echo=TRUE----------------------------------------------------------
cat(system("projinfo -s EPSG:22525 -t EPSG:31985", intern=TRUE), sep="\n")


## ---- echo=TRUE----------------------------------------------------------
olinda <- st_read("data/olinda.gpkg", quiet=TRUE)
st_crs(olinda)


## ---- echo=TRUE----------------------------------------------------------
xy_c <- st_centroid(st_geometry(olinda[  1,]))
st_coordinates(xy_c)


## ---- echo=TRUE----------------------------------------------------------
st_coordinates(st_transform(st_transform(xy_c, 4326), 31985))


## ---- echo=TRUE----------------------------------------------------------
# without CA7072_003.gsb
st_coordinates(st_transform(xy_c, 31985))


## ---- echo=TRUE, eval=FALSE----------------------------------------------
## # with CA7072_003.gsb
## st_coordinates(st_transform(xy_c, 31985))
## #          X       Y
## # 1 295489.3 9120352


## ---- echo=TRUE, eval=FALSE----------------------------------------------
## # with CA7072_003.gsb
## xy_c1 <- xy_c
## st_crs(xy_c1) <- "+proj=utm +zone=25 +south +ellps=intl +units=m +nadgrids=CA7072_003.gsb"
## print(st_coordinates(st_transform(xy_c1, 31985)), digits=9)
## #            X          Y
## # 1 295486.396 9120350.62


## ---- echo=TRUE, eval=FALSE----------------------------------------------
## # with CA7072_003.gsb
## cat(system(paste0("echo ", paste(xy, collapse=" "), " | cs2cs EPSG:22525 EPSG:31985"), intern=TRUE))
## # 295486.40 9120350.62 0.00


## ---- echo=TRUE----------------------------------------------------------
xy <- st_coordinates(xy_c)
# without CA7072_003.gsb
cat(system(paste0("echo ", paste(xy, collapse=" "), " | cs2cs EPSG:22525 EPSG:31985"), intern=TRUE))


## ---- echo=TRUE, warning=FALSE-------------------------------------------
# without CA7072_003.gsb
xy_c2 <- xy_c
st_crs(xy_c2) <- "+proj=utm +zone=25 +south +ellps=intl +units=m +towgs84=-206.05,168.28,-3.82,0,0,0,0"
st_coordinates(st_transform(xy_c2, 31985))


## ---- echo=TRUE----------------------------------------------------------
# without CA7072_003.gsb
# -DACCEPT_USE_OF_DEPRECATED_PROJ_API_H
st_coordinates(lwgeom::st_transform_proj(xy_c, 31985))


## ---- echo=TRUE----------------------------------------------------------
olinda <- st_read("output/olinda_sirgas2000.gpkg", quiet=TRUE)
xy_c <- st_centroid(st_geometry(olinda[  1,]))
st_coordinates(xy_c)


## ---- echo=TRUE----------------------------------------------------------
rgeos::version_GEOS0()


## ---- echo=TRUE, warning=FALSE-------------------------------------------
cV_old_default <- ifelse(rgeos::version_GEOS0() >= "3.7.2", 0L, FALSE)
yy <- rgeos::readWKT(readLines("data/invalid.wkt"))
rgeos::gIsValid(yy, byid=TRUE, reason=TRUE)


## ---- echo=TRUE----------------------------------------------------------
sf::sf_extSoftVersion()


## ---- echo=TRUE----------------------------------------------------------
sf::st_is_valid(sf::st_as_sf(yy), reason=TRUE)


## ---- echo=TRUE, warning=FALSE-------------------------------------------
ply <- rgeos::readWKT(readLines("data/ply.wkt"))
oo <- try(rgeos::gIntersection(yy, ply, byid=TRUE, checkValidity=cV_old_default), silent=TRUE)
print(attr(oo, "condition")$message)

## ---- echo=TRUE----------------------------------------------------------
ooo <- try(sf::st_intersection(sf::st_as_sf(yy), sf::st_as_sf(ply)), silent=TRUE)
print(attr(oo, "condition")$message)


## ---- echo=TRUE----------------------------------------------------------
cV_new_default <- ifelse(rgeos::version_GEOS0() >= "3.7.2", 1L, TRUE)
try(rgeos::gIntersection(yy, ply, byid=TRUE, checkValidity=cV_new_default), silent=TRUE)


## ---- echo=TRUE----------------------------------------------------------
oo <- rgeos::gIntersection(yy, ply, byid=TRUE, checkValidity=2L)
rgeos::gIsValid(oo)


## ---- echo=TRUE----------------------------------------------------------
oo <- rgeos::gIntersection(rgeos::gBuffer(yy, byid=TRUE, width=0), ply, byid=TRUE, checkValidity=1L)
rgeos::gIsValid(oo)


## ---- echo=TRUE----------------------------------------------------------
ooo <- sf::st_intersection(sf::st_buffer(sf::st_as_sf(yy), dist=0), sf::st_as_sf(ply))
all(sf::st_is_valid(ooo))


## ---- echo=TRUE----------------------------------------------------------
buildings <- sf::st_read("data/snow/buildings.gpkg", quiet=TRUE)
st_crs(buildings)


## ---- echo=TRUE----------------------------------------------------------
library(mapview)
mapview(buildings)


## ---- echo=TRUE----------------------------------------------------------
library(RSQLite)
db = dbConnect(SQLite(), dbname="data/snow/buildings.gpkg")
dbReadTable(db, "gpkg_spatial_ref_sys")$definition[4]
dbDisconnect(db)


## ---- echo=TRUE----------------------------------------------------------
buildings1 <- rgdal::readOGR("data/snow/buildings.shp", verbose=FALSE)
sp::proj4string(buildings1)


## ---- echo=TRUE----------------------------------------------------------
mapview(buildings1)


## ---- echo=TRUE, warning=FALSE-------------------------------------------
readLines("data/snow/buildings.prj")


## ---- echo=TRUE, warning=FALSE-------------------------------------------
fixed <- "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +nadgrids=OSTN15_NTv2_OSGBtoETRS.gsb +units=m +no_defs"
st_crs(buildings) <- fixed
sp::proj4string(buildings1) <- sp::CRS(fixed)


## ---- echo=TRUE----------------------------------------------------------
mapview(buildings)


## ---- echo=TRUE----------------------------------------------------------
mapview(buildings1)


## ---- echo=TRUE----------------------------------------------------------
library(sf)


## ---- echo=TRUE----------------------------------------------------------
olinda_sirgas2000 <- st_read("output/olinda_sirgas2000.gpkg", quiet=TRUE)
bounds <- st_sf(st_union(olinda_sirgas2000))
SG <- maptools::Sobj_SpatialGrid(as(bounds, "Spatial"), n=1000000)$SG


## ---- echo=TRUE----------------------------------------------------------
library(rgrass7)
packageVersion("rgrass7")
use_sp()
myGRASS <- "/home/rsb/topics/grass/g761/grass76"
myPROJSHARE <- "/usr/local/share/proj"
if (Sys.getenv("GRASS_PROJSHARE") == "") Sys.setenv(GRASS_PROJSHARE=myPROJSHARE)
loc <- initGRASS(myGRASS, tempdir(), SG=SG, override=TRUE)


## ---- echo=TRUE----------------------------------------------------------
execGRASS("g.mapset", mapset="PERMANENT", flag="quiet")
execGRASS("g.proj", flag="c", proj4=st_crs(bounds)$proj4string)
execGRASS("g.mapset", mapset=loc$MAPSET, flag="quiet")
execGRASS("g.region", flag="d")


## ---- echo=TRUE----------------------------------------------------------
execGRASS("r.in.gdal", flag=c("overwrite", "quiet"), input="output/elevation.tif", output="dem")
execGRASS("g.region", raster="dem")


## ---- echo=TRUE----------------------------------------------------------
execGRASS("r.watershed", flag=c("overwrite", "quiet"), elevation="dem", stream="stream", threshold=2500L, convergence=5L, memory=300L)
execGRASS("r.thin", flag=c("overwrite", "quiet"), input="stream", output="stream1", iterations=200L)


## ---- echo=TRUE----------------------------------------------------------
use_sf()
writeVECT(bounds, "bounds", v.in.ogr_flags=c("overwrite", "quiet"))
execGRASS("r.mask", vector="bounds", flag=c("overwrite", "quiet"))
execGRASS("r.to.vect", flag=c("overwrite", "quiet"), input="stream1", output="stream", type="line")
imputed_streams <- readVECT("stream", ignore.stderr=TRUE)


## ---- echo=TRUE, warning=FALSE-------------------------------------------
library(mapview)
mapview(imputed_streams)


## ---- echo=TRUE----------------------------------------------------------
execGRASS("r.slope.aspect", elevation="dem", slope="slope", aspect="aspect", flag=c("quiet", "overwrite"))
writeVECT(olinda_sirgas2000[, "SETOR_"], "olinda", ignore.stderr=TRUE, v.in.ogr_flags=c("overwrite", "quiet"))
execGRASS("v.rast.stats", map="olinda", raster=c("slope", "aspect"), method=c("first_quartile", "median", "third_quartile"), column_prefix=c("slope", "aspect"), flag=c("c", "quiet"))


## ---- echo=TRUE----------------------------------------------------------
execGRASS("r.in.gdal", flag=c("overwrite", "quiet"), input="output/L7_ndvi.tif", output="ndvi")
execGRASS("g.region", raster="ndvi")
execGRASS("v.rast.stats", map="olinda", raster="ndvi", method=c("first_quartile", "median", "third_quartile"), column_prefix="ndvi", flag=c("c", "quiet"))


## ---- echo=TRUE----------------------------------------------------------
olinda_gmm_ndvi <- readVECT("olinda", ignore.stderr=TRUE)
head(olinda_gmm_ndvi)


## ----echo=FALSE----------------------------------------------------------
knitr::include_graphics('snowmap.png')


## ----echo=FALSE----------------------------------------------------------
knitr::include_graphics('brodyetal00_fig1.png')


## ---- echo=TRUE----------------------------------------------------------
library(sf)
bbo <- st_read("data/snow/bbo.gpkg")
fixed <- "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +nadgrids=OSTN15_NTv2_OSGBtoETRS.gsb +units=m +no_defs"
st_crs(bbo) <- fixed


## ---- echo=TRUE----------------------------------------------------------
library(rgrass7)
myPROJSHARE <- "/usr/local/share/proj"
if (Sys.getenv("GRASS_PROJSHARE") == "") Sys.setenv(GRASS_PROJSHARE=myPROJSHARE)
myGRASS <- "/home/rsb/topics/grass/g761/grass76"
td <- tempdir()
SG <- maptools::Sobj_SpatialGrid(as(bbo, "Spatial"))$SG
use_sp()
soho <- initGRASS(gisBase=myGRASS, home=td, SG=SG, override=TRUE)
soho


## ---- echo=TRUE----------------------------------------------------------
MAPSET <- execGRASS("g.mapset", flags="p", intern=TRUE)
execGRASS("g.mapset", mapset="PERMANENT", flags="quiet")
execGRASS("g.proj", flags=c("p", "quiet"))
execGRASS("g.proj", proj4=st_crs(bbo)$proj4string, flags=c("c", "quiet"))


## ---- echo=TRUE----------------------------------------------------------
execGRASS("g.mapset", mapset=MAPSET, flags="quiet")
execGRASS("g.region", flags="p", intern=TRUE)[3:11]
execGRASS("g.region", flags="a", res="1")
execGRASS("g.region", flags="p", intern=TRUE)[3:11]


## ---- echo=TRUE, warning=FALSE-------------------------------------------
buildings <- st_read("data/snow/buildings.gpkg", quiet=TRUE)
st_crs(buildings) <- fixed
deaths <- st_read("data/snow/deaths.gpkg", quiet=TRUE)
st_crs(deaths) <- fixed
sum(deaths$Num_Css)
b_pump <- st_read("data/snow/b_pump.gpkg", quiet=TRUE)
st_crs(b_pump) <- fixed
nb_pump <- st_read("data/snow/nb_pump.gpkg", quiet=TRUE)
st_crs(nb_pump) <- fixed


## ---- echo=TRUE, warning=FALSE-------------------------------------------
use_sf()
fl <- c("overwrite", "quiet")
writeVECT(bbo, vname="bbo", v.in.ogr_flags=c("o", fl), ignore.stderr=TRUE)
writeVECT(buildings[,1], vname="buildings", v.in.ogr_flags=c("o", fl), ignore.stderr=TRUE)
writeVECT(b_pump, vname="b_pump", v.in.ogr_flags=c("o", fl), ignore.stderr=TRUE)
writeVECT(nb_pump, vname="nb_pump", v.in.ogr_flags=c("o", fl), ignore.stderr=TRUE)
writeVECT(deaths, vname="deaths", v.in.ogr_flags=c("o", fl), ignore.stderr=TRUE)
execGRASS("g.list", type="vector", intern=TRUE)


## ---- echo = TRUE, mysize=TRUE, size='\\tiny'----------------------------
execGRASS("v.overlay", ainput="buildings", binput="bbo", operator="xor", output="roads", flags=fl, ignore.stderr = TRUE)
execGRASS("v.to.rast", input="roads", output="rroads", use="val", value=1, flags=fl)
execGRASS("r.stats", input="rroads", flags=c("c", "quiet"))


## ---- echo = TRUE, mysize=TRUE, size='\\tiny'----------------------------
execGRASS("r.buffer", input="rroads", output="rroads4", distances=4, flags=fl)
execGRASS("r.stats", input="rroads4", flags=c("c", "quiet"))
tf <- tempfile()
cat("1 2 = 1\n", file=tf)
execGRASS("r.reclass", input="rroads4", output="rroads4a", rules=tf, flags=fl)
execGRASS("r.stats", input="rroads4a", flags=c("c", "quiet"))


## ---- echo = TRUE, mysize=TRUE, size='\\tiny'----------------------------
execGRASS("r.cost", input="rroads4a", output="dist_broad", start_points="b_pump", flags=fl)
execGRASS("r.cost", input="rroads4a", output="dist_not_broad", start_points="nb_pump", flags=fl)


## ---- echo = TRUE, mysize=TRUE, size='\\tiny'----------------------------
execGRASS("v.db.addcolumn", map="deaths", columns="broad double precision", flags="quiet")
execGRASS("v.what.rast", map="deaths", raster="dist_broad", column="broad", flags="quiet")
execGRASS("v.db.addcolumn", map="deaths", columns="not_broad double precision", flags="quiet")
execGRASS("v.what.rast", map="deaths", raster="dist_not_broad", column="not_broad", flags="quiet")


## ---- echo = TRUE, mysize=TRUE, size='\\tiny'----------------------------
deaths1 <- readVECT("deaths", ignore.stderr=TRUE)
deaths1$b_nearer <- deaths1$broad < deaths1$not_broad
by(deaths1$Num_Css, deaths1$b_nearer, sum)


## ---- echo=TRUE, warning=FALSE-------------------------------------------
water_buf_50 <- st_buffer(imputed_streams, dist=50)
setor_area <- st_area(olinda_sirgas2000)
near_water0 <- st_intersection(olinda_sirgas2000[,"SETOR_"], water_buf_50[,"cat"])
near_water <- aggregate(near_water0, by=list(near_water0$SETOR_), head, n=1)


## ---- echo=TRUE----------------------------------------------------------
area_near_water <- st_area(near_water)
olinda_sirgas2000$setor_area <- setor_area
o <- match(near_water$SETOR_, olinda_sirgas2000$SETOR_)
olinda_sirgas2000$area_near_water <- 0
olinda_sirgas2000$area_near_water[o] <- area_near_water
olinda_sirgas2000$prop_near_water <- olinda_sirgas2000$area_near_water/olinda_sirgas2000$setor_area
summary(olinda_sirgas2000$prop_near_water)


## ---- echo=TRUE----------------------------------------------------------
library(tmap)
tm_shape(olinda_sirgas2000) + tm_fill("prop_near_water", palette="Blues", style="fisher", n=5)


## ---- echo=TRUE, warning=FALSE-------------------------------------------
streams_by_setor <- st_intersection(olinda_sirgas2000[,"SETOR_"], imputed_streams[,"cat"])
lngths0 <- aggregate(streams_by_setor, by=list(streams_by_setor$SETOR_), head, n=1)
lngths <- st_length(lngths0)
units(lngths) <- "mm"
o <- match(lngths0$SETOR_, olinda_sirgas2000$SETOR_)
olinda_sirgas2000$lngths <- 0
olinda_sirgas2000$lngths[o] <- lngths
olinda_sirgas2000$lngth_area <- olinda_sirgas2000$lngths/olinda_sirgas2000$setor_area


## ---- echo=TRUE----------------------------------------------------------
tm_shape(olinda_sirgas2000) + tm_fill("lngth_area", palette="Blues", style="fisher", n=5)


## ---- echo=TRUE----------------------------------------------------------
cor(olinda_sirgas2000$lngth_area, olinda_sirgas2000$prop_near_water)


## ---- echo=TRUE----------------------------------------------------------
library(raster)
r <- raster("output/elevation.tif")
r


## ---- echo=TRUE, cache=TRUE----------------------------------------------
slope_aspect <- terrain(r, opt=c('slope','aspect'), unit='degrees', neighbors=8)


## ---- echo=TRUE, cache=TRUE----------------------------------------------
slopes <- extract(slope_aspect, olinda_sirgas2000, fun=median)


## ---- echo=TRUE, cache=TRUE----------------------------------------------
ndvi <- extract(raster("output/L7_ndvi.tif"), olinda_sirgas2000, fun=median)


## ---- echo=TRUE----------------------------------------------------------
summary(ndvi[,1])
if (exists("olinda_gmm_ndvi")) summary(olinda_gmm_ndvi$ndvi_median)


## ---- echo=TRUE----------------------------------------------------------
summary(slopes[,1])
if (exists("olinda_gmm_ndvi")) summary(olinda_gmm_ndvi$slope_median)


## ---- echo=TRUE----------------------------------------------------------
summary(slopes[,2])
if (exists("olinda_gmm_ndvi")) summary(olinda_gmm_ndvi$aspect_median)


## ---- echo=TRUE, warning=FALSE-------------------------------------------
library(sf)
buildings1 <- st_intersection(buildings, bbo)
buildings2 <- st_buffer(buildings1, dist=-4)


## ---- echo=TRUE, warning=FALSE-------------------------------------------
library(mapview)
mapview(buildings2)


## ---- echo=TRUE----------------------------------------------------------
library(raster)
resolution <- 1
r <- raster(extent(buildings2), resolution=resolution, crs=fixed)
r[] <- resolution
summary(r)


## ---- echo=TRUE, cache=TRUE----------------------------------------------
buildings3 <- as(buildings2[!st_is_empty(buildings2),], "Spatial")
cfp <- cellFromPolygon(r, buildings3)
is.na(r[]) <- unlist(cfp)
summary(r)


## ---- echo=TRUE, warning=FALSE-------------------------------------------
library(mapview)
mapview(r)


## ---- echo=TRUE, warning=FALSE, message=FALSE----------------------------
library(gdistance)


## ---- echo=TRUE, cache=TRUE----------------------------------------------
tr1 <- transition(r, transitionFunction=function(x) 1/mean(x), directions=8, symm=TRUE)


## ---- echo=TRUE, cache=TRUE----------------------------------------------
sp_deaths <- as(deaths, "Spatial")
d_b_pump <- st_length(st_as_sfc(shortestPath(tr1, as(b_pump, "Spatial"), sp_deaths, output="SpatialLines")))


## ---- echo=TRUE, cache=TRUE----------------------------------------------
res <- matrix(NA, ncol=nrow(nb_pump), nrow=nrow(deaths))
sp_nb_pump <- as(nb_pump, "Spatial")
for (i in 1:nrow(nb_pump)) res[,i] <- st_length(st_as_sfc(shortestPath(tr1, sp_nb_pump[i,], sp_deaths, output="SpatialLines")))
d_nb_pump <- apply(res, 1, min)


## ---- echo=TRUE----------------------------------------------------------
library(units)
units(d_nb_pump) <- "m"
deaths$b_nearer <- d_b_pump < d_nb_pump
by(deaths$Num_Css, deaths$b_nearer, sum)

