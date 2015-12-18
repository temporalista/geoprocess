##### KRIGING INTERPOLATION AND SIMULATION #####

rm(list = ls())  # clean memory
graphics.off()   # close graphic windows

# Load libraries
library(foreign)
library(RColorBrewer)
library(rgdal)
library(sp)
library(maptools)
library(gstat)

# read DEM Data
dem=read.table("Sample_error.txt",header=TRUE)
dim(dem)
names(dem)
summary(dem)
class(dem)

# histogram
windows(width=5, height=5)
hist(dem$resid, main="DEM sample error")

# make spatial
coordinates(dem) = ~x+y

# read boundary study area
dem.poly = readShapePoly("study_area.shp")

# plot observations
windows(width=5, height=5)
spplot(dem, zcol="resid", xlim=c(152600,187620), ylim=c(418050,435150),
       key.space=list(x=0.02, y=0.01, corner=c(0,1)), cex=1, main="Resid data", 
       sp.layout=list("sp.polygons", dem.poly), col.regions=brewer.pal(4, "Oranges"))

# bubble plot
windows(width=5, height=5)
bubble(dem, zcol="resid", xlim=c(152600,187620), ylim=c(418050,435150),
       maxsize=3, key.entries=8*2^(1:6),
       sp.layout=list("sp.polygons", dem.poly))

# define gstat object and compute variogram
demres = gstat(id = c("resid"), formula = resid~1, data = dem)
vdemres = variogram(demres)
vdemres = variogram(demres, boundaries = c(200,500,750,1000,1250,1500,2000,2500,3000,3500))
windows(width = 5, height = 4)
plot(vdemres, plot.numbers = TRUE)


# fit nested variogram

vgdemres <- vgm(psill=1000, range=800, model="Sph", 
                add.to = vgm(nugget=3500, psill=4000, range=400, model="Exp"))
vgdemres
vgdemres <- fit.variogram(vdemres,vgdemres)
plot(vdemres,vgdemres,plot.numbers=FALSE)
vgdemres


#xvalid
demres.xv = krige.cv(resid~1, dem, vgdemres,nmax=20)
str(demres.xv)
hist(demres.xv$var1.pred)
names(demres.xv)
mean(demres.xv$zscore)
sd(demres.xv$zscore)

names(demres.xv)

windows(width = 4, height = 4)
hist(demres.xv$zscore)

windows(width = 5, height = 6)
bubble(demres.xv, zcol="residual", xlim=c(152600,187620), ylim=c(418050,435150),
       maxsize=3, key.entries=c(-330,-200,-100,-50,0,50,100,200,330), 
       sp.layout=list("sp.polygons", dem.poly), main = "cross-validation residuals")

ahn <- readGDAL("ahn100_f.asc")

# local kriging
demres.krig = krige(resid~1, dem, newdata = ahn, vgdemres, maxdist=3500, nmax=20)
names(demres.krig)

windows(width = 5, height = 6)
spplot(demres.krig, zcol = "var1.pred", col.regions = bpy.colors())

demres.krig$var1.sd = sqrt(demres.krig$var1.var)
windows(width = 5, height = 6)
spplot(demres.krig, zcol = "var1.sd", col.regions = bpy.colors())

#DEM+error
newdem<-ahn@data + demres.krig$var1.pred
coordinates(newdem) <- coordinates(ahn)
gridded(newdem) = TRUE
write.asciigrid(newdem, "DEMdefault.asc")
summary(ahn@data)
summary(newdem)

# simulation
demres.sim = krige(resid~1, dem, newdata = ahn, vgdemres, nsim = 100, nmax = 20, maxdist=3500)
names(demres.sim)
#windows(width = 5, height = 6)
#spplot(demres.sim, zcol = "sim1", col.regions = bpy.colors())
#windows(width = 5, height = 6)
#spplot(demres.sim, col.regions = bpy.colors())

#summary (demres.sim[3])
#summary (demres.sim[1])

#Cycles
outroot="M:/Flooding_case/sim100/dem"

for(i in 1:100)
{
  outname <- paste(outroot, i, ".asc", sep = "")
  outdem <- ahn@data + demres.sim[i]@data
  # promote data.frame to SpatialPPixelsDataFrame
  coordinates(outdem) <- coordinates(ahn)
  gridded(outdem) <- T
  write.asciigrid(outdem, outname)
}

#####INFLUX########
influx.sim <- rnorm(100, mean = 2700, sd = 120)
influx.sim
write.table(influx.sim,file="M:/Flooding_case/sim100/Influx100")



