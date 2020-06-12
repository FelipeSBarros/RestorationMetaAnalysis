# rm(list=ls())

sp.corrAnalysis <- function(r.model, npts, ...){

  # loading library
  library(raster)
  library(readr)
  library(pgirmess)
  
  # output drectory
  outDir <- "./CorrelationResults"
  
  # testign if output dir exists
  # if not, will be created
  if (! dir.exists(outDir)){
    dir.create(outDir)
  }
  
  # generating random points
  df.pts <- as_tibble(
    sampleRandom(
      r.model, npts*5, xy = TRUE, sp = FALSE, na.rm = TRUE))
  
  # getting model name to use when saving as csv
  modelName <- colnames(df.pts)[ncol(df.pts)]
  
  # saving as csv
  write_csv(df.pts, paste0(outDir, "/", modelName, "_rndPts.csv"))
  
  # estimating correlog
  modelCorrelog <- pgirmess::correlog(coords = df.pts[,c(1:2)], 
                          z = df.pts[[ncol(df.pts)]])
 
  # saving correlog result as csv
  write_csv(as_tibble(modelCorrelog), paste0(outDir, "/", modelName, "_moranCorrelog.csv"))

  # saving correlog plot
  png(paste0(outDir, "/", modelName, "_moranCorrelog.png"))
  plot(modelCorrelog)
  dev.off()
  }

getDistance <- function(inputDir, modelName, suffix, ...){
  # inputDir <- "./CorrelationResults"
  # r.model1 = stack( list.files("../Data/Raster/Results", pattern = "Abundance_Norm.tif$|ness_Norm.tif$", recursive = TRUE, full.names = TRUE))
  # modelName <- r.model1@layers[[1]]@data@names
  # modelName1 <- r.model1@layers[[1]]@data@names
  # modelName2 <- r.model2@layers[[1]]@data@names
  # suffix = "_moranCorrelog.csv"
  
  # loading library
  library(readr)
  
  # loading spatial correlation data
  distCorr <- read_csv(
    list.files(
      inputDir, paste0(modelName, suffix), full.names = TRUE)
    )
  
  # identifying min dist w/out sp correlation
  distCorr <- distCorr[ min( which(distCorr$coef < 0)), 1 ]
  return(distCorr)
}

fcn.name <- function(r.model1, r.model2, rmhmodel = FALSE){
  # r.model1 = stack( list.files("../Data/Raster/Results", pattern = "Abundance_Norm.tif$|ness_Norm.tif$", recursive = TRUE, full.names = TRUE)[1])
  # r.model2 = stack( list.files("../Data/Raster/Results", pattern = "Abundance_Norm.tif$|ness_Norm.tif$", recursive = TRUE, full.names = TRUE)[2])
  
  # loading library
  library(dplyr)
  library(spatstat)
  
  # defining input dir and model name
  inputDir <- "./CorrelationResults"
  modelName1 <- r.model1@layers[[1]]@data@names
  modelName2 <- r.model2@layers[[1]]@data@names
  
  # loading spatial correlation data
  distCorr1 <- getDistance("./CorrelationResults", r.model1@layers[[1]]@data@names,  "_moranCorrelog.csv")
  distCorr2 <- getDistance("./CorrelationResults", r.model2@layers[[1]]@data@names,  "_moranCorrelog.csv")
  
  # identifying max dist w/out sp correlation
  distCorrFinal <- max(distCorr1, distCorr2)
  
  # laoding world shapelie
  library(sf)
  library(spatstat)
  
  # world <- read_sf("/media/felipe/DATA/Proyectos/Dados/IPUMSI_world_release2017/world_countries_2017.shp")
  
 # plot(
 #  st_bbox(
 #   st_geometry(world)))
  world.w <- owin(xrange=c(-122.1189, 180), yrange=c(-34.99998, 35), poly=NULL, mask=NULL,
       unitname=NULL, xy=NULL)
  
  # Create an owin:
  #world.w <- as.owin(
  #  as_Spatial(
  #    st_bbox(
  #      world)))
  # or
  #world.w <- as(as_Spatial(world),"owin")
  # plot(world.w)
  
  betaModel <- 10 # what is beta on hardcore model
  t <- distCorrFinal/1
  if ( rmhmodel ){
    model <- rmhmodel(cif = "hardcore", 
                      par = list(beta = betaModel, hc = t), 
                      w = world.w)
    X2 <- rmh(model)
    #> Checking arguments..determining simulation windows...Starting simulation.
    #> Initial state...Ready to simulate. Generating proposal points...Running Metropolis-Hastings.
    
  }
  
  
  # plot(X2, main = paste("Approx. sim. of hardcore model; beta =", betaModel, "and R =", t))
  # plot(r.model1, add=T)
  
  
  #gerando os ptos regulares considerando os pontos aleatoriogerados anteriormente 'random'
  # random <- rpoint(1000, giveup=2000, verbose=FALSE, nsim=1, drop=TRUE)
  # plot(random)
  
  inibition_random <- rSSI( r = distCorrFinal/10, 
                            n = 100, 
                            giveup = 10, 
                            nsim=1, 
                            drop=TRUE,
                            win = world.w)
  # plot(inibition_random)
  # plot(r.model1, add=T)
  return(inibition_random)
  }

world <- read_sf("/media/felipe/DATA/Proyectos/Dados/IPUMSI_world_release2017/world_countries_2017.shp")

fcn.name <- function(r.model1, r.model2, distance = 600000, meters = TRUE, windowSample = world, nSample = 2000){
  # r.model1 = stack( list.files("../Data/Raster/Results", pattern = "Abundance_Norm.tif$|ness_Norm.tif$", recursive = TRUE, full.names = TRUE)[1])
  # r.model2 = stack( list.files("../Data/Raster/Results", pattern = "Abundance_Norm.tif$|ness_Norm.tif$", recursive = TRUE, full.names = TRUE)[2])
  
  # loading library
  library(dplyr)
  library(spatstat)
  
  # changing the distance to decimal degrees
  if( meters){
    convertFactor = 111000
    distance = round(distance / convertFactor, 2)
  }
  
  # defining window to sample
  world.w <- owin(xrange = extent(
    windowSample)[1:2], 
    yrange = extent(
      windowSample)[3:4])
  
  # sampling random points with hinibition
  inibition_random <- rSSI( r = distance, 
                            n = nSample,
                            win = world.w)
  
  plot(inibition_random)
  
  inibition_random <- st_as_sf(inibition_random)
  inibition_random <- st_set_crs(inibition_random, 4326)
  
  library(sp)
  class(inibition_random)
  as(inibition_random, "Spatial")
  sf::as_Spatial(inibition_random$geom)
  class(st_cast(inibition_random)  )
  class(st_sfc(st_geometry(inibition_random)))
  class(st_sf(st_geometry(inibition_random)))
  
  raster::extract(r.model1[[1]], inibition_random)
  raster::extract(r.model1[[1]], st_sfc(st_geometry(inibition_random)))
  
  betaModel <- 10 # what is beta on hardcore model
  if ( rmhmodel ){
    model <- rmhmodel(cif = "hardcore", 
                      par = list(beta = betaModel, hc = distance), 
                      w = world.w)
    X2 <- rmh(model)
    }
  }
