###### Gippet et al., 2023,   
# SPECIES distirbution maps, human influence metrics and supplementary figures 

library(biomod2) #for making sdm maps
library(raster) #dealing with plotting and manipulating maps
library(ggplot2) #human density plots
library(spThin) #thinning function
library(RStoolbox) #rasterPCA function
library(climateStability)
library(viridis) #colour palette 
library(dunn.test) # pairwise dunn test 

##### #### #####
## load data 
##### #### #####

# load data 
load('Lfulica2023.RData')
#save(Lfulica , background, Human_density, bin, currentPred_ensem, w_rworldmap, file='Lfulica2023.RData')
## dataset explanation
    # Lfulica  - dataset for all L. fulica occurence points
    # background - full dataset of XXX species of XXX order for background pseudoabsence selection
    # human density - human density map downloaded from XX 
    # bin - binary output map from SDM code showing areas suitable for Lfulica invasion (binary)
    # currentPred_ensem - projection map output from biomod showing suitable areas for Lfulica invasion (raw)
    #w_rworldmap - worldmap for plots
    
#thin L.fulica dataset 
native <- Lfulica[which(Lfulica$range=='N'),]
native <-  thin(
  native,
  lat.col = "decimallatitude",
  long.col = "decimallongitude",
  spec.col = "species",
  thin.par=5,
  reps=1,
  locs.thinned.list.return = T,
  write.files = F,
  write.log.file = F,
)
native <- as.data.frame(native)
native$range <- 'N'

invasive <- Lfulica[which(Lfulica $range=='I'),]
invasive <-  thin(
  invasive, 
  lat.col = "decimallatitude",
  long.col = "decimallongitude",
  spec.col = "species",
  thin.par=5,
  reps=1,
  locs.thinned.list.return = T,
  write.files = F,
  write.log.file = F,
)
invasive <- as.data.frame(invasive)
invasive$range <- 'I'
Lfulica <- rbind(native, invasive)
colnames(Lfulica) <- c( 'decimallongitude', 'decimallatitude', "range")


worldclim_2_5 <- raster::getData("worldclim", var = "bio", res = 2.5)
pca_Worldclim <- rasterPCA(worldclim_2_5 ,spca=TRUE, nComp=5, maskCheck=F)
pca_worldclim  <- stack(pca_Worldclim$map)
myExpl <-  pca_worldclim 


#### ##### ###
#### FORMATING DATASET ####
#### #### ####

## Invasive + native dataset 
myRespXY <- Lfulica[c('decimallongitude', 'decimallatitude')]# the XY coordinates of species data
myResp <- rep(1, nrow(myRespXY ) ) # the presence/absences data for our species

#background points based on the known distirbution of other species 
PA <- background[c('Longitude', 'Latitude')]
colnames(PA) <- c('decimallongitude', 'decimallatitude')
myRespXY.PA <- rbind(myRespXY , PA)
myRespPA <- rep(NA, nrow(PA ) ) # the presence/absences data for our species
myResp.PA <- c(myResp, myRespPA ) 
## Next, format pseudoabsenses to sample for PA.table argument. 
myPAtable <- data.frame(PA1 = ifelse(is.na(myResp.PA), FALSE, TRUE),
                        PA2 = ifelse(is.na(myResp.PA), FALSE, TRUE), 
                        PA3 = ifelse(is.na(myResp.PA), FALSE, TRUE), 
                        PA4 = ifelse(is.na(myResp.PA), FALSE, TRUE), 
                        PA5 = ifelse(is.na(myResp.PA), FALSE, TRUE), 
                        PA6 = ifelse(is.na(myResp.PA), FALSE, TRUE), 
                        PA7 = ifelse(is.na(myResp.PA), FALSE, TRUE), 
                        PA8 = ifelse(is.na(myResp.PA), FALSE, TRUE), 
                        PA9 = ifelse(is.na(myResp.PA), FALSE, TRUE), 
                        PA10 = ifelse(is.na(myResp.PA), FALSE, TRUE)
)
#included in analysis 'making them true' )
for (i in 1:ncol(myPAtable)) myPAtable[sample(which(myPAtable[, i] == FALSE), nrow(myRespXY)), i] = TRUE  # tells to sample the length of the presence dataset
# then can format the biomod data 
myBiomodData_Inv <- BIOMOD_FormatingData(resp.var = myResp.PA, #shows whether occurence points are presences (1) or pseudoabsences (NA)
                                         expl.var = myExpl, #explanatory, the stack of raster layers from biomod. 
                                         resp.xy = myRespXY.PA, #occurence points
                                         resp.name = spe, #species name
                                         PA.strategy = 'user.defined',
                                         PA.user.table = myPAtable
) 


##### #### #####
### RUNNING SDM ### 
##### #### #####
myBiomodOption <- BIOMOD_ModelingOptions() # default parameters 
### Computing the models
myBiomodModelOut_E <- BIOMOD_Modeling(
  myBiomodData_Inv, #data just formatted - introduced 
  models = c('GLM','GBM','CTA','ANN','MARS','RF', 'MAXENT.Phillips.2'), 
  bm.options =  myBiomodOption,
  nb.rep = 10, #number of repeats for each model 
  data.split.perc = 70, #testing-evaluation split 
  data.split.table = NULL,
  do.full.models = TRUE,
  weights = NULL,
  prevalence = 0.5,
  metric.eval = c( "TSS"), #evaluation metric
  var.import = 3,
  save.output = TRUE,
  scale.models = FALSE,
  nb.cpu = 1,
  seed.val = NULL,
  do.progress = TRUE,
  modeling.id = paste(species_name,"Allrange_E",sep=""))
###  ensemble models - this is where we can get the model outputs out for 
myBiomodEM2_E <- BIOMOD_EnsembleModeling(
  bm.mod =  myBiomodModelOut_E,
  models.chosen = 'all',
  em.by='all',
  metric.eval = c('TSS'), #evaluation model used for the threshold. 
  metric.select.thresh = c(0.7), #threshold to exclude the model for TSS . 
  prob.mean = F,
  prob.cv = F,
  prob.ci = F,
  prob.ci.alpha = 0.05,
  prob.median = T,
  committee.averaging = F,
  prob.mean.weight = T,
  prob.mean.weight.decay = 'proportional')
# get evaluation scores
myBiomodModelEval_ensemble_invnat <-get_evaluations(myBiomodEM2_E) 

### projection of the individual models 
myBiomodProj_E <- BIOMOD_Projection(
  bm.mod = myBiomodModelOut_E,
  new.env = myExpl,
  proj.name = 'invasive', 
  models.chosen = 'all',
  metric.binary  = 'TSS',
  compress = T,
  build.clamping.mask = F, 
  output.format = '.grd', 
  on_0_1000=FALSE
)
# see predicitons 
myCurrentProj_E <- get_predictions(myBiomodProj_E)

### Projection of  ensemble model
myBiomodProjPresEnsemble_E <- BIOMOD_EnsembleForecasting(
  bm.em = myBiomodEM2_E,
  bm.proj = myBiomodProj_E, 
  metric.binary= c('TSS'))

currentPred_ensem<- get_predictions(myBiomodProjPresEnsemble_E) #map preictions
presentEnsembleResult_E_median <- currentPred_ensem[[1]] #This is the median model ensemble
presentEnsembleResult_E_wMean <- currentPred_ensem[[2]]


##### #### #####
## Figure 2A - plot the Projection
##### #### #####
currentPred_ensem <- projectRaster(currentPred_ensem, crs = '+proj=moll')
map <-currentPred_ensem[[1]] # plot the weighted mean
raster01 = function(r){
            # get the min max values
            minmax_r = range(values(r), na.rm=TRUE) 
            # rescale 
            return( (r-minmax_r[1]) / (diff(minmax_r)))
            }
map <- raster01(map)
map.p <- rasterToPoints(map)
df <- data.frame(map.p)
colnames(df) <- c("Longitude", "Latitude", "MAP")
df$MAP <- df$MAP/100
fill1 <- viridis_pal(option="mako", direction = -1, begin =0.1, end=0.9)(30)
w_rworldmap_sf <- st_as_sf(w_rworldmap)      

Fig2A <- tm_shape(w_rworldmap_sf, projection=st_crs('+proj=moll'), raster.warp = FALSE) + #fill="grey95", col="grey75",
          tm_polygons() + 
          tm_layout(title = "",
                    frame=F,
                    bg.color = "#f1fefcff",
                    earth.boundary=T,
                    earth.boundary.color = "black",
                    space.color = "white",
                    legend.outside = TRUE,
                    legend.outside.position = "right",
                    earth.boundary.lwd=0.1) +
          tm_shape(map , projection=st_crs('+proj=moll'), raster.warp=F) + #fill="grey95", col="grey75",
          tm_raster(palette=fill1, style = 'cont', 
                    breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                    title='Probability of \nOccurence') 
tmap_save(filename="map_fig2.pdf", tm=Fig2A, width=7.5, height= 6, units = "in")

##### #### #####
## Figure 2A - Human density overlaps 
##### #### #####

background$range <- "B"
background2 <- background[sample(nrow(background), size=90000), ] #downsampling
lfulica <- Lfulica[,c(3,4, 7)]
colnames( background2) <- c('x', "decimallongitude", "decimallatitude",  "range" )

#combine datasets
all <- rbind(lfulica, background2[,c(2:4)])
xy <- all[,1:2]
points <- SpatialPointsDataFrame(xy, data=all, proj4string = r@crs) #this converts the lat and longs into spatial data
values <- raster::extract(Human_desntiy,points)
points$values <- values
#extract human density for each occurence points 
points$constrained <- raster::extract(bin, points)
points <- points %>%  as.data.frame()


### plot an overlap map of the native and invasive range 
Fig2B <- ggplot(points[which(points$constrained==T),] , aes(y=values + 0.1, x=range, colour=range, fill=range)) + 
  ggdist::stat_halfeye(adjust = .5, width = .7, .width = 0, justification = -.2, point_colour = NA, alpha=0.7) + 
  geom_jitter(width = .1, alpha = .2) + 
  geom_boxplot(width = .08, outlier.shape = NA, alpha=0.7, fill='white', colour='black') + 
  xlab('Range') +
  ylab('Human Density') + 
  scale_x_discrete(labels=c("B" = "Background", "I" = "Invasive", "N" = "Native")) + 
  scale_y_continuous(trans='log10', breaks=c(1, 10, 100, 1000, 10000)) + 
  scale_colour_manual(values=c('darkgrey', "darkorchid4", "darkcyan"), labels=c('Background', 'Invasive', 'Native')) +
  scale_fill_manual(values=c('darkgrey', "darkorchid4", "darkcyan"), labels=c('Background', 'Invasive', 'Native')) +
  theme_classic()  + theme(legend.position='none') +
  coord_flip()
ggsave(file="rainplot.pdf", plot=Fig2B, width=7.5, height= 7.5, units = "in", device='pdf')

##statistically test 
hist(points$values[which(points$constrained==T)])
test <- points[which(points$constrained==T),]
kruskal.test(values ~ range, data = test) 
dunn.test(test$values,  as.factor(test$range), method='BH')

##minimaps 
##### ##### ##### ##### #####
### plot points  of ranges 
##### ##### ##### ##### #####
#native range
tm_shape(World, projection=st_crs('+proj=moll')) + #fill="grey95", col="grey75",
  tm_polygons() +
  tm_layout(title = "",
            frame=F,
            bg.color = "#f1fefcff",
            earth.boundary=T,
            earth.boundary.color = "black",
            space.color = "white",
            earth.boundary.lwd=1.5) + 
  tm_shape(sf_N) + 
  tm_dots(col = "darkcyan", alpha=0.8, size=0.08) 
#invasive range 
tmap_mode("view")
tmap_mode("plot") 
tm_shape(World, projection=st_crs('+proj=moll')) + #fill="grey95", col="grey75",
  tm_polygons() +
  tm_layout(title = "",
            frame=F,
            bg.color = "#f1fefcff",
            earth.boundary=T,
            earth.boundary.color = "black",
            space.color = "white",
            earth.boundary.lwd=1.5) + 
  tm_shape(sf_E) + 
  tm_dots(col="darkorchid4", alpha=0.8, size=0.08)
#background points
setwd("/Volumes/RECHERCHE/FAC/FBM/DEE/cbertels/default/D2c/obates/Data")
bin<- raster("L.fulica/proj_invasive/proj_invasive_L.fulica_TSSbin.grd") 
background  #all background points 
xy <- background[,3:4]
points <- SpatialPointsDataFrame(xy, data=background[,-(3:4)], proj4string = r@crs) 
values2 <- raster::extract(bin, points)
points$constrained <- values2
points <- points[which(points$constrained==T),]


tm_shape(World, projection=st_crs('+proj=moll')) + #fill="grey95", col="grey75",
  tm_polygons() +
  tm_layout(title = "",
            frame=F,
            bg.color = "#f1fefcff",
            earth.boundary=T,
            earth.boundary.color = "black",
            space.color = "white",
            earth.boundary.lwd=1.5) + 
  tm_shape(points) + 
  tm_dots(col="black", alpha=0.8, size=0.08)





### figS1 ###
##plot the worldclim maps 
setwd("/Volumes/RECHERCHE/FAC/FBM/DEE/cbertels/default/D2c/obates/Data")
datafiles <- Sys.glob("worldclim_pca_2_5/*.tif") 
pca_world <- stack()
for(i in 1:NROW(datafiles)){
  tempraster <- raster(datafiles[i])
  pca_world <- stack(pca_world,tempraster)
}
worldclim_2_5 <- raster::getData("worldclim", var = "bio", res = 2.5) # to get the right projection
crs( pca_world) <- crs(worldclim_2_5 )

pca_world$PCA_1 %>% values %>% hist
pca_world$PCA_2 %>% values %>% hist
pca_world$PCA_3 %>% values %>% hist

plot_list = list()
limits <- list()
limits[[1]] <- c(0,0.7)
limits[[2]] <- c(0.5,1)
limits[[3]] <- c(0.6,1)
limits[[4]] <- c(0,0.4)
limits[[5]] <- c(0.5,0.9)

for(i in 1:5) {
  plot <- pca_world[[i]]
  plot <- raster01(plot)
  
  map.p <- rasterToPoints( plot)
  df <- data.frame(map.p)
  colnames(df) <- c("Longitude", "Latitude", "MAP")
  #hist(df$MAP)
  map <-  ggplot(data=df, aes(y=Latitude, x=Longitude)) +
    geom_tile(aes(fill=(MAP))) +
    coord_map(
      projection = "moll") +
    # coord_map('+proj=moll')
    theme_void() +
    #coord_equal() +
    scale_fill_continuous(limit=limits[[i]], type="viridis") +  #I LIKE THIS BEST 
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x =  element_blank(),
          axis.text.y =  element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "right",
          legend.key = element_blank(),
          legend.title = element_blank() ) + 
    theme(text = element_text(family = "Helvetica", size=10), title=element_text(face="italic")) 
  
  plot_list[[i]] <- map 
}

figS1 <- ggarrange(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]], plot_list[[5]],scatter ,  labels=c('A', 'B', "C", "D", 'E', 'F'))
ggsave( figS1 , file="figS1 .pdf",  width=7.5, height= 7.5, units = "in", device='pdf')




##### ##### ##### #####
## plot FigS3
##### ##### ##### #####
load("assesmentsSnails.RData")
load("~/Desktop/Niche_shift/Worldmap_2019-04-03.RData")
Lfulica <- read.csv('occ_Lfulica2.csv')
nrow(Lfulica[which(Lfulica$range=='I'),])
nrow(Lfulica[which(Lfulica$range=='N'),])
background <- read.csv('background_25mm.csv')
nrow(background)
#background <- read.csv('occ_background2.csv')
tmap_mode("plot")
## L.fulcia map 
# prepare data 
xy <- Lfulica[,c(3:4)]
spdf <- SpatialPointsDataFrame(xy,  Lfulica, proj = r@crs)
sf <- st_as_sf(spdf, "+proj=moll")
sf_N <- sf[sf$range == "N",] # native range
sf_E <- sf[sf$range == "I",] #exotic range 
data("World")
# make map 
lful <- tm_shape(World, projection=st_crs('+proj=moll')) + #fill="grey95", col="grey75",
  tm_polygons() +
  tm_layout(title = "",
            frame=F,
            bg.color = "#f1fefcff",
            earth.boundary=T,
            earth.boundary.color = "black",
            space.color = "white",
            earth.boundary.lwd=1.5) + 
  tm_shape(sf_N) +  # add the native points
  tm_dots(col = "darkcyan", alpha=0.5, size=0.08) + 
  tm_shape(sf_E) +  # add the introduced points 
  tm_dots(col="darkorchid4", alpha=0.5, size=0.08) 

## background map 
# prepare data 
xy <- background[,c(2:3)]
spdf <- SpatialPointsDataFrame(xy,  background, proj = r@crs)
sf <- st_as_sf(spdf)
plot(sf)
# make map 
back  <- tm_shape(World, projection=st_crs('+proj=moll')) + #fill="grey95", col="grey75",
  tm_polygons() +
  tm_layout(title = "",
            frame=F,
            bg.color = "#f1fefcff",
            earth.boundary=T,
            earth.boundary.color = "black",
            space.color = "white",
            earth.boundary.lwd=1.5) + 
  tm_shape(sf) + 
  tm_dots(col = "black", alpha=0.4, size=0.08) 

#Psudo-Absence sammple example map 
#load in the data 
# load in the pseudoabsence table from the cluster then plot 
load("pseudoabsences.RData")
background <- read.csv('background_25mm.csv')
length(which(myBiomodData_Inv@PA.table$PA1==TRUE))

background <- background[complete.cases(background),]
index <- which(myBiomodData_Inv@PA.table$PA9==TRUE) # take one example of pseudoabsences 
#prepare the data 
xy <- background[index,c(2:3)]
spdf <- SpatialPoints(xy[complete.cases(xy),], proj = r@crs)
sf <- st_as_sf(spdf)
#make the map 
Pseudo  <- tm_shape(World, projection=st_crs('+proj=moll')) + #fill="grey95", col="grey75",
  tm_polygons() +
  tm_layout(title = "",
            frame=F,
            bg.color = "#f1fefcff",
            earth.boundary=T,
            earth.boundary.color = "black",
            space.color = "white",
            earth.boundary.lwd=1.5) + 
  tm_shape(sf) + 
  tm_dots(col = "deeppink4", alpha=0.4, size=0.08) 

#save figure 
figS3 <- tmap_arrange(lful, back, Pseudo)
tmap_save(filename="figS3.png", tm=figS3, width=7.5, height= 10, units = "in")
#ggsave(file="maps.pdf", plot=plot, width=7.5, height= 10, units = "in", device='pdf')




