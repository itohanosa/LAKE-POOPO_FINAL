library(raster)
library(rgdal)
library(rgeos)
x <- "F:/Itohan-Osa/Lake Poopo"
setwd(x)

LS_2001 <-brick("F:/Itohan-Osa/Lake Poopo/new_1a.tif")
LS_2001

LS_2018a <-brick("F:/Itohan-Osa/Lake Poopo/new_18a.tif")
LS_2018a

TD_01_18a <- shapefile("F:/Itohan-Osa/Lake Poopo/new_shp.shp")
TD_01_18a

LS_2017<- brick("F:/Itohan-Osa/Lake Poopo/2017CROP.tif")
LS_2017

LS_2018b <- brick("F:/Itohan-Osa/Lake Poopo/2018CROP.tif")
LS_2018b


TD_17_18b <- shapefile("F:/Itohan-Osa/shapefiles.shp")
TD_17_18b

nlayers(LS_2001)
nlayers(LS_2017)
nlayers(LS_2018a)
nlayers(LS_2018b)
nlayers(TD_01_18a)
nlayers(TD_17_18b)


crs(LS_2001)
crs(LS_2017)
crs(LS_2018a)
crs(LS_2018b)
crs(TD_01_18a)
crs(TD_17_18b)


ncell(LS_2001)
ncell(LS_2017)
ncell(LS_2018a)
ncell(LS_2018b)
ncell(TD_01_18a)
ncell(TD_17_18b)

dim(LS_2001)
dim(LS_2017)
dim(LS_2018a)
dim(LS_2018b)
dim(TD_01_18a)
dim(TD_17_18b)


res(LS_2001)
res(LS_2017)
res(LS_2018a)
res(LS_2018b)


plotRGB(LS_2001, r=5,g=4,b=3, stretch='lin')
plotRGB(LS_2017, r=5,g=4,b=3, stretch='lin')
plotRGB(LS_2018a, r=5,g=4,b=3, stretch='lin')
plotRGB(LS_2018b, r=5,g=4,b=3, stretch='lin')
plot(TD_01_18a,col="green",border="black", add=TRUE)
plot(TD_17_18b,col="blue",border="black", add=TRUE)

library(ggplot2)
library(maptools)
library(dplyr)
#to add the north arrow and a scale bar
library(ggsn)
options(stringsAsFactors = FALSE)

#Bolivia's shapefile
library(ggmap)
qmap("bolivia", maptype="hybrid", zoom = 6)

qmap("poopo", maptype="hybrid", zoom = 7)

bolivia <- getData("GADM", country= "BOL", level=1)
plot(bolivia)

#rainfall in bolivia
#prec <- getData("worldclim", var="bio12", res=2.5)
prec <- getData("worldclim", var="prec", res=2.5)
names(prec) <- c("Jan", "Feb", "Mar", "Apr", "May", "June", "July", "Aug", "Sept", "Oct", "Nov", "Dec")
prec
plot(prec)

#crop Bolivia
prec_Bol <- mask(crop(prec, bolivia),bolivia)
prec_Bol
plot(prec_Bol)
writeRaster(prec_Bol, dataType='INT1S', filename = 'prec_Bol.tif', format="GTiff", overwrite=TRUE)


#ndvi
LS_2001_ndvi <- (LS_2001[[4]]-LS_2001[[3]])/(LS_2001[[4]]+LS_2001[[3]])
LS_2001_ndvi
plot(LS_2001_ndvi)
writeRaster(LS_2001_ndvi, dataType='INT1S', filename = 'LS_2001_ndvi.tif', format="GTiff", overwrite=TRUE)

LS_2017_ndvi <- (LS_2017[[5]]-LS_2017[[4]])/(LS_2017[[5]]+LS_2017[[4]])
LS_2017_ndvi
plot(LS_2017_ndvi)
hist(LS_2017_ndvi)
writeRaster(LS_2017_ndvi, dataType='INT1S', filename = 'LS_2017_ndvi.tif', format="GTiff", overwrite=TRUE)

LS_2018a_ndvi <- (LS_2018a[[5]]-LS_2018a[[4]])/(LS_2018a[[5]]+LS_2018a[[4]])
LS_2018a_ndvi
plot(LS_2018a_ndvi)
#hist(LS_2018_ndvi)
writeRaster(LS_2018a_ndvi, dataType='INT1S', filename = 'LS_2018a_ndvi.tif', format="GTiff", overwrite=TRUE)


LS_2018b_ndvi <- (LS_2018b[[5]]-LS_2018b[[4]])/(LS_2018b[[5]]+LS_2018b[[4]])
LS_2018b_ndvi
plot(LS_2018b_ndvi)
#hist(LS_2018_ndvi)
writeRaster(LS_2018b_ndvi, dataType='INT1S', filename = 'LS_2018b_ndvi.tif', format="GTiff", overwrite=TRUE)


#NDVI DIFFERENCE MAP
LS_2001_ndvi <- raster("F:/Itohan-Osa/Lake Poopo/LS_2001_ndvi.tif")
LS_2017_ndvi <- raster("F:/Itohan-Osa/Lake Poopo/LS_2017_ndvi.tif")
LS_2018a_ndvi <- raster("F:/Itohan-Osa/Lake Poopo/LS_2018a_ndvi.tif")
LS_2018b_ndvi <- raster("F:/Itohan-Osa/Lake Poopo/LS_2018b_ndvi.tif")

LS_01_18_ndvi_diff <- ((LS_2001_ndvi)-(LS_2018a_ndvi))
LS_01_18_ndvi_diff
plot(LS_01_18_ndvi_diff)
hist(LS_01_18_ndvi_diff)
writeRaster(LS_01_18_ndvi_diff, dataType='INT1S', filename = 'LS_01_18_ndvi_diff.tif', format="GTiff", overwrite=TRUE)

LS_17_18_ndvi_diff <- LS_2017_ndvi-LS_2018b_ndvi
plot(LS_17_18_ndvi_diff)
hist(LS_17_18_ndvi_diff)
writeRaster(LS_17_18_ndvi_diff, dataType='INT1S', filename = 'LS_17_18_ndvi_diff.tif', format="GTiff", overwrite=TRUE)



#SR or RVI = NIR/RED by Jordan(1969)
#LS_2017_RVI <- ((LS_2017[[5]])/(LS_2017[[4]]))
#LS_2017_RVI

#normalized difference water index
#ndwi <- ((green-nir)/(green+nir)) defined by McFeeters(1996)
LS_2001_ndwi <- (LS_2001[[2]]-LS_2001[[4]])/(LS_2001[[2]]+LS_2001[[4]])
LS_2001_ndwi
plot(LS_2001_ndwi)
writeRaster(LS_2001_ndwi, dataType='INT1S', filename = 'LS_2001_ndwi.tif', format="GTiff", overwrite=TRUE)


LS_2017_ndwi <- (LS_2017[[3]]-LS_2017[[5]])/(LS_2017[[3]]+LS_2017[[5]])
LS_2017_ndwi
plot(LS_2017_ndwi)
hist(LS_2017_ndwi)
writeRaster(LS_2017_ndwi, dataType='INT1S', filename = 'LS_2017_ndwi.tif', format="GTiff", overwrite=TRUE)

LS_2018a_ndwi <- (LS_2018a[[3]]-LS_2018a[[5]])/(LS_2018a[[3]]+LS_2018a[[5]])
LS_2018a_ndwi
plot(LS_2018a_ndwi)
hist(LS_2018a_ndwi)
writeRaster(LS_2018a_ndwi, dataType='INT1S', filename = 'LS_2018a_ndwi.tif', format="GTiff", overwrite=TRUE)


LS_2018b_ndwi <- (LS_2018b[[3]]-LS_2018b[[5]])/(LS_2018b[[3]]+LS_2018b[[5]])
LS_2018b_ndwi
plot(LS_2018b_ndwi)
hist(LS_2018b_ndwi)
writeRaster(LS_2018b_ndwi, dataType='INT1S', filename = 'LS_2018b_ndwi.tif', format="GTiff", overwrite=TRUE)

#NDwI DIFFERENCE MAP

LS_2001_ndwi <- raster("F:/Itohan-Osa/Lake Poopo/LS_2001_ndwi.tif")
LS_2017_ndwi <- raster("F:/Itohan-Osa/Lake Poopo/LS_2017_ndwi.tif")
LS_2018a_ndwi <- raster("F:/Itohan-Osa/Lake Poopo/LS_2018a_ndwi.tif")
LS_2018b_ndwi <- raster("F:/Itohan-Osa/Lake Poopo/LS_2018b_ndwi.tif")

LS_01_18_ndwi_diff <- LS_2001_ndwi-LS_2018a_ndwi
plot(LS_01_18_ndwi_diff)
hist(LS_01_18_ndwi_diff)
writeRaster(LS_01_18_ndwi_diff, dataType='INT1S', filename = 'LS_01_18_ndwi_diff.tif', format="GTiff", overwrite=TRUE)

LS_17_18_ndwi_diff <- LS_2017_ndwi-LS_2018b_ndwi
plot(LS_17_18_ndwi_diff)
hist(LS_17_18_ndwi_diff)
writeRaster(LS_17_18_ndwi_diff, dataType='INT1S', filename = 'LS_17_18_ndwi_diff.tif', format="GTiff", overwrite=TRUE)

#spatial correlation of NDVI and NDWI
LS_2001_pca <- rasterPCA(LS_2001)
summary(LS_2001_pca$model)
loadings(LS_2001_pca$model)


#classification

library(caret)
library(randomForest)
library(e1071)
library(RStoolbox)
library(kernlab)

#rf

TD_2001 <- rgdal::readOGR("F:/Itohan-Osa/classification_2001.shp")
LS_2001_rf<- superClass(img = LS_2001, model="rf", nSamples = 1000, trainData = TD_2001, responseCol = "X2001", na.omit=TRUE, trainPartition = 0.7, verbose=TRUE)
head(LS_2001_rf)
LS_2001_rf
LS_2001_rf$classMapping

y_01 <- LS_2001_rf$map
ggR(y_01, geom_raster = TRUE)+
  theme(legend.title = element_text(size = 12, face = "bold"))+
  theme(legend.text = element_text(size = 10))+
  scale_fill_manual(values = c("grey", "blue"),
                    name="classification_2001")+
  xlab("")
ylab("")
north2(img_plot, x=0.67, y=0.9, scale = 0.1, symbol = 1)
barplot(y_01, breaks=2, col=c("blue", "green"), horiz=FALSE, digits=NULL, 
        las=1,names.arg=c("wetland", "fertile_land"))
writeRaster(y_01, dataType='FLT4S', filename = 'LS_01_class.tif', format="GTiff", overwrite=TRUE)

#nnet
LS_2001_nnet<- superClass(img = LS_2001, model="nnet", nSamples = 1000, trainData = TD_2001, responseCol = "X2001", na.omit=TRUE, trainPartition = 0.7, verbose = TRUE)
head(LS_2001_nnet)
LS_2001_nnet
LS_2001_nnet$classMapping

y_02a <- LS_2001_nnet$map
ggR(y_02a, geom_raster = TRUE)+
  theme(legend.title = element_text(size = 12, face = "bold"))+
  theme(legend.text = element_text(size = 10))+
  scale_fill_manual(values = c("grey", "blue"),
                    name="classification_2001")+
  xlab("")
ylab("")
north2(img_plot, x=0.67, y=0.9, scale = 0.1, symbol = 1)
barplot(y_02a, breaks=2, col=c("blue", "green"), horiz=FALSE, digits=NULL, 
        las=1,names.arg=c("wetland", "fertile_land"))
writeRaster(y_02a, dataType='FLT4S', filename = 'LS_01_nnet_class.tif', format="GTiff", overwrite=TRUE)
TRUE
LS_2017
TD_2017 <- rgdal::readOGR("F:/Itohan-Osa/classification_2017.shp")
TD_2017

LS_2017_rf<- superClass(img = LS_2017, model="rf", nSamples = 1000, trainData = TD_2017, na.omit=TRUE, trainPartition = 0.7, responseCol = "classifica", verbose = TRUE)
LS_2017_rf
LS_2017_rf$classMapping
y_17 <- LS_2017_rf$map
ggR(y_17, geom_raster = TRUE)+
  theme(legend.title = element_text(size = 12, face = "bold"))+
  theme(legend.text = element_text(size = 10))+
  scale_fill_manual(values = c("brown", "green", "blue"),
                    name="classification_17_rf")+
  xlab("")
ylab("")
north2(img_plot, x=0.67, y=0.9, scale = 0.1, symbol = 1)
#close(LS_2017_rf)
barplot(y_17, breaks=3, col=c("blue", "green","brown"), horiz=FALSE, digits=NULL, 
        las=1,names.arg=c("wetland", "vegetation", "barren land"))

writeRaster(y_17, dataType='FLT4S', filename = 'LS_17_class.tif', format="GTiff", overwrite=TRUE)

LS_2017_nnet<- superClass(img = LS_2017, model="nnet", nSamples = 1000, trainData = TD_2017, na.omit=TRUE, trainPartition = 0.7, responseCol = "classifica", verbose = TRUE)
LS_2017_nnet
LS_2017_nnet$classMapping
y_17nnet <- LS_2017_rf$map
ggR(y_17nnet, geom_raster = TRUE)+
  theme(legend.title = element_text(size = 12, face = "bold"))+
  theme(legend.text = element_text(size = 10))+
  scale_fill_manual(values = c("brown", "green", "blue"),
                    name="classification_17_rf")+
  xlab("")
ylab("")
north2(img_plot, x=0.67, y=0.9, scale = 0.1, symbol = 1)
#close(LS_2017_rf)
barplot(y_17nnet, breaks=3, col=c("blue", "green","brown"), horiz=FALSE, digits=NULL, 
        las=1,names.arg=c("wetland", "vegetation", "barren land"))

writeRaster(y_17nnet, dataType='FLT4S', filename = 'LS_17nnet_class.tif', format="GTiff", overwrite=TRUE)


LS_2018b
TD_2018 <- rgdal::readOGR("F:/Itohan-Osa/classification_2018.shp")
TD_2018


#rf
LS_2018_rf<- superClass(img = LS_2018b, model="rf", nSamples = 1000, trainData = TD_2018, trainPartition = 0.5, na.omit=TRUE, responseCol = "class2018", verbose = TRUE)
LS_2018_rf
LS_2018_rf$classMapping
y_18 <- LS_2018_rf$map
ggR(y_18, geom_raster = TRUE)+
  theme(legend.title = element_text(size = 12, face = "bold"))+
  theme(legend.text = element_text(size = 10))+
  scale_fill_manual(values = c("brown", "green", "blue"),
                    name="classification_18")+
  xlab("")
ylab("")
north2(img_plot, x=0.67, y=0.9, scale = 0.1, symbol = 1)

barplot(y_18, breaks=3, col=c("brown", "green","blue"), horiz=FALSE, digits=NULL, 
        las=1,names.arg=c("barren land", "vegetation", "wetland"))

#LS_2018_rf_1 <- superClass(img = LS_2018, trainData = TD_2018, valData = validationPolygons, responseCol = "class")

#LS_2018_rf_1$validation$performance

writeRaster(y_18, dataType='FLT4S', filename = 'LS_18_class.tif', format="GTiff", overwrite=TRUE)


#nnet
LS_2018_nnet<- superClass(img = LS_2018b, model="nnet", nSamples = 1000, trainData = TD_2018, trainPartition = 0.5, na.omit=TRUE, responseCol = "class2018", verbose = TRUE)
LS_2018_nnet
LS_2018_nnet$classMapping
y_18nnet <- LS_2018_nnet$map
ggR(y_18nnet, geom_raster = TRUE)+
  theme(legend.title = element_text(size = 12, face = "bold"))+
  theme(legend.text = element_text(size = 10))+
  scale_fill_manual(values = c("brown", "green", "blue"),
                    name="classification_18")+
  xlab("")
ylab("")
north2(img_plot, x=0.67, y=0.9, scale = 0.1, symbol = 1)

barplot(y_18nnet, breaks=3, col=c("brown", "green","blue"), horiz=FALSE, digits=NULL, 
        las=1,names.arg=c("barren land", "vegetation", "wetland"))

#LS_2018_rf_1 <- superClass(img = LS_2018, trainData = TD_2018, valData = validationPolygons, responseCol = "class")

#LS_2018_rf_1$validation$performance

writeRaster(y_18nnet, dataType='FLT4S', filename = 'LS_18_class.tif', format="GTiff", overwrite=TRUE)


#LS_2017 <- tasseledCap(dropLayer("F:/Itohan-Osa/Lake Poopo/2017CROP.tif", "B5_bt"), sat = "Landsat8OLI")

#CHANGE DETECTION - MULTI-DATE CLASSIFICATION
library(RStoolbox)
library(caret)
library(kernlab)
library(e1071)

#cV<- rasterCVA(LS_2001[[1:2]], LS_2018a[[1:2]])
LS_01_18 <- stack(LS_2001, LS_2018a)

TD1 <- rgdal::readOGR("F:/Itohan-Osa/Lake Poopo/new_shp.shp")

crs(TD1)
extent(TD1)
library(RStoolbox)

LS_Change1<- superClass(img = LS_01_18, model="rf", nSamples = 1000, trainData = TD1,trainPartition = 0.5, na.omit=TRUE, responseCol = "new_shp", verbose=TRUE)
LS_Change1
LS_Change1$classMapping
z_1 <- LS_Change1$map
ggR(z_1, geom_raster = TRUE)+
  theme(legend.title = element_text(size = 12, face = "bold"))+
  theme(legend.text = element_text(size = 10))+
  scale_fill_manual(values = c("brown","purple","blue"),
                    name="LS_CHANGE1_RF")+
  xlab("")
+ylab("")
writeRaster(z_1, dataType='FLT4S', filename = 'LS_01_18Change_SVM.tif', format="GTiff", overwrite=TRUE)

#SVM
LS_01_18 <- stack(LS_2001, LS_2018a)
TD1 <- rgdal::readOGR("F:/Itohan-Osa/Lake Poopo/new_shp.shp")

crs(TD1)
extent(TD1)
library(RStoolbox)
library(e1071)

LS_Change1a<- superClass(img = LS_01_18, model="svm", nSamples = 1000, trainData = TD1, na.omit=TRUE, responseCol = "new_shp", verbose = TRUE)

LS_Change1
LS_Change1$classMapping
z_svm_01_18 <- LS_Change1$map
ggR(z_svm_01_18, geom_raster = TRUE)+
  theme(legend.title = element_text(size = 12, face = "bold"))+
  theme(legend.text = element_text(size = 10))+
  scale_fill_manual(values = c("brown","purple","blue"),
                    name="LS_CHANGE_01_18_SVM")+
  xlab("")
+ylab("")
writeRaster(z_svm_01_18, dataType='FLT4S', filename = 'LS_01_18Change_SVM.tif', format="GTiff", overwrite=TRUE)

#NNET
LS_01_18 <- stack(LS_2001, LS_2018a)
TD1 <- rgdal::readOGR("F:/Itohan-Osa/Lake Poopo/new_shp.shp")

crs(TD1)
extent(TD1)
library(RStoolbox)
library(e1071)
LS_Change1<- superClass(img = LS_01_18, model="nnet", nSamples = 1000, trainData = TD1, na.omit=TRUE, responseCol = "new_shp")
LS_Change1
LS_Change1$classMapping
z_nnet_01_18 <- LS_Change1$map
ggR(z_nnet_01_18, geom_raster = TRUE)+
  theme(legend.title = element_text(size = 12, face = "bold"))+
  theme(legend.text = element_text(size = 10))+
  scale_fill_manual(values = c("brown","purple","blue"),
                    name="LS_CHANGE_01_18_NNET")+
  xlab("")
+ylab("")
writeRaster(z_nnet_01_18, dataType='FLT4S', filename = 'LS_01_18Change_NNET.tif', format="GTiff", overwrite=TRUE)
###############################################

#cV<- rasterCVA(LS_2017[[1:2]], LS_2018b[[1:2]])
LS_17_18 <- stack(LS_2017, LS_2018b)
TD <- rgdal::readOGR("F:/Itohan-Osa/shapefiles.shp")

TD

#RF
LS_Change<- superClass(img = LS_17_18, model="rf", nSamples = 1000, trainData = TD, na.omit=TRUE, responseCol = "class")
LS_Change
LS_Change$classMapping
y_rf_17_18 <- LS_Change$map
ggR(y_rf_17_18, geom_raster = TRUE)+
  theme(legend.title = element_text(size = 12, face = "bold"))+
  theme(legend.text = element_text(size = 10))+
  scale_fill_manual(values = c("brown","purple", "green", "pink", "yellow", "orange", "blue"),
                    name="change_map_RF")+
  xlab("")
+ylab("")
north2(img_plot, x=0.67, y=0.9, scale = 0.1, symbol = 1)

barplot(y_rf_17_18, breaks=10, col=c("brown","purple", "green", "pink", "yellow", "orange", "blue"), horiz=FALSE, digits=NULL, 
        las=1,names.arg=c("bar to bar","veg to bar","veg to veg","veg to wetland", "wetland to bar", "wetland to veg", "wetland to wetland"))

writeRaster(y_rf_17_18, dataType='FLT4S', filename = 'LS_17_18Change_RF.tif', format="GTiff", overwrite=TRUE)


#SVM
LS_Change<- superClass(img = LS_17_18, model="svm", nSamples = 1000, trainData = TD, na.omit=TRUE, responseCol = "class")
LS_Change
LS_Change$classMapping
y_svm_17_18 <- LS_Change$map
ggR(y_svm_17_18, geom_raster = TRUE)+
  theme(legend.title = element_text(size = 12, face = "bold"))+
  theme(legend.text = element_text(size = 10))+
  scale_fill_manual(values = c("brown","purple", "green", "pink", "yellow", "orange", "blue"),
                    name="change_map_SVM")+
  xlab("")
  +ylab("")
north2(img_plot, x=0.67, y=0.9, scale = 0.1, symbol = 1)

writeRaster(y_svm_17_18, dataType='FLT4S', filename = 'LS_Change_SVM.tif', format="GTiff", overwrite=TRUE)

#NNET
LS_Change<- superClass(img = LS_17_18, model="nnet", nSamples = 1000, trainData = TD, na.omit=TRUE,responseCol = "class")
LS_Change
LS_Change$classMapping
y_nnet_17_18 <- LS_Change$map
ggR(y_nnet_17_18, geom_raster = TRUE)+
  theme(legend.title = element_text(size = 12, face = "bold"))+
  theme(legend.text = element_text(size = 10))+
  scale_fill_manual(values = c("brown","purple", "green", "pink", "yellow", "orange", "blue"),
                    name="change_map_NNET")+
  xlab("")
+ylab("")
north2(img_plot, x=0.67, y=0.9, scale = 0.1, symbol = 1)

writeRaster(y_nnet_17_18, dataType='FLT4S', filename = 'LS_Change_NNET.tif', format="GTiff", overwrite=TRUE)


#comparison of indices for 2001
Ls_vi_01 <- spectralIndices(LS_2001, red= 3, nir = 4, green= 2,indices = c("NDVI", "DVI", "NDWI"))
LS_vi_sd <- calc(Ls_vi_01, fun=sd)
LS_vi_ndvi_ndwi_cor <- corLocal(LS_2001$NDVI, LS_2001$NDWI, ngb=11, method = "spearman")
LS_2001_pca <- rasterPCA(LS_2001)
summary(LS_2001_pca$model)
loadings(LS_2001_pca$model)

#comparison of indices for 2017
Ls_vi_17 <- spectralIndices(LS_2017, red= 4, nir = 5, green= 3,indices = c("NDVI", "DVI", "NDWI"))
LS_vi_sd <- calc(Ls_vi_17, fun=sd)
plot(LS_vi_sd)

LS_vi_ndvi_ndwi_cor <- corLocal(LS_2017$NDVI, LS_2017$NDWI, ngb=11, method = "spearman")
LS_2017_pca <- rasterPCA(LS_2017)

summary(LS_2017_pca$model)
loadings(LS_2017_pca$model)
plot(LS_2017_pca$model)

chkBands <- c("X2017CROP.1", "X2017CROP.2", "X2017CROP.3", "X2017CROP.4", "X2017CROP.5")
com_2017 <- tasseledCap(LS_2017[[chkBands]], sat = "Landsat8 OLI")

#comparison of indices for 2018
Ls_vi_18 <- spectralIndices(LS_2018, red= 4, nir = 5, green= 3,indices = c("NDVI", "DVI", "NDWI"))
LS_vi_sd <- calc(Ls_vi_18, fun=sd)
LS_vi_ndvi_ndwi_cor <- corLocal(LS_2018$NDVI, LS_2018$NDWI, ngb=11, method = "spearman")
LS_2018_pca <- rasterPCA(LS_2018)
summary(LS_2018_pca$model)
loadings(LS_2018_pca$model)
