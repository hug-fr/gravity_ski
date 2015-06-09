# TODO: clean script
# Focus watershed algorithm on focal regions


library(raster)
library(rgdal)
library(EBImage)
setwd("C:\\r_wd")
elev <- raster("MNT_7laux.tif")
plot(elev)
dim(elev)

lift_xy_tmp <- read.table("rm_points.csv", header = TRUE, sep = ";")

#create empty data for resort loop

lift_xy <- matrix(nrow=nrow(lift_xy_tmp), ncol=7)
colnames(lift_xy) <- c("gid", "upx", "upy", "upalt", "lowx", "lowy", "lowalt")

fromup <- elev
fromup[] <- NA
fromlow <- elev
fromlow[] <- NA


#first loop on ski-lift to retrieve bottom and  top alt
for (i in 1:nrow(lift_xy_tmp)) {
  xy1_alt <-  extract(elev, lift_xy_tmp[i, c("x1", "y1"), drop = FALSE])
  xy2_alt <-  extract(elev, lift_xy_tmp[i, c("x2", "y2"), drop = FALSE])
  
  if (xy1_alt > xy2_alt) {
    lift_xy[i,] <- c(lift_xy_tmp[i,1],
                    lift_xy_tmp[i,"x1"], lift_xy_tmp[i,"y1"], xy1_alt,
                    lift_xy_tmp[i,"x2"], lift_xy_tmp[i,"y2"], xy2_alt
      )
  }
  else
  {lift_xy[i,] <- c(lift_xy_tmp[i,1],
                    lift_xy_tmp[i,"x2"], lift_xy_tmp[i,"y2"], xy2_alt,
                    lift_xy_tmp[i,"x1"], lift_xy_tmp[i,"y1"], xy1_alt
      )
  }
}

#create wtershed from top

ws <- watershed(as.matrix(elev))
ws_raster <- elev
ws_raster[] <- ws

for (i in 1:nrow(lift_xy)) {
  #create slope down from top of ski-lift
  
  upper_xy <- lift_xy[i, c("upx", "upy"), drop = FALSE]
  
  upper_ws <- extract(ws_raster, upper_xy)
  
  piste <- elev
  piste[ws_raster != upper_ws] <- NA
  piste[elev < min(lift_xy[i,"lowalt"])] <- NA
  
  fromup <- merge(fromup, piste)
}

plot(fromup)


#create watershhed from bottom
wsm <- watershed(max(as.matrix(elev)) - as.matrix(elev))
wsm_raster <- elev
wsm_raster[] <- wsm


for (i in 1:nrow(lift_xy)) {
  lower_xy <- lift_xy[i, c("lowx", "lowy"), drop = FALSE]
  
  lower_ws <- extract(wsm_raster, lower_xy)
  
  piste <- elev
  piste[wsm_raster != lower_ws] <- NA
  piste[elev > min(lift_xy[i,"upalt"])] <- NA
  
  fromlow <- merge(fromlow, piste)
}

plot (fromlow)

#create ski area

ski_area <- fromup
ski_area[is.na(fromlow)==TRUE] <- NA

plot(ski_area)

for (i in 1:nrow(lift_xy)){
  points(lift_xy[i,"lowx"], lift_xy[i,"lowy"], col = "red")
  points(lift_xy[i, "upx"], lift_xy[i, "upy"], col = "blue")
}

if (FALSE){


plot(ws_raster)
plot(wsm_raster)

lower_ws <- extract(wsm_raster, lower_xy)


lower_ws
upper_ws

piste <- elev
piste[ws_raster != upper_ws | wsm_raster != lower_ws] <- NA
piste[elev > upper_elev] <- NA
piste[elev < lower_elev] <- NA
plot(piste)
points(lower_xy, col = "red")
points(upper_xy, col = "blue")

slope <- terrain(elev, opt='slope')
aspect <- terrain(elev, opt='aspect')
hill <- hillShade(slope, aspect, 40, 270)
plot(hill, col=grey(0:100/100), legend=FALSE)
plot(piste, add = TRUE)
points(lower_xy, col = "red")
points(upper_xy, col = "blue")

library("topmodel")

lowcell<-rowColFromCell(elev,cellFromXY(elev, lift_xy[1,c("lowx","lowy")]))

lowcatch<- subcatch(as.matrix(elev), lowcell)
piste<-elev
piste[lowcatch == 0]<-NA

plot(piste)

}

