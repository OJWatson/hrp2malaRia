
# figure 4 plot

# ------------------------------------------------------------------

# load packages
require(maptools)
require(ggplot2)
require(maptools)
require(rgeos)
require(Cairo)
require(ggmap)
require(scales)
require(RColorBrewer)
require(viridis)
require(plotrix)

# Extra functions
# ------------------------------------------------------------------

merge.SpatialPolygonsDataFrame <- function(shp,df) {
  #### reads in a shape file and a data frame. Merges data frame with the data already present in the shapefile while preserving the order of objects (ordinary merge operation causes polygons to become disassociated with data).

  # extract data from shapefile and add key
  shp_data <- shp@data
  shp_data$mergeKey <- 1:nrow(shp_data)

  # merge with df and sort based on key
  m <- merge(shp_data,df,all.x=T)
  m <- m[order(m$mergeKey),]
  m <- subset(m,select=-mergeKey)

  # fix row names
  row.names(m) <- row.names(shp)

  # make final SpatialPolygonsDataFrame object
  s <- SpatialPolygonsDataFrame(geometry(shp), m)
  return(s)
}

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# ------------------------------------------------------------------

# Load shape files
lakes <- readRDS(system.file("extdata/shape_files","lakes.rds",package="hrp2malaRia"))
regions <- readRDS(system.file("extdata/shape_files","regions.rds",package="hrp2malaRia"))
islands <- readRDS(system.file("extdata/shape_files","islands.rds",package="hrp2malaRia"))
countries <- readRDS(system.file("extdata/shape_files","countries.rds",package="hrp2malaRia"))

# Read in source data
fig4_data <- read.csv(system.file("extdata","fig4.csv",package="hrp2malaRia"))
# Bind source data of interest to shape file
regions@data[,c("PfPR_2010","fT","hrp2")] <- fig4_data[,6:8]

# Plot window
windows()


layout(matrix(c(1,1,1,1,4,4,0,2,2,2,2,5,5,0,0,3,3,3,3,3,3,3,6,6,6,0), 2, 13, byrow = TRUE),height=c(10,20))

pal <- colorRampPalette(c("blue", "cyan", "yellow", "red"), bias=1)
nocolor=rgb(0,0,0,alpha=0.001)

## PLOTTING
# ------------------------------------------------------------------

# PREVALENCE -------------------------------
p <- regions$PfPR_2010
colourLevel <- as.numeric(cut(p,breaks=seq(0,1,l=101)))
colVec <- pal(100)[colourLevel]

# plot shape file
par(mai=c(0,0,0,0))
plot(regions,col=colVec,border=F,cex=1.2)
plot(mad,lwd=0.1,add=TRUE)
plot(lakes,border=F,col="white",add=TRUE)
plot(countries,lwd=0.1,add=TRUE)
text(x = corner.label()$x+4,y=corner.label()$y-4,"a",cex=2)

# fT -------------------------------
p <- regions$ft
colourLevel <- as.numeric(cut(p,breaks=seq(0,1,l=101)))
colVec <- pal(100)[colourLevel]

# plot shape file
#par(oma=c(0,0,0,0))
plot(regions,col=colVec,border=F)
plot(mad,lwd=1,add=TRUE)
plot(lakes,border=F,col="white",add=TRUE)
plot(countries,lwd=0.1,add=TRUE)
text(x = corner.label()$x+4,y=corner.label()$y-4,"b",cex=2)

# hrp2 -------------------------------
p <- regions$hrp2
colourLevel <- as.numeric(cut(p,breaks=seq(0,max(p,na.rm=TRUE),l=5)))
colVec <- pal(4)[colourLevel]

# plot shape file
#par(oma=c(0,0,0,0))
plot(regions,col=colVec,border=F)
plot(mad,lwd=0.1,add=TRUE)
plot(lakes,border=F,col="white",add=TRUE)
plot(countries,lwd=0.1,add=TRUE)
text(x = corner.label()$x+4,y=corner.label()$y-4,"c",cex=2)

par(mar = c(2, 2, 2, 0) + 0.1)
legend_image <- rev(as.raster(pal(100), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Microscopy \nPfPR 2-10',cex.main=0.9)
text(x=0.8, y = seq(0,1,l=5), labels = paste(seq(0,1,l=5)*100,"%",sep=""),adj = c(0, NA))
rasterImage(legend_image, 0, 0, 0.5,1)

legend_image <- rev(as.raster(matrix(pal(100), ncol=1)))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = expression(bold("f"[T])),cex.main=0.9)
text(x=0.8, y = seq(0,1,l=5), labels = paste(seq(0,1,l=5)*100,"%",sep=""),adj = c(0, NA))
rasterImage(legend_image, 0, 0, 0.5,1)

legend_image <- rev(as.raster(matrix(pal(4), ncol=1)))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'HRP2 Concern',cex.main=0.9)
text(x=1, y = seq(0.125,0.875,l=4), labels = c("Marginal","Slight","Moderate","High"),adj = c(0, NA))
rasterImage(legend_image, 0, 0, 0.5,1,interpolate = F)


# ------------------------------------------------------------------
# Figure 4 supplement 2 plot

windows()

fig4_supp2a_data <- read.csv(system.file("extdata","fig4_s2.csv",package="hrp2malaRia"),header = T,row.names = 1)
fig4_supp2b_data <- fig4_supp2a_data[,1003:2002]
fig4_supp2a_data <- fig4_supp2a_data[,1:1000]
# fT values that correspond to source data matrix row names
ft <- as.numeric(unlist(lapply(strsplit(rownames(fig4_supp2a_data),"  ",fixed=T),function(x){return(x[2])})))

# microscopy prevalences that correspond to source data col names
prev <- as.numeric(unlist(lapply(strsplit(colnames(fig4_supp2a_data)[1:1000],"....",fixed=T),function(x){return(x[2])})))

fig4_supp2a <- list()
fig4_supp2a$x <- ft
fig4_supp2a$y <- prev
fig4_supp2a$z <- as.matrix(fig4_supp2a_data)


fig4_supp2b <- list()
fig4_supp2b$x <- ft
fig4_supp2b$y <- prev
fig4_supp2b$z <- as.matrix(fig4_supp2b_data)

layout(matrix(1:4,ncol=2), width = c(6,3),height = c(10,10))
image(fig4_supp2a,col=viridis(80),xlab = expression("f"[T]),ylab="Microscopy PfPR 2-10")
u <- par("usr")
rect(u[1], u[3], u[2], u[4], col = "grey", border = "black")
image(fig4_supp2a,col=viridis(100),xlab = expression("f"[T]),ylab="Microscopy PfPR 2-10",add = TRUE)
corner.label("a",cex=2)

### Epsilon = 0

image(fig4_supp2b,col=viridis(80),xlab = expression("f"[T]),ylab="Microscopy PfPR 2-10")
u <- par("usr")
rect(u[1], u[3], u[2], u[4], col = "grey", border = "black")
image(fig4_supp2b,col=viridis(100),xlab = expression("f"[T]),ylab="Microscopy PfPR 2-10",add = TRUE)
corner.label("b",cex=2)


legend_image <- rev(as.raster(c(rep("grey",10),rep("white",10),viridis(80))))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = "", ylab="",main = "Years")
text(x=1.75, y = c(0.05,0.20 ,0.4,0.6 ,0.8,1), labels = c("20+","0","5","10","15","20"))
rasterImage(legend_image, 0, 0, 1,1)

legend_image <- rev(as.raster(c(rep("grey",10),rep("white",10),viridis(80))))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = "", ylab="",main = "Years")
text(x=1.75, y = c(0.05,0.20 ,0.4,0.6 ,0.8,1), labels = c("20+","0","5","10","15","20"))
rasterImage(legend_image, 0, 0, 1,1)


# ------------------------------------------------------------------
# Figure 4 supplement 3 plot

# Load shape files
lakes <- readRDS(system.file("extdata/shape_files","lakes.rds",package="hrp2malaRia"))
regions <- readRDS(system.file("extdata/shape_files","regions.rds",package="hrp2malaRia"))
islands <- readRDS(system.file("extdata/shape_files","islands.rds",package="hrp2malaRia"))
countries <- readRDS(system.file("extdata/shape_files","countries.rds",package="hrp2malaRia"))


# Read in source data
fig4s3_data <- read.csv(system.file("extdata","fig4_s3.csv",package="hrp2malaRia"))
# Bind source data of interest to shape file
regions@data[,c("PfPR_2010","fT","hrp2")] <- fig4s3_data[,6:8]

# Plot window
windows()


layout(matrix(c(1,1,1,1,4,4,0,2,2,2,2,5,5,0,0,3,3,3,3,3,3,3,6,6,6,0), 2, 13, byrow = TRUE),height=c(10,20))

pal <- colorRampPalette(c("blue", "cyan", "yellow", "red"), bias=1)
nocolor=rgb(0,0,0,alpha=0.001)

## PLOTTING
# ------------------------------------------------------------------

# PREVALENCE -------------------------------
p <- regions$PfPR_2010
colourLevel <- as.numeric(cut(p,breaks=seq(0,1,l=101)))
colVec <- pal(100)[colourLevel]

# plot shape file
par(mai=c(0,0,0,0))
plot(regions,col=colVec,border=F,cex=1.2)
plot(mad,lwd=0.1,add=TRUE)
plot(lakes,border=F,col="white",add=TRUE)
plot(countries,lwd=0.1,add=TRUE)
text(x = corner.label()$x+4,y=corner.label()$y-4,"a",cex=2)

# fT -------------------------------
p <- regions$ft
colourLevel <- as.numeric(cut(p,breaks=seq(0,1,l=101)))
colVec <- pal(100)[colourLevel]

# plot shape file
#par(oma=c(0,0,0,0))
plot(regions,col=colVec,border=F)
plot(mad,lwd=1,add=TRUE)
plot(lakes,border=F,col="white",add=TRUE)
plot(countries,lwd=0.1,add=TRUE)
text(x = corner.label()$x+4,y=corner.label()$y-4,"b",cex=2)

# hrp2 -------------------------------
p <- regions$hrp2
colourLevel <- as.numeric(cut(p,breaks=seq(0,max(p,na.rm=TRUE),l=5)))
colVec <- pal(4)[colourLevel]

# plot shape file
#par(oma=c(0,0,0,0))
plot(regions,col=colVec,border=F)
plot(mad,lwd=0.1,add=TRUE)
plot(lakes,border=F,col="white",add=TRUE)
plot(countries,lwd=0.1,add=TRUE)
text(x = corner.label()$x+4,y=corner.label()$y-4,"c",cex=2)

par(mar = c(2, 2, 2, 0) + 0.1)
legend_image <- rev(as.raster(pal(100), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Microscopy \nPfPR 2-10',cex.main=0.9)
text(x=0.8, y = seq(0,1,l=5), labels = paste(seq(0,1,l=5)*100,"%",sep=""),adj = c(0, NA))
rasterImage(legend_image, 0, 0, 0.5,1)

legend_image <- rev(as.raster(matrix(pal(100), ncol=1)))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = expression(bold("f"[T])),cex.main=0.9)
text(x=0.8, y = seq(0,1,l=5), labels = paste(seq(0,1,l=5)*100,"%",sep=""),adj = c(0, NA))
rasterImage(legend_image, 0, 0, 0.5,1)

legend_image <- rev(as.raster(matrix(pal(4), ncol=1)))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'HRP2 Concern',cex.main=0.9)
text(x=1, y = seq(0.125,0.875,l=4), labels = c("Marginal","Slight","Moderate","High"),adj = c(0, NA))
rasterImage(legend_image, 0, 0, 0.5,1,interpolate = F)


