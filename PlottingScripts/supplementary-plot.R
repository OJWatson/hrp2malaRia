
# supplementary file 1 plot

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
supp1_data <- read.csv(system.file("extdata","Supplementary_File_1.csv",package="hrp2malaRia"))
supp1_data$X <- NULL

# merge source data with shape file
regions@data$id <- rownames(regions@data)
mergedf <- merge(regions, supp1_data, by="Country")

# merge years
yeardf <- data.frame("Year"=mergedf$Year,"id"=regions$id)
yeardf$Year[yeardf$Year==levels(yeardf$Year)[1]] <- NA

# collate for ggplot map
fort <- fortify(regions, region = "id")
mergedf <- merge(fort, yeardf, by="id")

# PLOTTING
# ------------------------------------------------------------------

windows()


ggplot(data = mergedf, aes(x=long, y=lat, group = group,
                           fill = Year)) +
  #ggtitle(paste(formatC(months[m],digits=1,flag="0",format = "fg"),".",years[m],sep="")) +
  #ggtitle(years[m]) +
  geom_polygon(data = countries, aes(x=long, y=lat, group = group),inherit.aes=F,color = "black",fill="grey80") +
  geom_polygon(data = islands, aes(x=long, y=lat, group = group),inherit.aes=F,color = "black",fill="grey80") +
  geom_polygon(data = mergedf, aes(x=long, y=lat, group = group,fill = Year))  +
  scale_fill_viridis(name = "Year RDT Use \nReported at \nCommunity Level",discrete = T) +
  geom_polygon(data = countries, aes(x=long, y=lat, group = group),inherit.aes=F,color = "black",fill=NA) +
  geom_polygon(data = lakes, aes(x=long, y=lat, group = group),fill = "white") +
  coord_equal() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.6,family = "Times New Roman"),
        legend.title = element_text(family = "Times New Roman",size = 16),
        legend.text = element_text(family = "Times New Roman",size = 16))
