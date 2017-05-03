
# figure 2 plot

# ------------------------------------------------------------------
# load packages
require(ggplot2)
require(grid)
require(gridExtra)
require(viridis)

## EXTRA FUNCTIONS

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


heatmap_plot <- function(data){

  df = data


    gp <- ggplot(df,aes(x=ft,y=PfPR_grouped,fill=Years.until.20..infected.with.only.hrp2.deletions))
    gp <- gp +
      scale_fill_viridis(direction=-1,name = "Years") +
      geom_tile(aes(fill = Years.until.20..infected.with.only.hrp2.deletions), color='white')


  gp <- gp + theme_bw() +
    xlab(expression("f"[T])) +
    ylab('PfPR') +
    coord_equal(ratio=0.1) +
    theme_classic() +
    theme_light() +
    theme(axis.line=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_line(color='#eeeeee'),
          axis.title.x = element_text(margin=margin(10)),
          axis.title.y = element_text(margin=margin(c(0,10))),
          text=element_text(size=10))

  gp <- gp + scale_fill_viridis(direction=-1,name = "Years",breaks=c(0,5,10,15,20),labels=c(0,5,10,15,20),limits=c(0,20)) +
    geom_point(data = df, aes(size="20+",shape=NA), colour = "grey50") +
    guides(size=guide_legend("", override.aes=list(shape=15, size = 5,legend.key = "grey50",label.hjust=0.8))) +
    theme(legend.key =  element_rect(fill = "grey50"))

  return(gp)

}

heatmap_allele_plot <- function(data){

  df = data

  gp <- ggplot(df,aes(x=ft,y=PfPR_grouped,fill=hrp2.deletion.frequency.after.20.years))
  gp <- gp + #ggtitle(paste("Frequency of pfhrp2 gene deletion \n After 20 Years. \u03B5 = ",epsilon)) +
    scale_fill_viridis(direction=1,name = "Frequency") +
    geom_tile(aes(fill = hrp2.deletion.frequency.after.20.years), color='white')


gp <- gp + theme_bw() +
  xlab(expression("f"[T])) +
  ylab('PfPR') +
  coord_equal(ratio=0.1) +
  theme_classic() +
  theme_light() +
  theme(axis.line=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_line(color='#eeeeee'),
        axis.title.x = element_text(margin=margin(10)),
        axis.title.y = element_text(margin=margin(c(0,10))),
        text=element_text(size=10, family="Times"),
        plot.title=element_text(size=10,family="Times"))


return(gp)

}

########

fig2a <- read.csv(system.file("extdata","fig2a.csv",package="hrp2malaRia"))
fig2a_heatmap <- heatmap_plot(fig2a)
leg <- get_legend(fig2a_heatmap)
########

fig2b <- read.csv(system.file("extdata","fig2b.csv",package="hrp2malaRia"))
fig2b_heatmap <- heatmap_plot(fig2b)

########

## PLOTTING
# ------------------------------------------------------------------

fig2a_heatmap <- fig2a_heatmap + theme(legend.position="none")  +
  geom_point(aes(size=ifelse(lessthan5, "dot", "no_dot"))) +
  scale_size_manual(values=c(dot=1, no_dot=NA,"20+"=NA), guide="none")

fig2b_heatmap <- fig2b_heatmap + theme(legend.position="none") +
  geom_point(aes(size=ifelse(lessthan5, "dot", "no_dot"))) +
  scale_size_manual(values=c(dot=1, no_dot=NA,"20+"=NA), guide="none")

myplot1 <- arrangeGrob(fig2a_heatmap, top = textGrob("a", x = unit(0, "npc")
                                                , y   = unit(1, "npc"), just=c("left","top"),
                                                gp=gpar(col="black", fontsize=18, fontfamily="Times New Roman")))

myplot2 <- arrangeGrob(fig2b_heatmap, top = textGrob("b", x = unit(0, "npc")
                                                   , y = unit(1, "npc"), just=c("left","top"),
                                                   gp=gpar(col="black", fontsize=18, fontfamily="Times New Roman")))

windows()
grid.arrange(myplot1,myplot2,  leg, ncol = 3,widths = c(2,2,1))


# figure 2 - figure supplement 1 plot
# ------------------------------------------------------------------

########

fig2s1a <- read.csv(system.file("extdata","fig2_s1a.csv",package="hrp2malaRia"))
fig2s1a_heatmap <- heatmap_allele_plot(fig2s1a)
leg <- get_legend(fig2s1a_heatmap)
########

fig2s1b <- read.csv(system.file("extdata","fig2_s1b.csv",package="hrp2malaRia"))
fig2s1b_heatmap <- heatmap_allele_plot(fig2s1b)

########

## PLOTTING
# ------------------------------------------------------------------

fig2s1a_heatmap <- fig2s1a_heatmap + theme(legend.position="none")

fig2s1b_heatmap <- fig2s1b_heatmap + theme(legend.position="none")

myplot1 <- arrangeGrob(fig2s1a_heatmap, top = textGrob("a", x = unit(0, "npc")
                                                     , y   = unit(1, "npc"), just=c("left","top"),
                                                     gp=gpar(col="black", fontsize=18, fontfamily="Times New Roman")))

myplot2 <- arrangeGrob(fig2s1b_heatmap, top = textGrob("b", x = unit(0, "npc")
                                                     , y = unit(1, "npc"), just=c("left","top"),
                                                     gp=gpar(col="black", fontsize=18, fontfamily="Times New Roman")))

windows()
grid.arrange(myplot1,myplot2,  leg, ncol = 3,widths = c(2,2,1))

