
# figure 3 plot

# ------------------------------------------------------------------

# load packages
require(ggplot2)
require(Rmisc)
require(ggplot2)
require(grid)
require(gridExtra)
require(viridis)

# Extra functions
# ------------------------------------------------------------------

fmt_dcimals <- function(decimals=0){
  function(x) format(x,nsmall = decimals,scientific = FALSE)
}

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# ------------------------------------------------------------------

# Read in observed source data - This is the observed data of hrp2 deletions collected by Parr at al 2016 (Pfhrp2-deleted Plasmodium falciparum parasites in the Democratic Republic of Congo: A national cross-sectional survey. J Infect Dis: 1–34. doi:10.1093/ecco-jcc/jjw024),
fig3_observed <- read.csv(system.file("extdata","fig3.csv",package="hrp2malaRia"))
fig3a_simulated <- read.csv(system.file("extdata","fig3a.csv",package="hrp2malaRia"))

df <- as.data.frame(list("Prevalence"=fig3a_simulated$PCR.PfPR.6.59.months,
                         "Dels"=fig3a_simulated$Proportion.Population.only.infected.with.pfhrp2.deleted.mutants,
                         "Province"=fig3a_simulated$Province))
dfmean <- as.data.frame(list("Prevalence"=colMeans(matrix(df$Prevalence,nrow=10,byrow=T)),
                             "Dels"=colMeans(matrix(df$Dels,nrow=10,byrow=T)),
                             "Dels_U95"=apply((matrix(df$Dels,nrow=10,byrow=T)),2,Rmisc::CI)[1,],
                             "Dels_L95"=apply((matrix(df$Dels,nrow=10,byrow=T)),2,Rmisc::CI)[3,],
                             "Prevalence_U95"=apply((matrix(df$Prevalence,nrow=10,byrow=T)),2,Rmisc::CI)[1,],
                             "Prevalence_L95"=apply((matrix(df$Prevalence,nrow=10,byrow=T)),2,Rmisc::CI)[3,],
                             "Province"=rep(df$Province,1)))
##
# ------------------------------------------------------------------
fig3a <- ggplot(data = fig3_observed,aes(x=PCR.05,y=propDel,label=Province)) +
  geom_point(colour="black",alpha=.33) +
  geom_line(data=data.frame(loess.smooth(fig3_observed$PCR.05,fig3_observed$propDel,span =0.7)),aes(x=x,y=y,colour="y"),inherit.aes = F) +
  geom_errorbar(data = fig3_observed, aes(x = PCR.05, y = propDel, ymin = propDel.L95, ymax = propDel.U95),
                colour = 'grey', width = 0.4,alpha=0.66,size=0.25) +
  geom_errorbarh(data = fig3_observed, aes(x = PCR.05, y = propDel, xmin = PCR.05.L95, xmax = PCR.05.U95),
                 colour = 'grey',alpha=0.66,size=0.25) +
  geom_point(data = dfmean, aes(x=Prevalence,y = Dels,colour="Dels"),inherit.aes = F,alpha=0.8) +
  geom_smooth(data = df,aes(x=Prevalence,y = Dels),inherit.aes = F,se=F,span=0.7,method = "loess",colour="red",alpha=0.75,size=0.5) +
  xlab("PCR PfPR 6-59 months") + ylab(expression(paste("Proportion infected with only ",italic(pfhrp2),"-deleted mutants"))) +
  scale_colour_manual("",breaks = c("y", "Dels"),values = c("y"="black", "Dels"="red"),labels=c("DHS Dataset","Simulation")) +
  geom_errorbar(data = dfmean, aes(x = Prevalence, y = Dels, ymin = Dels_L95, ymax = Dels_U95),colour = 'red', width = 0.4,size=0.25) +
  geom_errorbarh(data = dfmean, aes(x = Prevalence, y = Dels, xmin = Prevalence_L95, xmax = Prevalence_U95),colour = 'red',size=0.25) +
  #geom_vline(xintercept = DRC_Dels$Updated.PCR.05,linetype="dashed",color="grey",alpha=0.5) +
  theme_classic() + theme_light() + theme_linedraw() +
  theme(axis.title.x = element_text(margin=margin(10)), axis.title.y = element_text(margin=margin(c(0,10))),
        text=element_text(size=10, family="Times New Roman"),plot.title=element_text(size=12,family="Times New Roman"),
        panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank())


# Read in observed source data - This is the observed data of hrp2 deletions collected by Parr at al 2016 (Pfhrp2-deleted Plasmodium falciparum parasites in the Democratic Republic of Congo: A national cross-sectional survey. J Infect Dis: 1–34. doi:10.1093/ecco-jcc/jjw024),
fig3_observed <- read.csv(system.file("extdata","fig3.csv",package="hrp2malaRia"))
fig3b_simulated <- read.csv(system.file("extdata","fig3b.csv",package="hrp2malaRia"))

df <- as.data.frame(list("Prevalence"=fig3b_simulated$PCR.PfPR.6.59.months,
                         "Dels"=fig3b_simulated$Proportion.Population.only.infected.with.pfhrp2.deleted.mutants,
                         "Province"=fig3b_simulated$Province))
dfmean <- as.data.frame(list("Prevalence"=colMeans(matrix(df$Prevalence,nrow=10,byrow=T)),
                             "Dels"=colMeans(matrix(df$Dels,nrow=10,byrow=T)),
                             "Dels_U95"=apply((matrix(df$Dels,nrow=10,byrow=T)),2,Rmisc::CI)[1,],
                             "Dels_L95"=apply((matrix(df$Dels,nrow=10,byrow=T)),2,Rmisc::CI)[3,],
                             "Prevalence_U95"=apply((matrix(df$Prevalence,nrow=10,byrow=T)),2,Rmisc::CI)[1,],
                             "Prevalence_L95"=apply((matrix(df$Prevalence,nrow=10,byrow=T)),2,Rmisc::CI)[3,],
                             "Province"=rep(df$Province,1)))
##
# ------------------------------------------------------------------
fig3b <- ggplot(data = fig3_observed,aes(x=PCR.05,y=propDel,label=Province)) +
  geom_point(colour="black",alpha=.33) +
  geom_line(data=data.frame(loess.smooth(fig3_observed$PCR.05,fig3_observed$propDel,span =0.7)),aes(x=x,y=y,colour="y"),inherit.aes = F) +
  geom_errorbar(data = fig3_observed, aes(x = PCR.05, y = propDel, ymin = propDel.L95, ymax = propDel.U95),
                colour = 'grey', width = 0.4,alpha=0.66,size=0.25) +
  geom_errorbarh(data = fig3_observed, aes(x = PCR.05, y = propDel, xmin = PCR.05.L95, xmax = PCR.05.U95),
                 colour = 'grey',alpha=0.66,size=0.25) +
  geom_point(data = dfmean, aes(x=Prevalence,y = Dels,colour="Dels"),inherit.aes = F,alpha=0.8) +
  geom_smooth(data = df,aes(x=Prevalence,y = Dels),inherit.aes = F,se=F,span=0.7,method = "loess",colour="red",alpha=0.75,size=0.5) +
  xlab("PCR PfPR 6-59 months") + ylab(expression(paste("Proportion infected with only ",italic(pfhrp2),"-deleted mutants"))) +
  scale_colour_manual("",breaks = c("y", "Dels"),values = c("y"="black", "Dels"="red"),labels=c("DHS Dataset","Simulation")) +
  geom_errorbar(data = dfmean, aes(x = Prevalence, y = Dels, ymin = Dels_L95, ymax = Dels_U95),colour = 'red', width = 0.4,size=0.25) +
  geom_errorbarh(data = dfmean, aes(x = Prevalence, y = Dels, xmin = Prevalence_L95, xmax = Prevalence_U95),colour = 'red',size=0.25) +
  #geom_vline(xintercept = DRC_Dels$Updated.PCR.05,linetype="dashed",color="grey",alpha=0.5) +
  theme_classic() + theme_light() + theme_linedraw() +
  theme(axis.title.x = element_text(margin=margin(10)), axis.title.y = element_text(margin=margin(c(0,10))),
        text=element_text(size=10, family="Times New Roman"),plot.title=element_text(size=12,family="Times New Roman"),
        panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank())

## PLOTTING
# ------------------------------------------------------------------

g1 <- ggplotGrob(fig3a)
g2 <- ggplotGrob(fig3b)
g2$widths <- g1$widths
myplot1 <- arrangeGrob(g1, top = textGrob("a", x = unit(0, "npc")
                                          , y   = unit(1, "npc"), just=c("left","top"),
                                          gp=gpar(col="black", fontsize=18, fontfamily="Times New Roman")))

myplot2 <- arrangeGrob(g2, top = textGrob("b", x = unit(0, "npc")
                                          , y = unit(1, "npc"), just=c("left","top"),
                                          gp=gpar(col="black", fontsize=18, fontfamily="Times New Roman")))
windows()
gridExtra::grid.arrange(myplot1,myplot2,ncol=1)



# figure 3 - supplement figure 1 plot

# ------------------------------------------------------------------

# Read in observed source data - This is the observed data of hrp2 deletions collected by Parr at al 2016 (Pfhrp2-deleted Plasmodium falciparum parasites in the Democratic Republic of Congo: A national cross-sectional survey. J Infect Dis: 1–34. doi:10.1093/ecco-jcc/jjw024),
fig3_observed <- read.csv(system.file("extdata","fig3.csv",package="hrp2malaRia"))
figs3_s1_simulated <- read.csv(system.file("extdata","fig3_s1.csv",package="hrp2malaRia"))

df <- as.data.frame(list("Prevalence"=figs3_s1_simulated$PCR.PfPR.6.59.months,
                         "Dels"=figs3_s1_simulated$Proportion.Population.only.infected.with.pfhrp2.deleted.mutants,
                         "Province"=figs3_s1_simulated$Province))
dfmean <- as.data.frame(list("Prevalence"=colMeans(matrix(df$Prevalence,nrow=10,byrow=T)),
                             "Dels"=colMeans(matrix(df$Dels,nrow=10,byrow=T)),
                             "Dels_U95"=apply((matrix(df$Dels,nrow=10,byrow=T)),2,Rmisc::CI)[1,],
                             "Dels_L95"=apply((matrix(df$Dels,nrow=10,byrow=T)),2,Rmisc::CI)[3,],
                             "Prevalence_U95"=apply((matrix(df$Prevalence,nrow=10,byrow=T)),2,Rmisc::CI)[1,],
                             "Prevalence_L95"=apply((matrix(df$Prevalence,nrow=10,byrow=T)),2,Rmisc::CI)[3,],
                             "Province"=rep(df$Province,1)))

dfcomps <- as.data.frame(list("DHS" = fig3_observed$propDel,"Simulation" = dfmean$Dels,"Province"=fig3_observed$Province,
                              "Diff" = abs(fig3_observed$propDel-dfmean$Dels)))

fig3_s1a <- ggplot(data = fig3_observed,aes(x=PCR.05,y=propDel,label=Province)) +
  geom_point(colour="black",alpha=.33) +
  geom_line(data=data.frame(loess.smooth(fig3_observed$PCR.05,fig3_observed$propDel,span =0.6)),aes(x=x,y=y,colour="y"),inherit.aes = F) +
  geom_errorbar(data = fig3_observed, aes(x = PCR.05, y = propDel, ymin = propDel.L95, ymax = propDel.U95),
                colour = 'grey', width = 0.4,alpha=0.66,size=0.25) +
  geom_errorbarh(data = fig3_observed, aes(x = PCR.05, y = propDel, xmin = PCR.05.L95, xmax = PCR.05.U95),
                 colour = 'grey',alpha=0.66,size=0.25) +
  geom_point(data = dfmean, aes(x=Prevalence,y = Dels,colour="Dels"),inherit.aes = F,alpha=0.8) +
  geom_smooth(data = df,aes(x=Prevalence,y = Dels),inherit.aes = F,se=F,span=0.6,method = "loess",colour="red",alpha=0.75,size=0.5) +
  xlab("PCR PfPR 6-59 months") + ylab(expression(paste("Proportion infected with only ",italic(pfhrp2),"-deleted mutants"))) +
  scale_colour_manual("",breaks = c("y", "Dels"),values = c("y"="black", "Dels"="red"),labels=c("DHS Dataset","Simulation")) +
  geom_errorbar(data = dfmean, aes(x = Prevalence, y = Dels, ymin = Dels_L95, ymax = Dels_U95),colour = 'red', width = 0.4,size=0.25) +
  geom_errorbarh(data = dfmean, aes(x = Prevalence, y = Dels, xmin = Prevalence_L95, xmax = Prevalence_U95),colour = 'red',size=0.25) +
  #geom_vline(xintercept = DRC_Dels$Updated.PCR.05,linetype="dashed",color="grey",alpha=0.5) +
  theme_classic() + theme_light() + theme_linedraw() +
  theme(axis.title.x = element_text(margin=margin(10)), axis.title.y = element_text(margin=margin(c(0,10))),
        text=element_text(size=10, family="Times New Roman"),plot.title=element_text(size=12,family="Times New Roman"),
        panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank())

fig3_s1b <- ggplot(dfcomps, aes(x=DHS,y=Simulation))+
  #ggtitle("Simulation vs DHS") +
  theme_classic() + theme_light() + theme_linedraw() +
  geom_text(data = dfcomps[dfcomps$Diff>5,][c(2),],aes(x=DHS,y=Simulation,label=Province),hjust = 0.7, nudge_y = -1,nudge_x = -0.3) +
  geom_text(data = dfcomps[dfcomps$Diff>5,][c(3),],aes(x=DHS,y=Simulation,label=Province),hjust = 0.7, nudge_y = -1,nudge_x = -0.3) +
  geom_text(data = dfcomps[dfcomps$Diff>5,][c(1),],aes(x=DHS,y=Simulation,label=Province),hjust = 0.5, nudge_y = -1) +
  geom_text(data = dfcomps[dfcomps$Diff>5,][c(4),],aes(x=DHS,y=Simulation,label=Province),hjust = 1, nudge_y = -0,nudge_x = -0.5) +
  geom_text(data = dfcomps[dfcomps$Diff>5,][c(5),],aes(x=DHS,y=Simulation,label=Province),hjust = -0.2, nudge_y = -0) +
  geom_text(data = dfcomps[dfcomps$Diff>5,][c(6),],aes(x=DHS,y=Simulation,label=Province),hjust = -0.2, nudge_y = -0) +
  geom_text(data = dfcomps[dfcomps$Diff>5,][c(7),],aes(x=DHS,y=Simulation,label=Province),hjust = -0.1, nudge_y = -0.4) +
  geom_segment(aes(x = DHS, y = Simulation,xend = DHS, yend = DHS),colour="grey") +
  geom_abline(slope=1) +
  geom_point(colour="red") +
  theme(axis.title.x = element_text(margin=margin(10)), axis.title.y = element_text(margin=margin(c(0,10))),
        text=element_text(size=10, family="Times New Roman"),plot.title=element_text(size=12,family="Times New Roman"),
        panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank())

## PLOTTING
# ------------------------------------------------------------------

g1 <- ggplotGrob(fig3_s1a)
g2 <- ggplotGrob(fig3_s1b)
g2$widths <- g1$widths
myplot1 <- arrangeGrob(g1, top = textGrob("a", x = unit(0, "npc")
                                          , y   = unit(1, "npc"), just=c("left","top"),
                                          gp=gpar(col="black", fontsize=18, fontfamily="Times New Roman")))

myplot2 <- arrangeGrob(g2, top = textGrob("b", x = unit(0, "npc")
                                          , y = unit(1, "npc"), just=c("left","top"),
                                          gp=gpar(col="black", fontsize=18, fontfamily="Times New Roman")))
windows()
gridExtra::grid.arrange(myplot1,myplot2,ncol=1)

