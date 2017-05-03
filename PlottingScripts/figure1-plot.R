
# figure 1 plot

# ------------------------------------------------------------------

# load packages
require(ggplot2)
require(extrafont)
require(gridExtra)
require(grid)

# Extra functions
# ------------------------------------------------------------------

fmt_dcimals <- function(decimals=0){
  function(x) format(x,nsmall = decimals,scientific = FALSE)
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
# ------------------------------------------------------------------

# Read in simulation data
fig1a <- read.csv(system.file("extdata","fig1a.csv",package="hrp2malaRia"))
fig1a$Prevalence <- as.factor(fig1a$Prevalence)
fig1b <- read.csv(system.file("extdata","fig1b.csv",package="hrp2malaRia"))
fig1b$Starting.Frequency <- as.factor(fig1b$Starting.Frequency)
fig1cd <- read.csv(system.file("extdata","fig1cd.csv",package="hrp2malaRia"))

# fig1a
panels1 <- ggplot(fig1a,aes(x = Years,y=hrp2.Deletion.Allele.Frequency,color=Prevalence)) + geom_smooth(span=0.2,se=F) + ylim(0,1) +
  geom_vline(xintercept = 10) + geom_vline(xintercept = 15,linetype="dashed",color="black")
panels1 <- panels1 + theme_classic() + theme_light() + theme_linedraw() + xlab("Time (years)") + ylab(expression(paste("Proportion of strains that are ",italic(pfhrp2),"-deleted")))
panels1 <- panels1 + theme(axis.title.x = element_text(margin=margin(10)), axis.title.y = element_text(margin=margin(c(0,10))),
                           legend.justification = 'left',legend.position=c(0,0.85),legend.background = element_rect(colour = 'black')) +
  scale_color_manual(labels = c("10%", "25%","60%"),values = gg_color_hue(3)) + annotate("text", x = 7.5, y =0.56125, label = "RDT Introduction")

# fig1b
panels2 <- ggplot(fig1b,aes(x = Years,y=hrp2.Deletion.Allele.Frequency,color=Starting.Frequency)) + geom_smooth(span=0.2,se=F) + ylim(0,1) +
  geom_vline(xintercept = 10) + geom_vline(xintercept = 15,linetype="dashed",color="black")
panels2 <- panels2 + theme_classic() + theme_light() + theme_linedraw() + xlab("Time (years)") + ylab(expression(paste("Proportion of strains that are ",italic(pfhrp2),"-deleted")))
panels2 <- panels2 + theme(axis.title.x = element_text(margin=margin(10)), axis.title.y = element_text(margin=margin(c(0,10))),
                           legend.justification = 'left',legend.position=c(0,0.85),legend.background = element_rect(colour = 'black')) +
  scale_color_manual(labels = c("10%", "25%","60%"),values = gg_color_hue(3)) + annotate("text", x = 7.5, y =0.56125, label = "RDT Introduction")

#fig1c
panels3l <- ggplot(fig1cd,aes(x = PCR.Prevalence,y=hrp2.Deletion.Allele.Frequency)) + geom_point(size=1,alpha=0.25)
panels3l <- panels3l + geom_smooth(se=F,method = "loess") + xlab("PfPR") + ylab(expression(paste("Proportion of strains that are ",italic(pfhrp2),"-deleted"))) + theme_classic() + theme_light() + theme_linedraw()
panels3l <- panels3l + theme(axis.title.x = element_text(margin=margin(10)), axis.title.y = element_text(margin=margin(c(0,10)))) + scale_y_continuous(labels = fmt_dcimals(2)) + ylim(0,1)
panels3l

#fig1d
panels4l <- ggplot(fig1cd,aes(x = PCR.Prevalence,y=Proportion.of.Monoinfections.Due.to.hrp2.Deleted.Strains)) + geom_point(size=1,alpha=0.25)
panels4l <- panels4l + geom_smooth(se=F,method = "loess") + xlab("PfPR") + ylab(expression(paste("Population only infected with  ",italic(pfhrp2),"-deleted mutants"))) + theme_classic() + theme_light() + theme_linedraw()
panels4l <- panels4l + theme(axis.title.x = element_text(margin=margin(10)), axis.title.y = element_text(margin=margin(c(0,10)))) + scale_y_continuous(labels = fmt_dcimals(2)) + ylim(0,1)
panels4l


## PLOTTING
# ------------------------------------------------------------------

### Letter Arrange ###

mp1 <- arrangeGrob(panels1, top = textGrob("a", x = unit(0, "npc")
                                           , y   = unit(1, "npc"), just=c("left","top"),
                                           gp=gpar(col="black", fontsize=18, fontfamily="Times New Roman")))

mp2 <- arrangeGrob(panels2, top = textGrob("b", x = unit(0, "npc")
                                           , y = unit(1, "npc"), just=c("left","top"),
                                           gp=gpar(col="black", fontsize=18, fontfamily="Times New Roman")))

mp3 <- arrangeGrob(panels3l, top = textGrob("c", x = unit(0, "npc")
                                            , y   = unit(1, "npc"), just=c("left","top"),
                                            gp=gpar(col="black", fontsize=18, fontfamily="Times New Roman")))

mp4 <- arrangeGrob(panels4l, top = textGrob("d", x = unit(0, "npc")
                                            , y = unit(1, "npc"), just=c("left","top"),
                                            gp=gpar(col="black", fontsize=18, fontfamily="Times New Roman")))
windows()
gridExtra::grid.arrange(mp1,mp2,mp3,mp4,ncol=2)



# figure 1 - supplement figure 1 plot

# ------------------------------------------------------------------

fig_s1 <- read.csv(system.file("extdata","fig1_s1.csv",package="hrp2malaRia"))
names(fig_s1) <- c("X","hrp2D","prev","Years")
panels <- ggplot(fig_s1,aes(x = Years,y=hrp2D,color=`prev`)) + geom_point()+ ylim(0,1) +
  geom_vline(xintercept = 10) + geom_vline(xintercept = 15,linetype="dashed",color="black")
panels <- panels + theme_classic() + theme_light() + theme_linedraw() + xlab("Time (years)") + ylab("PCR PfPR All Ages")
panels <- panels + theme(axis.title.x = element_text(margin=margin(10)), axis.title.y = element_text(margin=margin(c(0,10))),
                         legend.justification = 'left',legend.position=c(0,0.8),legend.background = element_rect(colour = 'black')) +
  scale_color_manual(labels = c("10%","20%","40%"),values = gg_color_hue(3),name="PCR Prevalence") +
  annotate("text", x = 7.5, y =0.56125, label = "RDT Introduction")

###


mp <- arrangeGrob(panels, top = textGrob("", x = unit(0, "npc")
                                         , y = unit(1, "npc"), just=c("left","top"),
                                         gp=gpar(col="black", fontsize=24, fontfamily="Times New Roman")))

windows()
gridExtra::grid.arrange(mp,ncol=1)



# figure 1 - supplement figure 2 plot

# ------------------------------------------------------------------

fig_s2 <- read.csv(system.file("extdata","fig1_s2.csv",package="hrp2malaRia"))
fig_s2$Comparative.Fitness <- as.factor(fig_s2$Comparative.Fitness)


panel <- ggplot(fig_s2,aes(x = Years,y=hrp2.Deletion.Allele.Frequency,color=Comparative.Fitness)) + geom_smooth(span=0.2,se=F) + ylim(0,1) +
  geom_vline(xintercept = 10) + geom_vline(xintercept = 15,linetype="dashed",color="black")
panel <- panel + theme_classic() + theme_light() + theme_linedraw() + xlab("Time (years)") + ylab(expression(paste("Proportion of strains that are ",italic(pfhrp2),"-deleted")))
panel <- panel + theme(axis.title.x = element_text(margin=margin(10)), axis.title.y = element_text(margin=margin(c(0,10))),
                           legend.justification = 'left',legend.position=c(0,0.8),legend.background = element_rect(colour = 'black')) +
  scale_color_manual(labels = c("100%","97.5","95","92.5",
                                "90%","75%","50%","25%","5%"),values = gg_color_hue(9)) +
  annotate("text", x = 7.5, y =0.56125, label = "RDT Introduction") +
  geom_segment(aes(x = 9.1, y = 0.56, xend = 9.5, yend = 0.56), colour='black', size=1,arrow = arrow(length = unit(0.2, "cm")))


###


mp <- arrangeGrob(panel, top = textGrob("", x = unit(0, "npc")
                                           , y = unit(1, "npc"), just=c("left","top"),
                                           gp=gpar(col="black", fontsize=24, fontfamily="Times New Roman")))

windows()
gridExtra::grid.arrange(mp,ncol=1)


# figure 1 - supplement figure 3 plot

# ------------------------------------------------------------------

fig_s3 <- read.csv(system.file("extdata","fig1_s3.csv",package="hrp2malaRia"))
names(fig_s3) <- c("X","hrp2D","micro","Years")

panel <- ggplot(fig_s3,aes(x = Years,y=hrp2D,color=`micro`)) + geom_smooth(span=0.2,se=F) + ylim(0,1) +
  geom_vline(xintercept = 10) + geom_vline(xintercept = 15,linetype="dashed",color="black")
panel <- panel + theme_classic() + theme_light() + theme_linedraw() + xlab("Time (years)") + ylab(expression(paste("Proportion of strains that are ",italic(pfhrp2),"-deleted")))
panel <- panel + theme(axis.title.x = element_text(margin=margin(10)), axis.title.y = element_text(margin=margin(c(0,10))),
                           legend.justification = 'left',legend.position=c(0,0.8),legend.background = element_rect(colour = 'black')) +
  scale_color_manual(labels = c("0%","30% microscopy",
                                "10% nonadherence","30% microscopy & 10% nonadherence"),values = gg_color_hue(4),
                     name="Microscopy vs Nonadherence") +
  annotate("text", x = 7.5, y =0.56125, label = "RDT Introduction")


###


mp <- arrangeGrob(panel, top = textGrob("", x = unit(0, "npc")
                                        , y = unit(1, "npc"), just=c("left","top"),
                                        gp=gpar(col="black", fontsize=24, fontfamily="Times New Roman")))

windows()
gridExtra::grid.arrange(mp,ncol=1)


# figure 1 - supplement figure 4 plot

# ------------------------------------------------------------------

fig_s4 <- read.csv(system.file("extdata","fig1_s4.csv",package="hrp2malaRia"))
names(fig_s4) <- c("X","hrp2D","nmfv","Years")

panel <- ggplot(fig_s4,aes(x = Years,y=hrp2D,color=`nmfv`)) + geom_smooth(span=0.2,se=F) + ylim(0,1) +
  geom_vline(xintercept = 10) + geom_vline(xintercept = 15,linetype="dashed",color="black")
panel <- panel + theme_classic() + theme_light() + theme_linedraw() + xlab("Time (years)") + ylab(expression(paste("Proportion of strains that are ",italic(pfhrp2),"-deleted")))
panel <- panel + theme(axis.title.x = element_text(margin=margin(10)), axis.title.y = element_text(margin=margin(c(0,10))),
                       legend.justification = 'left',legend.position=c(0,0.8),legend.background = element_rect(colour = 'black')) +
  scale_color_manual(labels = c("0%","80% nmf",
                                "100% nmf","120% nmf"),values = gg_color_hue(4),
                     name="Relative non-hrp2malaRial fever (NMF)") +
  annotate("text", x = 7.5, y =0.56125, label = "RDT Introduction")


###


mp <- arrangeGrob(panel, top = textGrob("", x = unit(0, "npc")
                                        , y = unit(1, "npc"), just=c("left","top"),
                                        gp=gpar(col="black", fontsize=24, fontfamily="Times New Roman")))

windows()
gridExtra::grid.arrange(mp,ncol=1)


# figure 1 - supplement figure 5 plot

# ------------------------------------------------------------------

fig_s5 <- read.csv(system.file("extdata","fig1_s5.csv",package="hrp2malaRia"))
names(fig_s5) <- c("X","hrp2D","level","Years")

panel <- ggplot(fig_s5,aes(x = Years,y=level,color=`nmf`)) + geom_smooth(span=0.2,se=F) + ylim(0,1) +
  geom_vline(xintercept = 10) + geom_vline(xintercept = 15,linetype="dashed",color="black")
panel <- panel + theme_classic() + theme_light() + theme_linedraw() + xlab("Time (years)") + ylab(expression(paste("Proportion of strains that are ",italic(pfhrp2),"-deleted")))
panel <- panel + theme(axis.title.x = element_text(margin=margin(10)), axis.title.y = element_text(margin=margin(c(0,10))),
                           legend.justification = 'left',legend.position=c(0,0.8),legend.background = element_rect(colour = 'black')) +
  scale_color_manual(labels = c("Main Method",
                                "92.5% Comp. Fitness, 40% Microscopy Use, 20% RDT non-adherence, 1.2x NMF Rate",
                                "95% Comp. Fitness, 30% Microscopy Use, 10% RDT non-adherence, 1x NMF Rate",
                                "97.5% Comp. Fitness, 20% Microscopy Use, 5% RDT non-adherence, 0.8x NMF Rate"),
                     values = gg_color_hue(4),
                     name="Model Assumptions") +
  annotate("text", x = 7.5, y =0.56125, label = "RDT Introduction")




###


mp <- arrangeGrob(panel, top = textGrob("", x = unit(0, "npc")
                                        , y = unit(1, "npc"), just=c("left","top"),
                                        gp=gpar(col="black", fontsize=24, fontfamily="Times New Roman")))

windows()
gridExtra::grid.arrange(mp,ncol=1)


