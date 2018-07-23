library(shiny)
library(ggplot2)
library(dplyr)
library(markdown)

# ggplot prep

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

mytheme <- ggplot2::theme_classic() + ggplot2::theme_light() +
  ggplot2::theme(axis.text = ggplot2::element_text(size = 12, family = "Calibri"),
                 axis.title.y = ggplot2::element_text(margin=ggplot2::margin(c(0,10)),size = 12, family = "Calibri"),
                 axis.title.x = ggplot2::element_text(margin=ggplot2::margin(c(10)),size = 12, family = "Calibri"),
                 legend.text = ggplot2::element_text(size = 12, family = "Calibri"),
                 legend.title = ggplot2::element_text(size = 12, family = "Calibri"),
                 plot.title = ggplot2::element_text(size = 14, family = "Calibri",hjust = 0.5),
                 axis.line = ggplot2::element_line(),
                 strip.background = ggplot2::element_rect(fill = "white",color = "grey",size = .5),
                 strip.text.x = ggplot2::element_text(colour = 'black', size = 10, family = "Calibri"),
                 strip.text.y = ggplot2::element_blank(),
                 panel.spacing = unit(1, "lines")
  )

# data transform prep
dat <- readRDS("administrative_covariates.rds")
mal_dat <- which(dat$PfPR_2015>0 & dat$Treat_2015>0)
dat <- dat[mal_dat,]
admin_res <- readRDS("admin_final_stats.rds")

## add back in these columns that were removed to save upload space
for(i in 1:598){
  
  admin_res[[i]]$raw$sample <- c(rep("All",2400),rep("Under 5s",2400))
  admin_res[[i]]$raw$rep <- as.character(c(rep(1:200,24)))
  admin_res[[i]]$raw$include <- c(admin_res[[i]]$raw$relativelower < 0 & admin_res[[i]]$raw$relativeupper > 0)
  
}

#plotting function prep

fig1 <- function(country,admin){
  
  i <- which(dat$Country==country & dat$NAME == admin)
  # true_df_sum <- group_by(admin_res[[i]]$raw,`8-week period`,sample) %>% summarise(PointEst=mean(relative),
  #                                                                                  Lower=mean(relativelower),
  #                                                                                  Upper=mean(relativeupper),
  #                                                                                  inc=mean(inc))
  # true_df_sum <- group_by(admin_res[[i]]$raw,`8-week period`,sample) %>% summarise(PointEst=mean(relative),
  #                                                                                  Lower=mean(relative)-(2*sd(relative)),
  #                                                                                  Upper=mean(relative)+(2*sd(relative)),
  #                                                                                  inc=mean(inc))
  true_df_sum <- group_by(admin_res[[i]]$raw,`8-week period`,sample) %>% summarise(PointEst=mean(relative),
                                                                                   Lower=t.test(relativelower)$conf.int[1],
                                                                                   Upper=t.test(relativeupper)$conf.int[2])
  
  
  true_df2 <- group_by(true_df_sum,sample,`8-week period`) %>% summarise(mean=mean(PointEst),low=mean(Lower),high=mean(Upper))
  melt_good_true <- reshape2::melt(true_df2,id.vars=c("sample","8-week period","low","high"),measure.vars=c("mean"))
  names(melt_good_true)[2] <- "8-week period"
  melt_good_true[,c(3,4,6)] <- melt_good_true[,c(3,4,6)]*100
  
  
  submelt <- melt_good_true[melt_good_true$`8-week period`!="All-year",]
  
  ggplot(data = submelt,aes(x=`8-week period`,y=value, color = `8-week period` )) + 
    geom_hline(yintercept = 0,col="red",linetype="dashed") +
    geom_pointrange(aes(x=`8-week period`,y=value,ymax=high,ymin=low,color=`8-week period`),
                    position=position_dodge(width=0.75)) +
    geom_tile(data = submelt[best_months(country,admin),],aes(height=(high-low)*1.2) ,col="black",fill=NA,lwd=1) +
    facet_wrap(~sample) +
    scale_color_manual(name="Eight Week Interval",
                       values = c(gg_color_hue(length(levels(melt_good_true$`8-week period`))))) +
    xlab("Individuals Sampled") +
    ylab(expression(paste("Bias Microscopy +ve / RDT -ve due to ",italic(pfhrp2),"-deletions"))) + 
    ggtitle(paste0(dat$NAME[i],", ",dat$Country[i]))  + mytheme +
    theme(axis.text.x = element_text(angle=45,vjust = 0.5),
          axis.title.x = element_blank())
  
}

fig2 <- function(country,admin,month_pair){
  
  i <- which(dat$Country==country & dat$NAME == admin)
  
  ggplot(data = admin_res[[i]]$raw[admin_res[[i]]$raw$`8-week period`==month_pair,],aes(x=rep,y=relative*100, color = include )) + 
    geom_pointrange(aes(x=rep,y=relative*100,ymax=relativeupper*100,ymin=relativelower*100,color=include),
                    position=position_dodge(width=0.75)) +
    facet_wrap(~sample)+
    geom_hline(yintercept = 0,col="black",linetype="dashed",cex=1) +
    scale_color_manual(name="Confidence Interval\nIncludes True",
                       values = c(gg_color_hue(2))) +
    xlab("Simulation Repetitions") +
    ylab(expression(paste("Bias Microscopy +ve / RDT -ve due to ",italic(pfhrp2),"-deletions"))) + 
    mytheme + 
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          plot.title = element_text(hjust=0.5),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank()
    ) +
    ggtitle(month_pair)
  
  
}

best_months <- function(country, admin, month=FALSE){
  i <- which(dat$Country==country & dat$NAME == admin)
  true_df_sum <- group_by(admin_res[[i]]$raw,`8-week period`,sample) %>% summarise(PointEst=mean(relative),
                                                                                   Lower=t.test(relativelower)$conf.int[1],
                                                                                   Upper=t.test(relativeupper)$conf.int[2])
  true_df2 <- group_by(true_df_sum,sample,`8-week period`) %>% summarise(mean=mean(PointEst),low=mean(Lower),high=mean(Upper))
  melt_good_true <- reshape2::melt(true_df2,id.vars=c("sample","8-week period","low","high"),measure.vars=c("mean"))
  names(melt_good_true)[2] <- "8-week period"
  melt_good_true[,c(3,4,6)] <- melt_good_true[,c(3,4,6)]*100
  
  all <- which.min(abs(melt_good_true$value[1:12]))
  fives <- which.min(abs(melt_good_true$value[13:24]))+12
  if(month){
    all <- as.character(melt_good_true$`8-week period`[all])
    fives <- as.character(melt_good_true$`8-week period`[fives])
  }
  return(c(all,fives))
}

table3 <- function(country,admin,month_pair){
  
  i <- which(dat$Country==country & dat$NAME == admin)
  
  l <- length(admin_res[[i]]$raw[admin_res[[i]]$raw$`8-week period`==month_pair &
                                   admin_res[[i]]$raw$sample=="All"  ,]$relativelower)
  
  # above_all <- sum(admin_res[[i]]$raw[admin_res[[i]]$raw$`8-week period`==month_pair &
  #                                       admin_res[[i]]$raw$sample=="All"  ,]$relativelower>0)
  # above_kids <- sum(admin_res[[i]]$raw[admin_res[[i]]$raw$`8-week period`==month_pair &
  #                                        admin_res[[i]]$raw$sample=="Under 5s"  ,]$relativelower>0)
  # below_all <- sum(admin_res[[i]]$raw[admin_res[[i]]$raw$`8-week period`==month_pair &
  #                                       admin_res[[i]]$raw$sample=="All"  ,]$relativeupper<0)
  # below_kids <- sum(admin_res[[i]]$raw[admin_res[[i]]$raw$`8-week period`==month_pair &
  #                                        admin_res[[i]]$raw$sample=="Under 5s"  ,]$relativeupper<0)
  correct_all <- sum(admin_res[[i]]$raw[admin_res[[i]]$raw$`8-week period`==month_pair &
                                          admin_res[[i]]$raw$sample=="All"  ,]$include)
  correct_kids <- sum(admin_res[[i]]$raw[admin_res[[i]]$raw$`8-week period`==month_pair &
                                           admin_res[[i]]$raw$sample=="Under 5s"  ,]$include)
  
  bms <- best_months(country,admin,TRUE)
  best_all <- sum(admin_res[[i]]$raw[admin_res[[i]]$raw$`8-week period`==bms[1] &
                                          admin_res[[i]]$raw$sample=="All"  ,]$include)
  best_kids <- sum(admin_res[[i]]$raw[admin_res[[i]]$raw$`8-week period`==bms[2] &
                                           admin_res[[i]]$raw$sample=="Under 5s"  ,]$include)
  
  df <- c(100*correct_all/l,100*correct_kids/l,100*best_all/l,100*best_kids/l)
  df
}

admin_country_pause <- function(country,admin){
  
  ads <- dat$NAME[which(dat$Country==country)]
  if(is.element(admin,ads)){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

ui <- shinyUI(
  fluidPage(
  includeCSS("style.css"),
  br(),
  sidebarPanel(width = 3,
               h3('Seasonal Impacts on pfhrp2-deletions'),
               br(),
               uiOutput("firstSelection"),
               uiOutput("secondSelection"),
               uiOutput("thirdSelection"),
               hr(),
               uiOutput("tableheader"),
               hr(id="a"),
               uiOutput("table"),
               hr(),
               uiOutput("table2"),
               hr(),
               uiOutput("tableheader2"),
               hr(id="a"),
               uiOutput("table3"),
               hr(),
               uiOutput("table4")
  ),
  mainPanel(
    tabsetPanel(
      tabPanel("Description",
               includeMarkdown('blurb.md')
      ),
      tabPanel("Plots",
               br(),
               plotOutput('plot1', height = "350px"),
               plotOutput('plot2', height = "350px")
      )
    )
  )
)
)

server <- function(input, output) {
  
  output$firstSelection <- renderUI({
    selectInput('country', 'Country', dat$Country %>% unique %>% sort, 
                selected="Angola")
  })
  
  
  output$secondSelection <- renderUI({
    selectInput("admin", "Administrative Unit", dat$NAME[which(dat$Country==input$country)],
                selected = "Cabinda")
  })
  
  output$thirdSelection <- renderUI({
    selectInput('month', '8-week Interval',
                c("Jan-Mar","Feb-Apr","Mar-May","Apr-Jun","May-Jul",
                  "Jun-Aug","Jul-Sep","Aug-Oct","Sep-Nov","Oct-Dec",
                  "Nov-Jan","Dec-Feb"),
                selected = "Jan-Mar")
  })
  
  output$plot1 <- renderPlot({
    if(!is.null(input$admin)){
      if(admin_country_pause(input$country,input$admin)){
        fig1(input$country,input$admin)
      }
    }
  })
  
  output$plot2 <- renderPlot({
    if(!is.null(input$admin)){
      if(admin_country_pause(input$country,input$admin)){
        fig2(input$country,input$admin,input$month)
      }
    }
  })
  
  output$tableheader <- renderUI({
    if(!is.null(input$admin)){
      if(admin_country_pause(input$country,input$admin)){
        tagList(h5(paste("Sampling Accuracy for",input$month)))
      }
    }
  })
  
  output$table <- renderUI({
    if(!is.null(input$admin)){
      if(admin_country_pause(input$country,input$admin)){
        tagList(h5(paste("Correct Prediction:")),
                h5(paste("\n",table3(input$country,input$admin,input$month)[1],"%")
        ))
      }
    }
  })
  
  output$table2 <- renderUI({
    if(!is.null(input$admin)){
      if(admin_country_pause(input$country,input$admin)){
        tagList(h5(paste("Correct Prediction Under 5s:")),
                h5(paste("\n",table3(input$country,input$admin,input$month)[2],"%")
                ))
      }
    }
  })
  
  output$tableheader2 <- renderUI({
    if(!is.null(input$admin)){
      if(admin_country_pause(input$country,input$admin)){
        tagList(h5(paste("Best 8-week Interval for ",input$admin)))
      }
    }
  })
  
  output$table3 <- renderUI({
    if(!is.null(input$admin)){
      if(admin_country_pause(input$country,input$admin)){
        tagList(h5(paste0("Whole Population: ",best_months(input$country,input$admin,TRUE)[1])),
                h5(paste0(table3(input$country,input$admin,input$month)[3],"% Correct Prediction")
        ))
      }
    }
  })
  
  output$table4 <- renderUI({
    if(!is.null(input$admin)){
      if(admin_country_pause(input$country,input$admin)){
        tagList(h5(paste0("Under 5s: ",best_months(input$country,input$admin,TRUE)[2])),
                h5(paste0(table3(input$country,input$admin,input$month)[4],"% Correct Prediction")
                ))
      }
    }
  })
  
  
}


shinyApp(ui, server)
