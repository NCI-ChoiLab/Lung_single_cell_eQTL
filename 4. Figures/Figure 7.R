#Date: Jan 23rd 2026
#Author: Oscar Florez Vargas
#Modified: Thong Luong

rm(list=ls(all=TRUE))

## ----------- ##
## -- XCell -- ##
## ----------- ##

library(openxlsx)
library(plater)
library(dplyr)
library(tidyr)
library(plotrix)
library(ggplot2)



setwd('/Users/luongtt/Desktop/xCelligence/')


path_xcell <- "/Users/luongtt/Desktop/xCelligence"

make_df = function(f){
  df_xcell.key <- read.xlsx(paste0(path_xcell, "/",f,".xlsx"), sheet = "Layout")
  df_xcell.key = df_xcell.key[c(1:16),]
  df_xcell.out <- read.xlsx(paste0(path_xcell, "/",f,".xlsx"), sheet = "Cell Index")
  
  
  # format
  df_xcell.key <- df_xcell.key[ complete.cases(df_xcell.key$`Target.Cell`), ]
  row.names(df_xcell.key) <- NULL
  df_xcell.key$Well <- paste0(gsub("[0-9]", "", df_xcell.key$Well),
                              "0",
                              gsub("[A-Z]", "", df_xcell.key$Well))
  
  # re-format: from Excel to CSV
  tmp.fileXLSX <- rbind("", df_xcell.out)
  tmp.fileXLSX[1,1] <- names(tmp.fileXLSX[1])
  names(tmp.fileXLSX)[1] <- "X1"
  
  
  df_plate.time <- data.frame()
  n_count = 1
  
  for ( i in unique(grep("Cell", tmp.fileXLSX$X1, value = TRUE)) ) {
    # i <- "Cell.Index.at:.0:00:00"
    # i <- "Cell Index at: 0:00:11"
    # i <- gsub("[.]", " ", i)
    
    tmp.i <- tmp.fileXLSX[ n_count:(n_count+9), ]
    tmp.i[2,1] <- tmp.i[1,1] # add time
    tmp.i <- cbind(tmp.i, "", "", "", "", "", "")
    names(tmp.i) <- c("time", 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)
    
    tmp.i[2,2:13] <- names(tmp.i[2:13])
    tmp.i <- tmp.i[2:nrow(tmp.i),]
    tmp.i <- rbind(tmp.i, "")
    
    df_plate.time <- rbind(df_plate.time, tmp.i)
    row.names(df_plate.time) <- NULL
    
    n_count = n_count + 10
  }
  
  names(df_plate.time) <- df_plate.time[1,]
  df_plate.time <- df_plate.time[2:nrow(df_plate.time), ]
  row.names(df_plate.time) <- NULL
  df_plate.time[is.na(df_plate.time)] = ''
  
  write.csv(df_plate.time, file = paste0(path_xcell, "/",f,".csv"), row.names = FALSE, quote = F) 
  
  
  
  # remember copy-paste on the CSV file
  fileName <- paste0("/",f,".csv") 
  
  
  xcell_df <- read_plate(
    file = paste0(path_xcell, fileName),
    well_ids_column = "Wells",
    sep = ","
  )
  
  
  
  xcell_df.all <- merge(df_xcell.key, xcell_df, by.x = "Well", "Wells") 
  str(xcell_df.all)
  names(xcell_df.all) <- gsub("Cell.Index.at:.", "", names(xcell_df.all))
  
  time_string <- data.frame(names(xcell_df.all[8:ncol(xcell_df.all)]))
  names(time_string) <- "time"
  time_string$timeID <- time_string$time
  
  time_string <- time_string %>% separate(timeID, into = paste("ID", 1:3, sep = "_"))
  
  names(time_string)[2:4] <- c("hours", "minutes", "seconds")
  time_string$hours <- as.numeric(as.character(time_string$hours))
  time_string$minutes <- as.numeric(as.character(time_string$minutes))
  time_string$seconds <- as.numeric(as.character(time_string$seconds))
  
  time_string$total_minutes <- time_string$hours * 60 + time_string$minutes + time_string$seconds / 60
  time_string$total_minutes <- round(time_string$total_minutes, digits = 2)
  
  names(xcell_df.all)[8:ncol(xcell_df.all)] <- time_string$total_minutes
  
  xcell_df.all <- xcell_df.all %>%
    group_by(Target.Cell,Well.Type) %>%
    mutate(replicate = row_number())
  
  xcell_df.all$replicate <- paste0("R", xcell_df.all$replicate)
  xcell_df.all <- xcell_df.all[,c(length(xcell_df.all), 1:length(xcell_df.all)-1) ]
  xcell_df.all <- xcell_df.all %>% unite("concat", c(1,3), remove = FALSE)
  
  
  # format for plotting and analysis
  xcell_format <- reshape2::melt(xcell_df.all[,c(2,4,10:ncol(xcell_df.all))], id.vars = c('Target.Cell', "replicate"),
                                 variable.name = "time", 
                                 value.name = "index")
  
  xcell_format <- reshape2::dcast(xcell_format, Target.Cell + time ~ replicate,
                                  value.var = "index")
  
  xcell_format$mean <- rowMeans(xcell_format[, grep("R[0-9]", names(xcell_format), value = TRUE) ], na.rm = TRUE)
  xcell_format$sd <- apply(xcell_format[, grep("R[0-9]", names(xcell_format), value = TRUE) ], 1, sd)
  xcell_format$se <- apply(xcell_format[, grep("R[0-9]", names(xcell_format), value = TRUE) ], 1, std.error)
  
  
  tmp.plot <- xcell_format
  tmp.plot$time <- as.numeric(as.character(tmp.plot$time))
  
  # plot
  tmp.plot$time <- tmp.plot$time/60

  
  return(tmp.plot)
}
graph_plot = function(e){
  p = ggplot(e, aes(x=time, y=mean, colour=Target.Cell, group=Target.Cell)) +
    geom_line(aes(colour=Target.Cell), size = 1) +
    # scale_colour_discrete(name = "class") +
    scale_color_manual(values=c("#FBB4AE", "#E41A1C", "#B3CDE3", "#377EB8")) + # red and blue
    labs(x= "Time (hours)",
         y = "Cell Index, xCELLigence",
         colour = 'Condition') +
    
    theme_bw() +
    theme(panel.border = element_rect(linetype = "solid",
                                      colour = "black", size = 0.25),
          legend.position = "right",
          axis.line.x = element_line(size = 0.25, colour = "black"),
          axis.line.y = element_line(size = 0.25, colour = "black"),
          axis.line = element_line(size = 1, color = "black"),
          axis.ticks.length = unit(0.15, "cm"),
          axis.ticks = element_line(size = 0.25, color = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "white",
                                          size = 0, colour = "black"),
          axis.text.x = element_text(colour = "black", size = 12), 
          axis.text.y = element_text(colour = "black", size = 12),
          axis.title.x = element_text(colour = "black", face = "plain", size = 14,
                                      margin = margin(10,0,0,0)),
          axis.title.y = element_text(colour = "black", face = "plain", size = 14,
                                      margin = margin(0,10,0,0)))
  p 
}
graph_plot2 = function(e){
  p= ggplot(e, aes(x=time, y=mean, colour=Target.Cell, group=Target.Cell)) +
    geom_line(aes(colour=Target.Cell), size = 1) +
    geom_ribbon(aes(ymin = mean - se, ymax = mean + se, fill = Target.Cell), size = .1, show.legend = F, alpha = .8) +
    # scale_colour_discrete(name = "class") +
    scale_color_manual(values=c("#FBB4AE", "#E41A1C", "#B3CDE3", "#377EB8")) + # red and blue
    labs(x= "Time (hours)",
         y = "Cell Index, xCELLigence",
         colour = 'Condition') +
    
    theme_bw() +
    theme(panel.border = element_rect(linetype = "solid",
                                      colour = "black", size = 0.25),
          legend.position = "right",
          axis.line.x = element_line(size = 0.25, colour = "black"),
          axis.line.y = element_line(size = 0.25, colour = "black"),
          axis.line = element_line(size = 1, color = "black"),
          axis.ticks.length = unit(0.15, "cm"),
          axis.ticks = element_line(size = 0.25, color = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "white",
                                          size = 0, colour = "black"),
          axis.text.x = element_text(colour = "black", size = 12), 
          axis.text.y = element_text(colour = "black", size = 12),
          axis.title.x = element_text(colour = "black", face = "plain", size = 14,
                                      margin = margin(10,0,0,0)),
          axis.title.y = element_text(colour = "black", face = "plain", size = 14,
                                      margin = margin(0,10,0,0)))
  p 
}
exp1 = make_df('Attempt_1')
graph_plot(exp1)

exp2 = make_df('Attempt_2')
graph_plot(exp2)
graph_plot2(exp2)

exp3 = make_df('Attempt_3')
graph_plot(exp3)
graph_plot2(exp3)

exp4 = make_df('Attempt_4')
graph_plot(exp4)

exp5 = make_df('Attempt_5')
graph_plot(exp5)
graph_plot2(exp5)


exp2 = subset(exp2, time < 100)
graph_plot2(exp2)
exp3 = subset(exp3, time < 100)
graph_plot2(exp3)
exp5 = subset(exp5, time < 100)
graph_plot2(exp5)

exp2.1 = exp2[,c("Target.Cell","time","R1",'R2','R3','R4')]
exp3.1 = exp3[,c("Target.Cell","time","R1",'R2','R3','R4')]
exp5.1 = exp5[,c("Target.Cell","time","R1",'R2','R3','R4')]

colnames(exp3.1)[3:6] = c("R5",'R6','R7','R8')
colnames(exp5.1)[3:6] = c("R9",'R10','R11','R12')

t = merge(exp2.1,exp3.1, by = c('Target.Cell','time'), all = T)
t = merge(t,exp5.1, by = c('Target.Cell','time'), all = T)
t$mean <- rowMeans(t[, grep("R[0-12]", names(t), value = TRUE) ], na.rm = TRUE)
t$sd <- apply(t[, grep("R[0-12]", names(t), value = TRUE) ], 1, sd)
t$se <- apply(t[, grep("R[0-12]", names(t), value = TRUE) ], 1, std.error)
graph_plot(t)



## ----- statistical analysis ----- ##
## -------------------------------- ##

library(nlme)

status.pairs <- data.frame(pair1 = c("pLKO.5", "shRNA_43"),
                           pair2 = c("pLKO.5", "shRNA_46"),
                           pair3 = c("pLKO.5", "shRNA_48"))



stat_test = function(df){
  lme_df.byStatus <- data_frame()
  for ( pair in names(status.pairs) ) {
    df2 = df
    df2 = subset(df2, time > 24)
    df2$time = df2$time*60
    
    model.tmp <- df2[ which(df2$Target.Cell %in% status.pairs[[pair]] ), ]
    names(model.tmp)
    
    model.tmp <- reshape2::melt(model.tmp, id.vars = c("mean", "sd", "se", "time", "Target.Cell"),
                                variable.name = "ID", 
                                value.name = "count")
    
    table(model.tmp$ID)
    
    ## Effect of treatment on cell count trajectory
    fita <- lme(count ~ time,
                random = ~ 1 | ID, method = "ML", data=model.tmp)
    lme.tmp <- summary(fita)
    print(lme.tmp$tTable)
    #sjPlot::tab_model(fita)
    
    fitb <- lme(count ~ time + Target.Cell,
                random = ~ 1 | ID, method = "ML", data=model.tmp)
    lme.tmp <- summary(fitb)
    print(lme.tmp$tTable)
    # sjPlot::tab_model(fitb)
    
    # fitc <- lme(count ~ time * status,
    #             random = ~ 1 | ID, method = "ML", data=model.tmp)
    # lme.tmp <- summary(fitc)
    # print(lme.tmp$tTable)
    # # sjPlot::tab_model(fitc)
    
    # fitb <- lme(count ~ time * status,
    #             random = ~ 1 | ID, method = "ML", data=model.tmp)
    # lme.tmp <- summary(fitb)
    # print(lme.tmp$tTable)
    # # sjPlot::tab_model(fitb)
    
    anova.fit <- data.frame(anova(fita,fitb))
    anova.fit
    
    # create temporary data frame
    df.lm <- data.frame(cell_line = 'A549',
                        pair = paste0(status.pairs[,pair][1], " and ", status.pairs[,pair][2]),
                        
                        b_fita_time = coef(summary(fita))[2,1],
                        SE_fita_time = coef(summary(fita))[2,2],
                        P_fita_time = coef(summary(fita))[2,5],
                        
                        b_fitb_time = coef(summary(fitb))[2,1],
                        SE_fitb_time = coef(summary(fitb))[2,2],
                        P_fitb_time = coef(summary(fitb))[2,5],
                        
                        b_fitb_trea = coef(summary(fitb))[3,1],
                        SE_fitb_trea = coef(summary(fitb))[3,2],
                        P_fitb_trea = coef(summary(fitb))[3,5],

                        
                        model_fita_AIC = anova.fit[1,4],
                        model_fitb_AIC = anova.fit[2,4],
                        model_fita_BIC = anova.fit[1,5],
                        model_fitb_BIC = anova.fit[2,5],
                        P_anova = anova.fit[2,9],
                        
                        stringsAsFactors = FALSE)
    
    lme_df.byStatus <- rbind(lme_df.byStatus, df.lm)
    
  }
  return(lme_df.byStatus)
}

exp2_stat = stat_test(exp2)
exp3_stat = stat_test(exp3)
exp5_stat = stat_test(exp5)

# # for Experiment 1
# write.csv(lme_df.byStatus, file = "LME_xCELLigence_byStatus_FBS.experiment_time_all_01162023.csv", row.names = FALSE)
v1 = c(exp2_stat[1,11],exp3_stat[1,11],exp5_stat[1,11])
v2 = c(exp2_stat[2,11],exp3_stat[2,11],exp5_stat[2,11])
v3 = c(exp2_stat[3,11],exp3_stat[3,11],exp5_stat[3,11])
fisherIntegration  <- function (vector){
  my_length=length(vector)
  deg_free=my_length*2
  y=-2*sum(log(vector))
  p.val <- 1-pchisq(y, df = deg_free);
  p.val=as.numeric(p.val);
  return(p.val)
}

graph_cellindex = function(e){
  ggplot(e, aes(x=time, y=mean, colour=Target.Cell, group=Target.Cell)) +
  geom_line(aes(colour=Target.Cell), size = .5) +
  #geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, fill = Target.Cell), size = .5) +
  # scale_colour_discrete(name = "class") +
  scale_color_manual(values=c("#FBB4AE", "#E41A1C", "#B3CDE3", "#377EB8")) + # red and blue
  labs(x= "Time (hours)",
       y = "Cell Index, xCELLigence",
       colour = 'Condition') +
  
  # scale_y_continuous(breaks = seq(0,4,0.5),
  #                    labels = seq(0,4,0.5),
  #                    limits = c(-0.05,4)) +
  # scale_x_continuous(breaks = seq(0,300,50),
  #                    labels = seq(0,300,50),
  #                    limits = c(0,300)) +
  # scale_y_log10(breaks = 10,
  #               expand = c(0, 0),
  #               limits = c(0.02, 100)) +
  # scale_y_log10() +
  
  theme_bw() +
  theme(panel.border = element_rect(linetype = "solid",
                                    colour = "black", size = 0.25),
        legend.position = "right",
        axis.line.x = element_line(size = 0.25, colour = "black"),
        axis.line.y = element_line(size = 0.25, colour = "black"),
        axis.line = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.ticks = element_line(size = 0.25, color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white",
                                        size = 0, colour = "black"),
        axis.text.x = element_text(colour = "black", size = 12), 
        axis.text.y = element_text(colour = "black", size = 12),
        axis.title.x = element_text(colour = "black", face = "plain", size = 14,
                                    margin = margin(10,0,0,0)),
        axis.title.y = element_text(colour = "black", face = "plain", size = 14,
                                    margin = margin(0,10,0,0)))
}




