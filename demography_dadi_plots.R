library("data.table")
library("dplyr")
library("tidyr")
library("ggplot2")
library(gridExtra)
library(openxlsx)



### Demographic parameters ###
p <- vector('list', 24)
mig <- vector('list', 24)
pre <- vector('list', 24)
i=0
wb <- createWorkbook()

for (POSSEL in c(0, 2, 20, 200))
{
  for (MIGRATION in c(0, 0.01))
  {
    for (CHANGE in c(0, 0.05, 0.25))
    {
      i = i+1
      print(paste(i, MIGRATION, POSSEL, CHANGE))
      dt = fread(paste0("../results/dadi_outputs/sim_", MIGRATION, "-", POSSEL, "-", CHANGE, "_demographic_confidence_intervals.txt"))
      colnames(dt) <- c("bound", "step_size", "Tpre", "nuPre", "Tdiv", "nu1div", "nu2div", "T1F", "T2F", "nu1F", "nu2F", "mw2d", "md2w", "missid")

      dt_step = data.frame(t(subset(dt, step_size == 1e-4)))
      #dt_step = subset(dt_step, X1 != "Lower_bound")
      
      dt_step$X1 = as.numeric(dt_step$X1)
      dt_step$X2 = as.numeric(dt_step$X2)

      dt_step$mean = ((dt_step$X1)+(dt_step$X2))/2
      
      max_TF = max(dt_step[c(8,9),]$mean)

      dt_step[dt_step < 0] <- 0

      wild_Ne_vector =   rev(c(dt_step[10,]$mean,dt_step[10,]$mean, dt_step[6,]$mean, dt_step[4,]$mean,1))
      dom_Ne_vector =    rev(c(dt_step[11,]$mean,dt_step[11,]$mean, dt_step[7,]$mean, dt_step[4,]$mean,1))

      wild_Ne_vector_low =   rev(c(dt_step[10,]$X1,dt_step[10,]$X1, dt_step[6,]$X1, dt_step[4,]$X1,1))
      wild_Ne_vector_up  =   rev(c(dt_step[10,]$X2,dt_step[10,]$X2, dt_step[6,]$X2, dt_step[4,]$X2,1))

      dom_Ne_vector_low  =    rev(c(dt_step[11,]$X1,dt_step[11,]$X1, dt_step[7,]$X1, dt_step[4,]$X1,1))
      dom_Ne_vector_up   =    rev(c(dt_step[11,]$X2,dt_step[11,]$X2, dt_step[7,]$X2, dt_step[4,]$X2,1))

      wild_Ne_true_vector = c(1,1)
      wild_time_true_vector = c(0,1)

      if (dt_step[8,]$mean == max_TF) {
        wild_time_vector = rev(c(0,dt_step[8,]$mean,  dt_step[8,]$mean+dt_step[5,]$mean,  dt_step[8,]$mean+dt_step[5,]$mean+dt_step[3,]$mean, (dt_step[8,]$mean+dt_step[5,]$mean+dt_step[3,]$mean)*1.5))
        wild_time_vector_low = rev(c(0,dt_step[8,]$X1,  dt_step[8,]$X1+dt_step[5,]$X1,  dt_step[8,]$X1+dt_step[5,]$X1+dt_step[3,]$X1, (dt_step[8,]$X1+dt_step[5,]$X1+dt_step[3,]$X1)*1.5))
        wild_time_vector_up  = rev(c(0,dt_step[8,]$X2,  dt_step[8,]$X2+dt_step[5,]$X2,  dt_step[8,]$X2+dt_step[5,]$X2+dt_step[3,]$X2, (dt_step[8,]$X2+dt_step[5,]$X2+dt_step[3,]$X2)*1.5))
        
        dom_time_vector =  rev(c(0,dt_step[9,]$mean,  dt_step[9,]$mean+dt_step[5,]$mean+(dt_step[8,]$mean-dt_step[9,]$mean), dt_step[9,]$mean+dt_step[5,]$mean+(dt_step[8,]$mean-dt_step[9,]$mean)+dt_step[3,]$mean,(dt_step[9,]$mean+dt_step[5,]$mean+(dt_step[8,]$mean-dt_step[9,]$mean)+dt_step[3,]$mean)*1.5))
        dom_time_vector_low =  rev(c(0,dt_step[9,]$X1,  dt_step[9,]$X1+dt_step[5,]$X1+(dt_step[8,]$X1-dt_step[9,]$X1), dt_step[9,]$X1+dt_step[5,]$X1+(dt_step[8,]$X1-dt_step[9,]$X1)+dt_step[3,]$X1,(dt_step[9,]$X1+dt_step[5,]$X1+(dt_step[8,]$X1-dt_step[9,]$X1)+dt_step[3,]$X1)*1.5))
        dom_time_vector_up  =  rev(c(0,dt_step[9,]$X2,  dt_step[9,]$X2+dt_step[5,]$X2+(dt_step[8,]$X2-dt_step[9,]$X2), dt_step[9,]$X2+dt_step[5,]$X2+(dt_step[8,]$X2-dt_step[9,]$X2)+dt_step[3,]$X2,(dt_step[9,]$X2+dt_step[5,]$X2+(dt_step[8,]$X2-dt_step[9,]$X2)+dt_step[3,]$X2)*1.5))
        
        } else {
        wild_time_vector = rev(c(0,dt_step[8,]$mean,  dt_step[8,]$mean+dt_step[5,]$mean+(dt_step[9,]$mean-dt_step[8,]$mean), dt_step[8,]$mean+dt_step[5,]$mean+(dt_step[9,]$mean-dt_step[8,]$mean)+dt_step[3,]$mean,(dt_step[8,]$mean+dt_step[5,]$mean+(dt_step[9,]$mean-dt_step[8,]$mean)+dt_step[3,]$mean)*1.5))
        wild_time_vector_low = rev(c(0,dt_step[8,]$X1,  dt_step[8,]$X1+dt_step[5,]$X1+(dt_step[9,]$X1-dt_step[8,]$X1), dt_step[8,]$X1+dt_step[5,]$X1+(dt_step[9,]$X1-dt_step[8,]$X1)+dt_step[3,]$X1,(dt_step[8,]$X1+dt_step[5,]$X1+(dt_step[9,]$X1-dt_step[8,]$X1)+dt_step[3,]$X1)*1.5))
        wild_time_vector_up  = rev(c(0,dt_step[8,]$X2,  dt_step[8,]$X2+dt_step[5,]$X2+(dt_step[9,]$X2-dt_step[8,]$X2), dt_step[8,]$X2+dt_step[5,]$X2+(dt_step[9,]$X2-dt_step[8,]$X2)+dt_step[3,]$X2, (dt_step[8,]$X2+dt_step[5,]$X2+(dt_step[9,]$X2-dt_step[8,]$X2)+dt_step[3,]$X2)*1.5))
        
        dom_time_vector  = rev(c(0,dt_step[9,]$mean,  dt_step[9,]$mean+dt_step[5,]$mean,  dt_step[9,]$mean+dt_step[5,]$mean+dt_step[3,]$mean, ( dt_step[9,]$mean+dt_step[5,]$mean+dt_step[3,]$mean)*1.5))
        dom_time_vector_low  = rev(c(0,dt_step[9,]$X1,  dt_step[9,]$X1+dt_step[5,]$X1,  dt_step[9,]$X1+dt_step[5,]$X1+dt_step[3,]$X1, (dt_step[9,]$X1+dt_step[5,]$X1+dt_step[3,]$X1)*1.5))
        dom_time_vector_up   = rev(c(0,dt_step[9,]$X2,  dt_step[9,]$X2+dt_step[5,]$X2,  dt_step[9,]$X2+dt_step[5,]$X2+dt_step[3,]$X2, (dt_step[9,]$X2+dt_step[5,]$X2+dt_step[3,]$X2)*1.5))
      }

      yLIM = max(c(dom_Ne_vector,wild_Ne_vector))+0.1
      xLIM = max(c(dom_time_vector,wild_time_vector))

      dom_Ne_true_vector =   rev(c(1,   1,   0.04,1   ,1))
      dom_time_true_vector = rev(c(xLIM,0.20,0.20,0.18,0))
      
      wild = cbind(wild_time_vector, wild_Ne_vector, "Wild", "Inferred")
      dom = cbind(dom_time_vector, dom_Ne_vector, "Domestic", "Inferred")
      dom2 = cbind(dom_time_true_vector, dom_Ne_true_vector, "Domestic", "True")
      dt_plot = data.frame(rbind(wild, dom, dom2))
      colnames(dt_plot) <- c("Time", "Ne", "Population", "State")
      dt_plot$Time = as.numeric(dt_plot$Time)
      dt_plot$Ne = as.numeric(dt_plot$Ne)
      
      dt_plot$State <- factor(dt_plot$State, levels = c("Inferred", "True"), ordered = TRUE)
      
      p[[i]] <- ggplot(subset(dt_plot), aes(x=Time, y=Ne, color=paste0(Population, "-", State))) +
        geom_step(direction = "vh", size=1.5, alpha=0.6) +
        scale_color_manual(values = c("#1b9e77", "black", "#d95f02")) +
        theme_bw() +
        xlab("") + ylab("") +
        theme(legend.position="none", plot.title = element_text(hjust = 0.5, size = 8)) +
        ggtitle(paste0("m = ", MIGRATION, ", Sb = ", POSSEL/2, " and pc = ", CHANGE))
    
      dt_step = dt_step[-c(1,2),]
      colnames(dt_step) <- c("CI_2.5%", "CI_97.5%", "value")
      
      sheet_name = paste0("SCENARIO_ID=", i)
      addWorksheet(wb,sheet_name)
      writeData(wb, sheet_name, dt_step, rowNames=TRUE)
      
      dt_step$params = rownames(dt_step)

      ### Pre Ne plots ###
      dt_step_pre <- subset(dt_step, params == "Tpre" | params == "nuPre") %>% 
        mutate(epoch = case_when(`CI_2.5%` == 0 & params == "Tpre" ~  "Not Significant",
                                (`CI_2.5%` <= 1 & `CI_97.5%` >= 1) & params == "nuPre" ~  "Not Significant"))
      dt_step_pre["epoch"][is.na(dt_step_pre["epoch"])] <- "Significant"

      pre[[i]] <- ggplot(dt_step_pre, aes(x=params, y=value, color=epoch)) +
        geom_point() +
        geom_errorbar(aes(ymin=`CI_2.5%`, ymax=`CI_97.5%`), width=.2,
                      position=position_dodge(0.05)) +
        #ylim(0,2) +
        xlab("") + ylab("") +
        facet_wrap(~params, scales = "free") +
        scale_color_manual(breaks = c("Significant", "Not Significant"),
                           values=c("red", "grey")) +
        theme_bw() +
        theme(legend.position="none", plot.title = element_text(hjust = 0.5, size = 8)) +
        ggtitle(paste0("m = ", MIGRATION, ", Sb = ", POSSEL/2, " and pc = ", CHANGE))
      
      
      ### Migration plots ###
      dt_step_mig <- subset(dt_step, params == "md2w" | params == "mw2d") %>% 
        mutate(epoch = case_when(`CI_2.5%` == 0 & params == "md2w" ~  "Not Significant",
                                 `CI_2.5%` == 0 & params == "mw2d" ~  "Not Significant"))
      dt_step_mig["epoch"][is.na(dt_step_mig["epoch"])] <- "Significant"
      
      mig[[i]] <- ggplot(dt_step_mig, aes(x=params, y=value, color=epoch)) +
        geom_point() +
        geom_errorbar(aes(ymin=`CI_2.5%`, ymax=`CI_97.5%`), width=.2,
                      position=position_dodge(0.05)) +
        #ylim(0,2) +
        xlab("") + ylab("") +
        #facet_wrap(~params) +
        scale_color_manual(breaks = c("Significant", "Not Significant"),
                           values=c("red", "grey")) +
        theme_bw() +
        geom_hline(yintercept = 0.01) +
        geom_hline(yintercept = 0, color = "white") +
        
        theme(legend.position="none", plot.title = element_text(hjust = 0.5, size = 8)) +
        ggtitle(paste0("m = ", MIGRATION, ", Sb = ", POSSEL/2, " and pc = ", CHANGE))
    }
  }
}

saveWorkbook(wb, "../results/Supp_Table1.xlsx", overwrite = TRUE)

all = grid.arrange(p[[1]], p[[2]], p[[3]], p[[4]], p[[5]], p[[6]], p[[7]], p[[8]], p[[9]], p[[10]], p[[11]], p[[12]], p[[13]], p[[14]], p[[15]], p[[16]], p[[17]], p[[18]], p[[19]], p[[20]], p[[21]], p[[22]], p[[23]], p[[24]], ncol = 3)
ggsave(all, filename = "../plots/Figure2.png",width = 8.5, height = 12, dpi = 300)

all = grid.arrange(pre[[1]], pre[[2]], pre[[3]], pre[[4]], pre[[5]], pre[[6]], pre[[7]], pre[[8]], pre[[9]], pre[[10]], pre[[11]], pre[[12]], pre[[13]], pre[[14]], pre[[15]], pre[[16]], pre[[17]], pre[[18]], pre[[19]], pre[[20]], pre[[21]], pre[[22]], pre[[23]], pre[[24]], ncol = 3)
ggsave(all, filename = "../plots/Supp_Figure3.png",width = 8.5, height = 12, dpi = 300)

all = grid.arrange(mig[[1]], mig[[2]], mig[[3]], mig[[4]], mig[[5]], mig[[6]], mig[[7]], mig[[8]], mig[[9]], mig[[10]], mig[[11]], mig[[12]], mig[[13]], mig[[14]], mig[[15]], mig[[16]], mig[[17]], mig[[18]], mig[[19]], mig[[20]], mig[[21]], mig[[22]], mig[[23]], mig[[24]], ncol = 3)
ggsave(all, filename = "../plots/Supp_Figure4.png",width = 8.5, height = 12, dpi = 300)

