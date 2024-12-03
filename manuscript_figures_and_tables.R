library("data.table")
library("dplyr")
library("tidyr")
library("ggplot2")
library(gridExtra)
library(openxlsx)


### MAIN FIGURE 2 AND SUPP FIGURES 2 AND 3 ###

p <- vector('list', 24)
mig <- vector('list', 24)
pre <- vector('list', 24)
i=0
wb <- createWorkbook()

for (MIGRATION in c(0, 0.01))
{
  for (POSSEL in c(0, 2, 20, 200))
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
      #if (i == 13){
      #  dt_step$mean = c(NA, 0.0001, 0.073896373549883,	1.70449312822617,	0.511609212705555,	0.710288356088611,	0.505914407426143,	0.445655843912499,	0.421018882179885,	4.57471975300398,	4.65693251839453,	0.0147527753761547,	0.0109875137620283,	0.00247505568469394)
      #}
      
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
ggsave(all, filename = "../plots/Supp_Figure2.png",width = 8.5, height = 12, dpi = 300)

all = grid.arrange(mig[[1]], mig[[2]], mig[[3]], mig[[4]], mig[[5]], mig[[6]], mig[[7]], mig[[8]], mig[[9]], mig[[10]], mig[[11]], mig[[12]], mig[[13]], mig[[14]], mig[[15]], mig[[16]], mig[[17]], mig[[18]], mig[[19]], mig[[20]], mig[[21]], mig[[22]], mig[[23]], mig[[24]], ncol = 3)
ggsave(all, filename = "../plots/Supp_Figure3.png",width = 8.5, height = 12, dpi = 300)



### MAIN FIGURES 3 AND 4 ###
polyDFE = fread("../results/data_frames/polyDFE_best_fits.txt", fill = T)
dadi = fread("../results/data_frames/dadi_best_fits.txt", fill = T)

### Shape parameter ###
dadi_b = dplyr::select(dadi, params, shape)
dadi_b$method = "dadi"
dadi_b$population = "Ancestral"

polyDFE_b = dplyr::select(polyDFE, b, Population, Simulation)
polyDFE_b$method = "polyDFE"
#polyDFE_b$Simulation = gsub("X", "", polyDFE_b$Simulation)
#polyDFE_b$Simulation = gsub("\\.", "-", polyDFE_b$Simulation)
polyDFE_b$Simulation <- factor(polyDFE_b$Simulation, levels = c('X0.0.0', 'X0.0.0.05', 'X0.0.0.25', 'X0.01.0.0', 'X0.01.0.0.05', 'X0.01.0.0.25',  'X0.2.0', 'X0.2.0.05', 'X0.2.0.25', 'X0.01.2.0', 'X0.01.2.0.05', 'X0.01.2.0.25', 'X0.20.0', 'X0.20.0.05', 'X0.20.0.25', 'X0.01.20.0', 'X0.01.20.0.05', 'X0.01.20.0.25', 'X0.200.0', 'X0.200.0.05', 'X0.200.0.25', 'X0.01.200.0', 'X0.01.200.0.05', 'X0.01.200.0.25'), labels = c('0-0-0','0-0-0.05', '0-0-0.25','0.01-0-0','0.01-0-0.05', '0.01-0-0.25', '0-2-0','0-2-0.05', '0-2-0.25', '0.01-2-0','0.01-2-0.05', '0.01-2-0.25', '0-20-0','0-20-0.05', '0-20-0.25', '0.01-20-0','0.01-20-0.05', '0.01-20-0.25', '0-200-0','0-200-0.05', '0-200-0.25', '0.01-200-0','0.01-200-0.05', '0.01-200-0.25'), ordered = TRUE)

colnames(polyDFE_b) <- c("shape", "population", "params", "method")

b = rbind(dadi_b, polyDFE_b)
b$params2 = b$params
b = b %>%
  separate(params2, c("mig", "possel", "pchange"), "-")
b$params = paste0("m = ", b$mig, ", Sb = ", as.numeric(b$possel)/2, " and pc = ", b$pchange)

myplot = ggplot(b, aes(x=shape, fill=population)) +
  scale_fill_manual(values = c("black","#1b9e77", "#d95f02")) +
  geom_density(adjust=2, alpha=0.5) +
  facet_wrap(~ params, ncol = 3) +
  scale_x_log10() +
  geom_vline(xintercept = 0.3, linetype="dashed") +
  xlab("") + ylab("") +
  theme_bw() +
  theme(text = element_text(size=12),legend.position="none") 
ggsave(myplot, file="../plots/Figure3.png", width = 8.5, height = 12, dpi = 300)

### Mean parameter ###
dadi_Sd = dplyr::select(dadi, params, shape, scale, theta_nsyn)
dadi_Sd$Sd = 40*((dadi_Sd$shape*dadi_Sd$scale)/dadi_Sd$theta_nsyn)
dadi_Sd = dplyr::select(dadi_Sd, params, Sd)
dadi_Sd$method = "dadi"
dadi_Sd$population = "Ancestral"

polyDFE_Sd = dplyr::select(polyDFE, Sd, Population, Simulation)
polyDFE_Sd$method = "polyDFE"
polyDFE_Sd$Simulation <- factor(polyDFE_Sd$Simulation, levels = c('X0.0.0', 'X0.0.0.05', 'X0.0.0.25', 'X0.01.0.0', 'X0.01.0.0.05', 'X0.01.0.0.25',  'X0.2.0', 'X0.2.0.05', 'X0.2.0.25', 'X0.01.2.0', 'X0.01.2.0.05', 'X0.01.2.0.25', 'X0.20.0', 'X0.20.0.05', 'X0.20.0.25', 'X0.01.20.0', 'X0.01.20.0.05', 'X0.01.20.0.25', 'X0.200.0', 'X0.200.0.05', 'X0.200.0.25', 'X0.01.200.0', 'X0.01.200.0.05', 'X0.01.200.0.25'), labels = c('0-0-0','0-0-0.05', '0-0-0.25','0.01-0-0','0.01-0-0.05', '0.01-0-0.25', '0-2-0','0-2-0.05', '0-2-0.25', '0.01-2-0','0.01-2-0.05', '0.01-2-0.25', '0-20-0','0-20-0.05', '0-20-0.25', '0.01-20-0','0.01-20-0.05', '0.01-20-0.25', '0-200-0','0-200-0.05', '0-200-0.25', '0.01-200-0','0.01-200-0.05', '0.01-200-0.25'), ordered = TRUE)

colnames(polyDFE_Sd) <- c("Sd", "population", "params", "method")

Sd = rbind(dadi_Sd, polyDFE_Sd)
Sd$params2 = Sd$params
Sd = Sd %>%
  separate(params2, c("mig", "possel", "pchange"), "-")
Sd$params = paste0("m = ", Sd$mig, ", Sb = ", as.numeric(Sd$possel)/2, " and pc = ", Sd$pchange)


myplot = ggplot(Sd, aes(x=abs(Sd), fill=population)) +
  scale_fill_manual(values = c("black","#1b9e77", "#d95f02")) +
  geom_density(adjust=2, alpha=0.5) +
  facet_wrap(~ params, ncol = 3) +
  geom_vline(xintercept = 0.01, linetype="dashed") +
  xlab("") + ylab("") +
  theme_bw() +
  theme(text = element_text(size=12),legend.position="none") 
ggsave(myplot, file="../plots/Figure4.png", width = 8.5, height = 12, dpi = 300)


### MAIN FIGURE 5 ###
dadi2 <- dadi %>% 
  mutate(scenario = case_when(true_possel == 0 ~   "0% beneficial mutations",
                              true_possel == 2 ~   "10% nearly neutral beneficial mutations",
                              true_possel == 20 ~  "1% weakly beneficial mutations",
                              true_possel == 200 ~ "0.1% strongly beneficial mutations"))
dadi2$scenario <- factor(dadi2$scenario, levels = c("0% beneficial mutations", "10% nearly neutral beneficial mutations", "1% weakly beneficial mutations", "0.1% strongly beneficial mutations"), ordered = TRUE)

dummy = fread("dummy_fig5.txt")
dummy$true_pchange = as.factor(dummy$true_pchange)
dummy$mig = as.factor(dummy$mig)
dummy$value = "True"

dt_fig5 = dplyr::select(dadi2, scenario, true_pchange, mig, pchange)
dt_fig5$value = "Inferred"
dt_fig5 = rbind(dt_fig5, dummy)

ggplot() +
  geom_violin( data=subset(dt_fig5, value == "Inferred"), aes(x = true_pchange, y = as.numeric(pchange), fill =mig)) +
  geom_boxplot(data=subset(dt_fig5, value == "True"),     aes(x = true_pchange, y = as.numeric(pchange), color="black")) +
  theme_bw() +
  scale_fill_brewer(palette = "Greens") +
  facet_wrap(~scenario, ncol=1) +
  ylab("Inferred pc") + xlab("Simulated pc") + 
  #theme_classic() +
  theme(text = element_text(size=25), legend.position = "none")
ggsave(file="../plots/Figure5.png", width=10, height=12, dpi=300)


### SUPP FIGURE 4 ###

### pb parameter ###
dadi_pb = dplyr::select(dadi, params, pb_wild, pb_domesticated)
dadi_pb$method = "dadi"
# change from wide to long
dadi_pb <- gather(dadi_pb, population, pb, pb_wild:pb_domesticated, factor_key=TRUE)
dadi_pb$population = gsub("pb_wild", "Wild", dadi_pb$population)
dadi_pb$population = gsub("pb_domesticated", "Domesticated", dadi_pb$population)

polyDFE_pb = dplyr::select(polyDFE, pb, Population, Simulation)
polyDFE_pb$method = "polyDFE"
polyDFE_pb$Simulation <- factor(polyDFE_pb$Simulation, levels = c('X0.0.0', 'X0.0.0.05', 'X0.0.0.25', 'X0.01.0.0', 'X0.01.0.0.05', 'X0.01.0.0.25',  'X0.2.0', 'X0.2.0.05', 'X0.2.0.25', 'X0.01.2.0', 'X0.01.2.0.05', 'X0.01.2.0.25', 'X0.20.0', 'X0.20.0.05', 'X0.20.0.25', 'X0.01.20.0', 'X0.01.20.0.05', 'X0.01.20.0.25', 'X0.200.0', 'X0.200.0.05', 'X0.200.0.25', 'X0.01.200.0', 'X0.01.200.0.05', 'X0.01.200.0.25'), labels = c('0-0-0','0-0-0.05', '0-0-0.25','0.01-0-0','0.01-0-0.05', '0.01-0-0.25', '0-2-0','0-2-0.05', '0-2-0.25', '0.01-2-0','0.01-2-0.05', '0.01-2-0.25', '0-20-0','0-20-0.05', '0-20-0.25', '0.01-20-0','0.01-20-0.05', '0.01-20-0.25', '0-200-0','0-200-0.05', '0-200-0.25', '0.01-200-0','0.01-200-0.05', '0.01-200-0.25'), ordered = TRUE)

colnames(polyDFE_pb) <- c("pb", "population", "params", "method")

pb = rbind(dadi_pb, polyDFE_pb)
pb$params2 = pb$params
pb = pb %>%
  separate(params2, c("mig", "possel", "pchange"), "-")
pb$params = paste0("m = ", pb$mig, ", Sb = ", as.numeric(pb$possel)/2, " and pc = ", pb$pchange)

A = ggplot(pb, aes(x=pb+0.0001, y=population, fill=population, linetype=method)) +
  ggtitle("A") +
  scale_fill_manual(values = c("#1b9e77", "#d95f02")) +
  #geom_density(adjust=2, alpha=0.5) +
  geom_violin() +
  facet_wrap(~ params, ncol = 3) +
  scale_x_log10() +
  #scale_x_continuous(trans=scales::pseudo_log_trans(base = 10)) +
  #geom_vline(xintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0.001, linetype="dashed") +
  geom_vline(xintercept = 0.01, linetype="dashed") +  
  geom_vline(xintercept = 0.1, linetype="dashed") +  
  xlab("") + ylab("") +
  theme_bw() +
  theme(text = element_text(size=12), axis.text.y=element_blank(), axis.title.y=element_blank(),legend.position="left")
A


polyDFE_Sb = dplyr::select(polyDFE, Sb, Population, Simulation)
polyDFE_Sb$method = "polyDFE"
polyDFE_Sb$Simulation <- factor(polyDFE_Sb$Simulation, levels = c('X0.0.0', 'X0.0.0.05', 'X0.0.0.25', 'X0.01.0.0', 'X0.01.0.0.05', 'X0.01.0.0.25',  'X0.2.0', 'X0.2.0.05', 'X0.2.0.25', 'X0.01.2.0', 'X0.01.2.0.05', 'X0.01.2.0.25', 'X0.20.0', 'X0.20.0.05', 'X0.20.0.25', 'X0.01.20.0', 'X0.01.20.0.05', 'X0.01.20.0.25', 'X0.200.0', 'X0.200.0.05', 'X0.200.0.25', 'X0.01.200.0', 'X0.01.200.0.05', 'X0.01.200.0.25'), labels = c('0-0-0','0-0-0.05', '0-0-0.25','0.01-0-0','0.01-0-0.05', '0.01-0-0.25', '0-2-0','0-2-0.05', '0-2-0.25', '0.01-2-0','0.01-2-0.05', '0.01-2-0.25', '0-20-0','0-20-0.05', '0-20-0.25', '0.01-20-0','0.01-20-0.05', '0.01-20-0.25', '0-200-0','0-200-0.05', '0-200-0.25', '0.01-200-0','0.01-200-0.05', '0.01-200-0.25'), ordered = TRUE)
polyDFE_Sb$params2 = polyDFE_Sb$Simulation
polyDFE_Sb = polyDFE_Sb %>%
  separate(params2, c("mig", "possel", "pchange"), "-")

polyDFE_Sb$params = paste0("m = ", polyDFE_Sb$mig, ", Sb = ", as.numeric(polyDFE_Sb$possel)/2, " and pc = ", polyDFE_Sb$pchange)

B = ggplot(polyDFE_Sb, aes(x=Sb, y=Population, fill=Population)) +
  ggtitle("B") +
  scale_fill_manual(values = c("#1b9e77", "#d95f02")) +
  #geom_density(adjust=2, alpha=0.5) +
  geom_violin(linetype = "dashed") +
  facet_wrap(~ params, ncol = 3) +
  scale_x_log10() +
  geom_vline(aes(xintercept=(1/10000)),
             color="black", linetype="dotted", size=1) +
  geom_vline(aes(xintercept=(10/10000)),
             color="black", linetype="dotted", size=1) +
  geom_vline(aes(xintercept=(100/10000)),
             color="black", linetype="dotted", size=1) +
  xlab("") + ylab("") +
  theme_bw() +
  theme(text = element_text(size=12), axis.text.y=element_blank(), axis.title.y=element_blank(),legend.position="none")
B

myplot = grid.arrange(A, B, ncol=2)
ggsave(myplot, file="../plots/SuppFigure4.png", width = 17, height = 12, dpi = 300)


dadi_Sb = dplyr::select(dadi, params, Sb, true_possel, pb_wild, pb_domesticated, params)
dadi_pb$method = "dadi"

dadi_Sb$Sb_wild <- ifelse(dadi_Sb$pb_wild == "0", # condition
                          0,    # what if condition is TRUE
                          dadi_Sb$Sb          # what if condition is FALSE
)

dadi_Sb$Sb_domestic <- ifelse(dadi_Sb$pb_domesticated == "0", # condition
                          0,    # what if condition is TRUE
                          dadi_Sb$Sb          # what if condition is FALSE
)

dadi_Sb = dadi_Sb %>%
  separate(params, c("mig", "possel", "pchange"), "-")


wb <- createWorkbook()
for (MIGRATION in c(0, 0.01))
{
  for (POSSEL in c(0, 2, 20, 200))
  {
    for (CHANGE in c(0, 0.05, 0.25))
    {
      print(paste(MIGRATION, POSSEL, CHANGE))
      scenario = subset(dadi_Sb, mig == MIGRATION & possel == POSSEL & pchange == CHANGE)
      
      sheet_name = paste0("D, m=", MIGRATION, ", Sb=", POSSEL/2, " and pc=", CHANGE)
      addWorksheet(wb,sheet_name)
      
      result = table(list(scenario$Sb_domestic, scenario$true_possel/2))
      
      writeData(wb, sheet_name, result, rowNames=TRUE, colNames = TRUE)
      
      
      sheet_name = paste0("W, m=", MIGRATION, ", Sb=", POSSEL/2, " and pc=", CHANGE)
      addWorksheet(wb,sheet_name)
      
      result = table(list(scenario$Sb_wild, scenario$true_possel/2))
      
      writeData(wb, sheet_name, result, rowNames=TRUE, colNames = TRUE)
      
      
    }
  }
}

saveWorkbook(wb, "../plots/Supp_Table3.xlsx", overwrite = TRUE)


### Supplementary Table 2 ###
data = read.delim("../results/data_frames/polyDFE_lrt_replicate.txt", sep = " ", header = F)
data = unique(data)
colnames(data) <- c("pvalue_model1_vs_model10", "pvalue_model1_vs_model2", "pvalue_model2_vs_model20", "pvalue_model10_vs_model20", "pvalue_model20_vs_model30", "mig", "possel", "pchange", "it")
data$params = paste0(data$mig, "-", data$possel, "-", data$pchange)
table(data$params) #number of replicates per scenario

data$params <- factor(data$params, levels = c('0-0-0','0-0-0.05', '0-0-0.25','0.01-0-0','0.01-0-0.05', '0.01-0-0.25', '0-2-0','0-2-0.05', '0-2-0.25', '0.01-2-0','0.01-2-0.05', '0.01-2-0.25', '0-20-0','0-20-0.05', '0-20-0.25', '0.01-20-0','0.01-20-0.05', '0.01-20-0.25', '0-200-0','0-200-0.05', '0-200-0.25', '0.01-200-0','0.01-200-0.05', '0.01-200-0.25'), labels = c("m=0, Sb=0, pchange=0%", "Mig=0%, Sb=0, pb=0%, Change=5%", "Mig=0%, Sb=0, pb=0%, Change=25%", "Mig=1%, Sb=0, pb=0%, Change=0%", "Mig=1%, Sb=0, pb=0%, Change=5%", "Mig=1%, Sb=0, pb=0%, Change=25%", "Mig=0%, Sb=1, pb=10%, Change=0%", "Mig=0%, Sb=1, pb=10%, Change=5%", "Mig=0%, Sb=1, pb=10%, Change=25%", "Mig=1%, Sb=1, pb=10%, Change=0%", "Mig=1%, Sb=1, pb=10%, Change=5%", "Mig=1%, Sb=1, pb=10%, Change=25%", "Mig=0%, Sb=10, pb=1%, Change=0%", "Mig=0%, Sb=10, pb=1%, Change=5%", "Mig=0%, Sb=10, pb=1%, Change=25%", "Mig=1%, Sb=10, pb=1%, Change=0%", "Mig=1%, Sb=10, pb=1%, Change=5%", "Mig=1%, Sb=10, pb=1%, Change=25%", "Mig=0%, Sb=100, pb=0.1%, Change=0%", "Mig=0%, Sb=100, pb=0.1%, Change=5%", "Mig=0%, Sb=100, pb=0.1%, Change=25%", "Mig=1%, Sb=100, pb=0.1%, Change=0%", "Mig=1%, Sb=100, pb=0.1%, Change=5%", "Mig=1%, Sb=100, pb=0.1%, Change=25%"), ordered = TRUE)

median_pvalue = data.frame(t(rbind(tapply(data$pvalue_model1_vs_model10, data$params, median),
                                 tapply(data$pvalue_model2_vs_model20, data$params, median),
                                 tapply(data$pvalue_model1_vs_model2, data$params, median),
                                 tapply(data$pvalue_model10_vs_model20, data$params, median),
                                 tapply(data$pvalue_model20_vs_model30, data$params, median))))
colnames(median_pvalue) <- c("M1 vs M10", "M2 vs M20", "M1 vs M2", "M10 vs M20", "M20 vs M30")
median_pvalue$params = row.names(median_pvalue)

median_pvalue$params2 = median_pvalue$params
median_pvalue = median_pvalue %>%
  separate(params2, c("mig", "possel", "pchange"), "-")
median_pvalue$params = paste0("m = ", median_pvalue$mig, ", Sb = ", as.numeric(median_pvalue$possel)/2, " and pc = ", median_pvalue$pchange)

median_pvalue = dplyr::select(median_pvalue, "params", "M1 vs M10", "M2 vs M20", "M1 vs M2", "M10 vs M20", "M20 vs M30")

write.table(median_pvalue, file="../plots/Supplementary_Table2.xlsx", quote = F, sep = "\t", append = FALSE, row.names = F)


### Supplementary Figure 5 ###
param_list =  fread("parameters_list.txt", fill = T)
MEAN_SD = 0.01
POP = "BOTH"

expectation = data.frame()
for (i in seq(24))
{
  # i = 18
  ROW = slice(param_list, i)
  POSSEL = ROW$possel
  SCALING = ROW$scaling
  
  if (POSSEL == 0)
  {
    Sb = MEAN_SD/100
    pb = 0
  }
  
  
  if (POSSEL == 2)
  {
    Sb = MEAN_SD/100
    pb = 0.1
  }
  
  if (POSSEL == 20)
  {
    Sb = MEAN_SD/10
    pb = 0.01
  }
  
  if (POSSEL == 200)
  {
    Sb = MEAN_SD/1
    pb = 0.001
  }
  
  neu = pgamma(q=c(0.001), shape=0.3, scale=(MEAN_SD/0.3)*SCALING)
  strong_del = 1-pgamma(q=c(0.01), shape=0.3, scale=(MEAN_SD/0.3)*SCALING)
  weak_del = 1 - neu - strong_del
  neu2 = pexp(q=c(0),rate=1/(Sb*SCALING))
  ben = 1-neu2
  
  BENEFICIAL = pb*ben
  NEUTRAL = (1-pb)*neu
  WEAKLY_DEL = (1-pb)*weak_del
  STRONG_DEL = (1-pb)*strong_del
  
  BENEFICIAL+NEUTRAL+WEAKLY_DEL+STRONG_DEL
  
  exp_1 = c(paste0(ROW$mig,"-",ROW$possel,"-",ROW$pchange), POP, "< -0.01", STRONG_DEL, STRONG_DEL, STRONG_DEL, SCALING)
  exp_2 = c(paste0(ROW$mig,"-",ROW$possel,"-",ROW$pchange), POP, "(-0.01, -0.001]", WEAKLY_DEL, WEAKLY_DEL, WEAKLY_DEL, SCALING)
  exp_3 = c(paste0(ROW$mig,"-",ROW$possel,"-",ROW$pchange), POP, "(-0.001, 0]", NEUTRAL, NEUTRAL, NEUTRAL, SCALING)
  exp_4 = c(paste0(ROW$mig,"-",ROW$possel,"-",ROW$pchange), POP, "> 0", BENEFICIAL, BENEFICIAL, BENEFICIAL, SCALING)
  
  results = rbind(exp_1, exp_2, exp_3, exp_4)
  expectation = rbind(expectation, results)
}


colnames(expectation) <- c("scenario", "population", "coeff", "value", "CI5", "CI95", "scaling")
expectation$Value = "True"
expectation$value = as.numeric(expectation$value)
expectation$CI5 = as.numeric(expectation$CI5)
expectation$CI95 = as.numeric(expectation$CI95)
expectation$Population = expectation$population
expectation$coeff <- factor(expectation$coeff, levels = c("< -0.01", "(-0.01, -0.001]", "(-0.001, 0]", "> 0"), ordered = TRUE)

options(scipen=10000)
ggplot(expectation, aes(y=as.numeric(value), x=as.factor(coeff), shape=Population, color=Value)) +
  geom_point(position=position_dodge(width = .5),size=3) +
  # geom_errorbar(aes(ymin = CI5, ymax = CI95)) + 
  geom_pointrange(aes(ymin=CI5, ymax=CI95), position=position_dodge(width = .5)) +
  scale_color_grey() +
  xlab("s") +
  ylab("Fraction") +
  ylim(0,0.5) +
  #scale_y_log10(breaks=c(0.0001, 0.001, 0.01, 0.1, 0.5, 1)) +
  theme_bw() +  
  facet_wrap(~scenario, ncol = 3) +
  theme(text = element_text(size=20), axis.title.y=element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "bottom")

tapply(expectation$value, list(paste0(expectation$scenario,"-",expectation$Value), expectation$population), sum, na.rm=T) 


# dadi #

dt_int = fread("../results/data_frames/dadi_best_fits.txt", fill = T)

dt_int$sd = 2*20*((dt_int$shape*dt_int$scale)/dt_int$theta_nsyn)
dt_int$neu = pgamma(q=c(0.001), shape=dt_int$shape, scale=(dt_int$sd/dt_int$shape))
dt_int$strong_del = 1-pgamma(q=c(0.01), shape=dt_int$shape, scale=(dt_int$sd/dt_int$shape))
dt_int$weak_del = 1 - dt_int$neu - dt_int$strong_del
dt_int$ben_dom = dt_int$pb_domesticated*dt_int$pchange+(1-dt_int$pchange)*dt_int$pb_wild
dt_int$ben_wild = dt_int$pb_wild

dt_int$BENEFICIAL_Domesticated = dt_int$ben_dom
dt_int$NEUTRAL_Domesticated = (1-dt_int$ben_dom)*dt_int$neu
dt_int$WEAKLY_Domesticated = (1-dt_int$ben_dom)*dt_int$weak_del
dt_int$STRONG_Domesticated = (1-dt_int$ben_dom)*dt_int$strong_del

dt_int$BENEFICIAL_Wild = dt_int$ben_wild
dt_int$NEUTRAL_Wild = (1-dt_int$ben_wild)*dt_int$neu
dt_int$WEAKLY_Wild = (1-dt_int$ben_wild)*dt_int$weak_del
dt_int$STRONG_Wild = (1-dt_int$ben_wild)*dt_int$strong_del

dt_int = dplyr::select(dt_int, params, BENEFICIAL_Domesticated:STRONG_Wild)
dt_int_long <- gather(dt_int, Coeff_Population, Value, BENEFICIAL_Domesticated:STRONG_Wild, factor_key=TRUE)
dt_int_long = dt_int_long %>%
  separate(Coeff_Population, c("coeff", "Population"), "_")

value =  tapply(dt_int_long$Value, paste0(dt_int_long$params,"_",dt_int_long$Population,"_", dt_int_long$coeff),  function(x) quantile(unlist(x), probs=c(0.5) ) ) 
CI5 =  tapply(dt_int_long$Value, paste0(dt_int_long$params,"_",dt_int_long$Population,"_", dt_int_long$coeff),  function(x) quantile(unlist(x), probs=c(0.025) ) ) 
CI95 =  tapply(dt_int_long$Value, paste0(dt_int_long$params,"_",dt_int_long$Population,"_", dt_int_long$coeff),  function(x) quantile(unlist(x), probs=c(0.975) ) ) 

dt_int = data.frame(cbind(value, CI5, CI95))
dt_int$info = row.names(dt_int)

dt_int = dt_int %>%
  separate(info, c("scenario", "Population", "coeff"), "_")

dt_int$coeff = gsub("BENEFICIAL", "> 0", dt_int$coeff)
dt_int$coeff = gsub("NEUTRAL", "(-0.001, 0]", dt_int$coeff)
dt_int$coeff = gsub("WEAKLY", "(-0.01, -0.001]", dt_int$coeff)
dt_int$coeff = gsub("STRONG", "< -0.01", dt_int$coeff)
#dt_int$scenario <- factor(dt_int$scenario, levels = c('0-2-0','0-2-0.05', '0-2-0.25', '0.01-2-0','0.01-2-0.05', '0.01-2-0.25', '0-20-0','0-20-0.05', '0-20-0.25', '0.01-20-0','0.01-20-0.05', '0.01-20-0.25', '0-200-0','0-200-0.05', '0-200-0.25', '0.01-200-0','0.01-200-0.05', '0.01-200-0.25'), labels = c("Mig=0%, Sb=1, pb=10%, Change=0%", "Mig=0%, Sb=1, pb=10%, Change=5%", "Mig=0%, Sb=1, pb=10%, Change=25%", "Mig=1%, Sb=1, pb=10%, Change=0%", "Mig=1%, Sb=1, pb=10%, Change=5%", "Mig=1%, Sb=1, pb=10%, Change=25%", "Mig=0%, Sb=10, pb=1%, Change=0%", "Mig=0%, Sb=10, pb=1%, Change=5%", "Mig=0%, Sb=10, pb=1%, Change=25%", "Mig=1%, Sb=10, pb=1%, Change=0%", "Mig=1%, Sb=10, pb=1%, Change=5%", "Mig=1%, Sb=10, pb=1%, Change=25%", "Mig=0%, Sb=100, pb=0.1%, Change=0%", "Mig=0%, Sb=100, pb=0.1%, Change=5%", "Mig=0%, Sb=100, pb=0.1%, Change=25%", "Mig=1%, Sb=100, pb=0.1%, Change=0%", "Mig=1%, Sb=100, pb=0.1%, Change=5%", "Mig=1%, Sb=100, pb=0.1%, Change=25%"), ordered = TRUE)

dt_int$scaling = 1
dt_int$Value = "Inferred"
dt_int$population = dt_int$Population

expectation = rbind(expectation, dt_int)
expectation$value = as.numeric(expectation$value)
expectation$CI5 = as.numeric(expectation$CI5)
expectation$CI95 = as.numeric(expectation$CI95)
expectation$coeff <- factor(expectation$coeff, levels = c("< -0.01", "(-0.01, -0.001]", "(-0.001, 0]", "> 0"), ordered = TRUE)

write.table(expectation, file="../results/data_frames/discretized_DFE_dadi_vs_true.txt", col.names = T, row.names = F, sep = "\t", quote = F, append = F)


# polyDFE #
dt_int = fread("../results/data_frames/polyDFE_best_fits.txt", fill = T)
dt_int$Simulation <- factor(dt_int$Simulation, levels = c('X0.0.0', 'X0.0.0.05', 'X0.0.0.25', 'X0.01.0.0', 'X0.01.0.0.05', 'X0.01.0.0.25',  'X0.2.0', 'X0.2.0.05', 'X0.2.0.25', 'X0.01.2.0', 'X0.01.2.0.05', 'X0.01.2.0.25', 'X0.20.0', 'X0.20.0.05', 'X0.20.0.25', 'X0.01.20.0', 'X0.01.20.0.05', 'X0.01.20.0.25', 'X0.200.0', 'X0.200.0.05', 'X0.200.0.25', 'X0.01.200.0', 'X0.01.200.0.05', 'X0.01.200.0.25'), labels = c('0-0-0','0-0-0.05', '0-0-0.25','0.01-0-0','0.01-0-0.05', '0.01-0-0.25', '0-2-0','0-2-0.05', '0-2-0.25', '0.01-2-0','0.01-2-0.05', '0.01-2-0.25', '0-20-0','0-20-0.05', '0-20-0.25', '0.01-20-0','0.01-20-0.05', '0.01-20-0.25', '0-200-0','0-200-0.05', '0-200-0.25', '0.01-200-0','0.01-200-0.05', '0.01-200-0.25'), ordered = TRUE)

colnames(dt_int) <- c("params", "Replicate", "Population", "sd", "shape", "sb", "pb")
dt_int$sd = -dt_int$sd

#dt_int$sd = 2*20*((dt_int$shape*dt_int$scale)/dt_int$theta_nsyn)
dt_int$neu = pgamma(q=c(0.001), shape=dt_int$shape, scale=(dt_int$sd/dt_int$shape))
dt_int$strong_del = 1-pgamma(q=c(0.01), shape=dt_int$shape, scale=(dt_int$sd/dt_int$shape))
dt_int$weak_del = 1 - dt_int$neu - dt_int$strong_del
dt_int$ben = dt_int$pb

dt_int$BENEFICIAL = dt_int$ben
dt_int$NEUTRAL = (1-dt_int$ben)*dt_int$neu
dt_int$WEAKLY = (1-dt_int$ben)*dt_int$weak_del
dt_int$STRONG = (1-dt_int$ben)*dt_int$strong_del


dt_int = dplyr::select(dt_int, params, Population, BENEFICIAL:STRONG)
dt_int_long <- gather(dt_int, coeff, Value, BENEFICIAL:STRONG, factor_key=TRUE)

value =  tapply(dt_int_long$Value, paste0(dt_int_long$params,"_",dt_int_long$Population,"_", dt_int_long$coeff),  function(x) quantile(unlist(x), probs=c(0.5) ) ) 
CI5 =  tapply(dt_int_long$Value, paste0(dt_int_long$params,"_",dt_int_long$Population,"_", dt_int_long$coeff),  function(x) quantile(unlist(x), probs=c(0.025) ) ) 
CI95 =  tapply(dt_int_long$Value, paste0(dt_int_long$params,"_",dt_int_long$Population,"_", dt_int_long$coeff),  function(x) quantile(unlist(x), probs=c(0.975) ) ) 


dt_int = data.frame(cbind(value, CI5, CI95))
dt_int$info = row.names(dt_int)

dt_int = dt_int %>%
  separate(info, c("scenario", "Population", "coeff"), "_")

dt_int$coeff = gsub("BENEFICIAL", "> 0", dt_int$coeff)
dt_int$coeff = gsub("NEUTRAL", "(-0.001, 0]", dt_int$coeff)
dt_int$coeff = gsub("WEAKLY", "(-0.01, -0.001]", dt_int$coeff)
dt_int$coeff = gsub("STRONG", "< -0.01", dt_int$coeff)

dt_int$scaling = 1
dt_int$Value = "Inferred"
dt_int$population = dt_int$Population
dt_int$Method = "polyDFE"

write.table(dt_int, file="../results/data_frames/discretized_DFE_polyDFE.txt", col.names = T, row.names = F, sep = "\t", quote = F, append = F)


# Merge discretized DFEs datasets #

dt_int_dadi = fread(file = "../results/data_frames/discretized_DFE_dadi_vs_true.txt", header = T)
dt_int_polyDFE = fread(file = "../results/data_frames/discretized_DFE_polyDFE.txt", header = T)
dt_int_dadi$Method = "dadi"

expectation = rbind(dt_int_dadi, dt_int_polyDFE)

expectation$Population = gsub("BOTH", "Ancestral", expectation$Population)
expectation$Method[expectation$population == "BOTH" & expectation$Method == 'dadi'] <- "True Value"
expectation$coeff <- factor(expectation$coeff, levels = c("< -0.01", "(-0.01, -0.001]", "(-0.001, 0]", "> 0"), ordered = TRUE)

expectation$params2 = expectation$scenario
expectation = expectation %>%
  separate(params2, c("mig", "possel", "pchange"), "-")
expectation$scenario = paste0("m = ", expectation$mig, ", Sb = ", as.numeric(expectation$possel)/2, " and pc = ", expectation$pchange)

options(scipen=10000)
dfe_plot = ggplot(expectation, aes(y=as.numeric(value), x=as.factor(coeff), shape=Method, color=Population)) +
  geom_point(position=position_dodge(width = .5),size=3) +
  geom_pointrange(aes(ymin=CI5, ymax=CI95), position=position_dodge(width = .5)) +
  scale_color_manual(values = c("black","#1b9e77", "#d95f02")) +
  xlab("s") +
  ylab("Fraction") +
  theme_bw() +  
  facet_wrap(~scenario, ncol = 3) +
  theme(text = element_text(size=20), axis.title.y=element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "bottom")
ggsave(plot = dfe_plot, "../plots/Supp_Figure5.png", width=20, height=15, dpi="retina")

