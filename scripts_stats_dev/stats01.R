library("ggpubr")

library("readxl")

table = read_excel("C:/Users/jxe8989/My Drive/Research/diss/tests_wiiiii/SimProt02/EvalMetrics_SimProt02.xlsx")
#table = read_excel("D:/Gugl/Research/diss/tests_wiiiii/test03/EvalMetrics_test03.xlsx")

table$SNR = factor(table$SNR, 
                   levels = c( "Inf", "35", "30", "25", "20", "15", "10", "5" ) 
)

table$Kappa = factor(table$Kappa)

table_regions = table[ table$Solver == "RegionPrior", ]

p <- ggviolin(table_regions, 
              x = "SNR", y = "LocalizationError",
              add = "mean_sd") +
  xlab("SNR [dB]") +
  ylab("Localization Eror [mm]") +
  ggtitle("Proposed Method") +
  grids()
p

p <- ggviolin(table_regions, 
              x = "Kappa", y = "Depth",
              add = "mean_sd") +
  xlab("Spread [mm]") +
  ylab("Depth [mm]") +
  grids()
p

p <- ggviolin(table_regions, 
              x = "Kappa", y = "LocalizationError",
              add = "mean_sd") +
  xlab("Spread [mm]") +
  ylab("Localizatio Error [mm]") +
  grids()
p

table_regions_20 = table_regions[ table_regions$SNR == 10, ]

p <- ggscatter(table_regions_20, 
              x = "Depth", y = "LocalizationError",
              color="Kappa") +
  xlab("Depth [mm]") +
  ylab("Localization Error [mm]") +
  grids()
p


##############

p <- ggviolin(table_regions, 
               x = "Kappa", y = "LocalizationError",
               color="SNR",
              add = "mean_sd") +
  xlab("Spread [mm]") +
  ylab("Localizatio Error [mm]") +
  ggtitle("Proposed Method") +
  grids()
p

table_sloreta = table[ table$Solver == "sLORETA", ]

p <- ggviolin(table_sloreta, 
              x = "Kappa", y = "LocalizationError",
              color="SNR",
              add = "mean_sd") +
  xlab("Spread [mm]") +
  ylab("Localizatio Error [mm]") +
  ggtitle("Proposed Method") +
  grids()
p

#p <- ggboxplot(table, 
#               x = "SNR", y = "LocalizationError",
#               fill = "Solver",
#)+ grids(linetype = "solid")
#
#p

##
ggboxplot(table, 
          x = "SNR", y = "LocalizationError",
          fill = "Solver",
)+ grids(linetype = "solid")

##
ggboxplot(table, 
          x = "SNR", y = "HalfMax",
          fill = "Solver",
)+ grids(linetype = "solid")

##
ggboxplot(table, 
          x = "SNR", y = "Algorithm Time",
          fill = "Solver",
)+ grids(linetype = "solid") +
  scale_y_continuous(trans='log10')


##
ggboxplot(table, 
          x = "SNR", y = "Parameter Tuning Time",
          fill = "Solver",
)+ grids(linetype = "solid") +
  scale_y_continuous(trans='log10')

