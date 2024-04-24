library("ggpubr")
library("readxl")

# hard-coded data-path
data_path = "C:/Users/jxe8989/OneDrive - University of Texas at Arlington/Research/diss/synth_data/stats - Copy/"

# SQUARE PROFILE
table = read_excel(paste(data_path, 
                         "EvalMetrics_protocol04_30_square.xlsx", sep = "") )
title_text = "Square profile"

# GAUSS PROFILE
table = read_excel(paste(data_path, 
                         "EvalMetrics_protocol04_30_gauss.xlsx", sep = "") )
title_text = "Gaussian profile"

# EXPONENTIAL PROFILE
table = read_excel(paste(data_path, 
                         "EvalMetrics_protocol04_30_exp.xlsx", sep = "") )
title_text = "Exponential profile"

################################################################################
# FORMATTING

# format (hard-coded for now)
table$SNR = factor(table$SNR, 
                   levels = c( "Inf", "30", "20", "10", "0", "-10" ) 
)
table$HalfMax = table$HalfMax/100

################################################################################
# WHICH LEVEL OF NOISE IS OK ?

# noiseless
p <- ggviolin(table[ table$SNR=="Inf", ], 
              x = "Solver",
              y = "DLE1",
              add = "mean_sd") +
  ylab("Distance Localization Eror [mm]") +
  ggtitle("Noiseless case, SNR = Inf dB") +
  grids()
p

# SNR = 30
p <- ggviolin(table[ table$SNR=="30", ], 
              x = "Solver",
              y = "DLE1",
              add = "mean_sd") +
  ylab("Distance Localization Eror [mm]") +
  ggtitle("Noiseless case, SNR = 30 dB") +
  grids()
p

# SNR = 20
p <- ggviolin(table[ table$SNR=="20", ], 
              x = "Solver",
              y = "DLE1",
              add = "mean_sd") +
  ylab("Distance Localization Eror [mm]") +
  ggtitle("Noiseless case, SNR = 20 dB") +
  grids()
p

# SNR = 10
p <- ggviolin(table[ table$SNR=="10", ], 
              x = "Solver",
              y = "DLE1",
              add = "mean_sd") +
  ylab("Distance Localization Eror [mm]") +
  ggtitle("Noiseless case, SNR = 10 dB") +
  grids()
p

# SNR = 0
p <- ggviolin(table[ table$SNR=="0", ], 
              x = "Solver",
              y = "DLE1",
              add = "mean_sd") +
  ylab("Distance Localization Eror [mm]") +
  ggtitle("Noiseless case, SNR = 0 dB") +
  grids()
p

# SNR = -10
p <- ggviolin(table[ table$SNR=="-10", ], 
              x = "Solver",
              y = "DLE1",
              add = "mean_sd") +
  ylab("Distance Localization Eror [mm]") +
  ggtitle("Noiseless case, SNR = -10 dB") +
  grids()
p

################################################################################
# WILL KEEP SNR UP TO 10

# Distance Localization Error
p <- ggviolin(table[ table$SNR %in% c("Inf","30","20","10"), ], 
              x = "SNR",
              y = "DLE1",
              fill = "Solver",
              add = "mean_sd") +
  ylab("Distance Localization Eror [mm]") +
  ggtitle(title_text) +
  grids()
p

# Spatial Dispersion
p <- ggviolin(table[ table$SNR %in% c("Inf","30","20","10"), ], 
              x = "SNR",
              y = "SpaDis2",
              fill = "Solver",
              add = "mean_sd") +
  ylab("Spatial Dispersion [mm]") +
  ggtitle(title_text) +
  grids() +
  geom_hline(yintercept=10*sqrt(5/pi)*((1/2)^(1/2)), 
             linetype="dashed", color = "blue")
p

# Spatial Dispersion, modified
p <- ggviolin(table[ table$SNR %in% c("Inf","30","20","10"), ], 
              x = "SNR",
              y = "SpaDis1",
              fill = "Solver",
              add = "mean_sd") +
  ylab("Spatial Dispersion (norm1) [mm]") +
  ggtitle(title_text) +
  grids() +
  geom_hline(yintercept=10*sqrt(5/pi)*((1/2)^(1/2)), 
             linetype="dashed", color = "blue")
p

# Area Under ROC
p <- ggviolin(table[ table$SNR %in% c("Inf","30","20","10"), ], 
              x = "SNR",
              y = "AUROC_loc_w",
              fill = "Solver",
              add = "mean_sd") +
  ylab("AUROC (local)") +
  ggtitle(title_text) +
  grids()
p

# Avg Precision
p <- ggviolin(table[ table$SNR %in% c("Inf","30","20","10"), ], 
              x = "SNR",
              y = "AP_loc_w",
              fill = "Solver",
              add = "mean_sd") +
  ylab("Average precision (local)") +
  ggtitle(title_text) +
  grids()
p

# Half-Max
p <- ggviolin(table[ table$SNR %in% c("Inf","30","20","10"), ], 
              x = "SNR",
              y = "HalfMax",
              fill = "Solver") +
  ylab("Half-Max Area [cm^2]") +
  ggtitle(title_text) +
  grids()
p

# RMSE
p <- ggviolin(table[ table$SNR %in% c("Inf","30","20","10"), ], 
              x = "SNR",
              y = "RMSE",
              fill = "Solver") +
  ylab("Relative Mean-Square Error") +
  ggtitle(title_text) +
  grids()
p

################################################################################
#
#
#
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

