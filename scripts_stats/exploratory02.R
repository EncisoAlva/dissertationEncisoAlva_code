library("ggpubr")
library("readxl")
library("tidyverse")
library("plyr")
library("latex2exp")

# dataset to analyze
tagname = 'EvalMetrics_Protocol04_30'
#tagname = 'EvalMetrics_protocol04_vol5k_pig'

# current directory
script_dir = getwd()

# data path
setwd('..')
data_path = paste0( getwd(), '/stats - Copy/' )
setwd( script_dir )

################################################################################

# ALL PROFILES COMBINED
table = read_excel( paste(data_path, tagname, "_ALL.xlsx", sep = "") )
title_text = "Polynomial profile"

# SQUARE PROFILE
table = read_excel( paste(data_path, tagname, "_square.xlsx", sep = "") )
title_text = "Square profile"

# GAUSS PROFILE
table = read_excel( paste(data_path, tagname, "_gauss.xlsx", sep = "") )
title_text = "Gaussian profile"

# EXPONENTIAL PROFILE
table = read_excel( paste(data_path, tagname, "_exp.xlsx", sep = "") )
title_text = "Exponential profile"

# POLYNOMIAL PROFILE
table = read_excel( paste(data_path, tagname, "_circ.xlsx", sep = "") )
title_text = "Polynomial profile"

################################################################################
# FORMATTING

# format (hard-coded for now)
table$SNR = factor(table$SNR, 
                   levels = c( "Inf", "30", "20", "10", "0", "-10" ) 
)
table$HalfMax = table$HalfMax/1000

table$Solver = factor(table$Solver, 
                   levels = c( "Tikhonov", "MSP", "sLORETA", "zSISSY", "SingleRegionPrior") 
)
table$Solver = revalue(table$Solver, 
                      c( "Tikhonov"= "WMNE", "zSISSY"="SISSY", 
                         "SingleRegionPrior" = "Proposed") 
)

################################################################################

table %>%
  #drop_na() %>%
  filter( SNR=="30" ) %>%
  ggviolin(x = "Solver",
           y = "DLE1",
           add = "mean_sd") +
  ylab("Distance Localization Eror [mm]") +
  #ggtitle("Noiseless case, SNR = Inf dB") +
  grids() +
  theme_bw()

################################################################################

setwd(paste0( script_dir, "/img/" ))

table %>%
  #drop_na() %>%
  filter( SNR=="30" ) %>% 
  filter( Solver != "Proposed" ) %>%
  ggboxplot(x = "Solver",
           y = "DLE1") +
  ylab("Distance Localization Eror [mm]") +
  #ggtitle("Noiseless case, SNR = Inf dB") +
  #labs(subtitle = ("SNR = 30 [dB]" ) ) +
  #labs(subtitle = expression( paste("SNR =", infinity, "[dB]") ) ) +
  grids() +
  theme_bw()
ggsave( paste0("plot_",tagname, "DLE", ".pdf"),
        width = 4.5, height = 3, units = "in")

table %>%
  #drop_na() %>%
  filter( SNR=="30" ) %>%
  filter( Solver != "Proposed" ) %>%
  ggboxplot(x = "Solver",
            y = "SpaDis2") +
  ylab("Spatial Dispersion [mm]") +
  #ggtitle("Noiseless case, SNR = Inf dB") +
  #labs(subtitle = ("SNR = 30 [dB]" ) ) +
  #labs(subtitle = expression( paste("SNR =", infinity, "[dB]") ) ) +
  grids() +
  theme_bw()
ggsave( paste0("plot_",tagname, "SpaDis", ".pdf"),
        width = 4.5, height = 3, units = "in")

table %>%
  #drop_na() %>%
  filter( SNR=="30" ) %>%
  filter( Solver != "Proposed" ) %>%
  filter( AUROC_loc>0 ) %>%
  ggboxplot(x = "Solver",
            y = "AUROC_loc") +
  ylab("AUROC (local)") +
  #ggtitle("Noiseless case, SNR = Inf dB") +
  #labs(subtitle = ("SNR = 30 [dB]" ) ) +
  #labs(subtitle = expression( paste("SNR =", infinity, "[dB]") ) ) +
  grids() +
  theme_bw()
ggsave( paste0("plot_",tagname, "AUROC", ".pdf"),
        width = 4.5, height = 3, units = "in")

table %>%
  #drop_na() %>%
  filter( SNR=="30" ) %>%
  filter( Solver != "Proposed" ) %>%
  filter( AP_loc>0 ) %>%
  ggboxplot(x = "Solver",
            y = "AP_loc") +
  ylab("Average Precision (local)") +
  #ggtitle("Noiseless case, SNR = Inf dB") +
  #labs(subtitle = ("SNR = 30 [dB]" ) ) +
  #labs(subtitle = expression( paste("SNR =", infinity, "[dB]") ) ) +
  grids() +
  theme_bw()
ggsave( paste0("plot_",tagname, "AvgPre", ".pdf"),
        width = 4.5, height = 3, units = "in")

################################################################################



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

p <- ggbarplot(table[ table$SNR %in% c("Inf","30","20","10"), ], 
              x = "SNR",
              y = "DLE1",
              fill = "Solver") +
  ylab("Distance Localization Eror [mm]") +
  ggtitle(title_text) +
  grids()

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

