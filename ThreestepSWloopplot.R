######################################################
## R code to accompany "Efficient designs for three-
## sequence stepped wedge trials with continuous recruitment"
## by R Hooper, O Quintin, J Kasza
##
## This R file includes code to generate a plot of 
## the generated results. Note that results were generated
## by splitting the computation between two machines.
######################################################

library("readxl")
library("ggpubr")
#library("devtools")
#devtools::install_github("matherealize/looplot")
library("looplot")


##Read in results file:
results <- read_excel("Simulation_results.xlsx")
##Read in the theoretical power
theoreticalpower <- read.csv(file="scenarios_power_20230214.csv") 

##Changes to ensure the nested loop plot looks good
resultssim <- cbind(theoreticalpower, as.numeric(results$V5))
resultssim <- as.data.frame(resultssim)
names(resultssim)[17] <- "RejPercentage"
names(resultssim)[3] <- "Design"
resultssim$rho <- as.factor(resultssim$rho)
resultssim$tau <- as.factor(resultssim$tau)

##Monte Carlo standard errors on empirical power
resultssim$powerSE <- sqrt((resultssim$RejPercentage)*(1-resultssim$RejPercentage)/1000)
resultssim$powerminusSE <- resultssim$RejPercentage - 2*sqrt((resultssim$RejPercentage)*(1-resultssim$RejPercentage)/1000)
resultssim$powerplusSE <- resultssim$RejPercentage + 2*sqrt((resultssim$RejPercentage)*(1-resultssim$RejPercentage)/1000)
resultssim$powerp <- resultssim$RejPercentage


NLP_plot <- nested_loop_plot(resdf =  resultssim, 
                                 x = "rho", steps = c("theta", "tau"), 
                                 grid_rows = "Design",
                                 steps_y_base = 0.65, steps_y_height = 0.05, steps_y_shift = 0.05,
                                 methods = c("powerp", "powerminusSE", "powerplusSE", "power"),
                                 line_linetypes = c(1, 1, 1, 2),
                                 colors = c("black", "darkgrey", "darkgrey", "darkred"),
                                 point_shapes = c("powerp" =16, "powerminusSE" = 32, "powerplusSE" = 32, "power" =17),
                                 point_size = 2, line_size = 0.75,
                                 steps_values_annotate = TRUE, steps_annotation_size = 3,
                                 steps_annotation_nudge = 0.2,
                                 y_name = "Power",
                                 x_name = "ICC",
                                 y_expand_add = c(0.1, NULL), 
                                 legend_name = "Legend", steps_names=c("Difference", "tau"),
                             hline_intercept = 0.8,
                                 legend_labels=c("Simulated power", "Sim. power -2se", "Sim. power +2se", "Theoretical power"),
                             replace_labels = list(
                               Design = c("n" = "Non-standard", 
                                          "y" = "Standard")
                             ),
                             post_processing = list(
                                   add_custom_theme = list(
                                     legend.position="bottom",
                                     axis.text.x = element_text(angle = -90, 
                                                                vjust = 0.5, 
                                                                size = 8) 
                                   )))


ggsave(NLP_plot, file= "NLP_plot.pdf", device = "pdf", width = 20, height = 15, units="cm")
ggsave(NLP_plot, file= "NLP_plot.png", device = "png", width = 20, height = 15, units="cm")