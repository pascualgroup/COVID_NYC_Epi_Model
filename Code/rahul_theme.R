##Default plot theme
rahul_theme = theme(
  axis.title.x = element_text(size = 16,
                              face = "bold",
                              color = "black"),
  axis.text.x = element_text(size = 14,
                             face = "bold",
                             color = "black"),
  axis.title.y = element_text(size = 16,
                              face = "bold",
                              color = "black"),
  legend.title = element_text(size = 16,
                              face = "bold",
                              color = "black"),
  legend.text = element_text(size = 14,
                             face = "bold",
                             color = "black"),
  axis.text.y = element_text(size = 14,
                             face = "bold",
                             color = "black"),
  strip.text = element_text(size = 16,
                            face = "bold",
                            color = "black")
)

rahul_man_figure_theme = theme(
  axis.title.x = element_text(size = 18,
                              face = "bold",
                              color = "black"),
  axis.text.x = element_text(size = 16,
                             face = "bold",
                             color = "black"),
  axis.title.y = element_text(size = 18,
                              face = "bold",
                              color = "black"),
  legend.title = element_text(size = 16,
                              face = "bold",
                              color = "black"),
  legend.text = element_text(size = 16,
                             face = "bold",
                             color = "black"),
  axis.text.y = element_text(size = 16,
                             face = "bold",
                             color = "black"),
  strip.text = element_text(size = 16,
                            face = "bold",
                            color = "black")
)

theme_white_background = theme(axis.line = element_line(colour = "black"),
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               panel.border = element_blank(),
                               panel.background = element_blank()) 
mean_legend_lab = scale_color_manual(name = "Legend", values = c("red","blue"),
                                     labels =
                                       c("Simulation Mean \n (Shaded Region: \n Mean +/- 2 SD)",
                                                "Observed")) 
median_legend_lab = scale_color_manual(name = "Legend", values = c("red","blue"),
                                     labels =
                                       c("Simulation Median \n (Shaded Region: \n 95% Quantiles)",
                                         "Observed")) 


median_legend_lab_with_fit_data = scale_color_manual(name = "Legend", values = c("red","blue", "darkgreen"),
                                       labels =
                                         c("Simulation Median \n (Shaded Region: \n 95% Quantiles)",
                                           "Observed", "Data Used For Fitting")) 