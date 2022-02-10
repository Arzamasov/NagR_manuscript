# this script contains functions for quantitative analysis of EMSA gel imaging

# 4-parameter logistic (4PL) model function
# in the 4PL model, lower limit is fixed to 0 and upper limit is fixed to 1
# input: data frame with EMSA data
fourPL_model <- function(site_data){
  drm(percent_shift~conc, 
      data = site_data, 
      fct= LL.2(names = c("Slope", "EC50")))
  }

# plotting function
# plots: (1) EMSA quantification data, (2) fit from the 4PL model with 95% confidence intervals
# input:
# (1) drc model
# (2) tibble with EMSA data
# (3) max protein concentration used (nM)
# example: plot_Blon_0879 <- plot_fit(model_Blon_0879, data_Blon_0879, 1000)
plot_fit <- function(model, site_data, limit) {
  # new concentration levels as support for the line
  conf_site <- expand.grid(conc=exp(seq(log(0.5), log(limit), length=100)))
  # calculate predicted values and respective 95% confidence intervals from the model using the "predict" function
  pm_site <- predict(model, newdata=conf_site, interval="confidence")
  # new data with predictions/confidence intervals
  conf_site$p <- pm_site[,1]
  conf_site$pmin <- pm_site[,2]
  conf_site$pmax <- pm_site[,3]
  EC50 <- model$coefficients[2]
  # plot the graph
  graph <- ggplot(data = site_data, aes(x = conc, y = percent_shift)) +
    geom_line(data=conf_site, aes(x=conc, y=p)) +
    geom_ribbon(data=conf_site, aes(x=conc, y=p, ymin=pmin, ymax=pmax), alpha=0.2) +
    geom_point(size=4, pch=21, fill="orangered2") +
    geom_point(aes(x = EC50, y = 0.5), size = 7, pch=21, fill="lightskyblue") +
    annotate("text", x=200, y=0.12, 
             label= paste("EC50 =",signif(EC50, 4), "nM"),
             size = 4.5) +
    scale_x_continuous(trans=scales::pseudo_log_trans(base=10),
                       limits=c(0, limit),
                       breaks=c(0, 5, 10, 50, 100, 500, 1000)) +
    xlab("NagR (nM)") +
    ylab("Shift (%)") + 
    theme_minimal_grid(12)
  return(graph)
}