library(OmixM6A)
glori_se <- readRDS("/Users/weizhen/Documents/GitHub/OmixM6A/inst/extdata/glori_se.rds")
glori_se_bu <- OmixM6A(se = glori_se[,c(1,2)]) #Recommended method, a core idea to present in paper
metadata(glori_se_bu)
glori_se_bn <- classify_m6A_states(glori_se, method = "bmix")
metadata(glori_se_bn)
set.seed(123)
glori_se_bb <- classify_m6A_states(glori_se[sample.int(nrow(glori_se), 30000),7], method = "bbmix") #Very slow method, recommend to fit only on subset of rows
glori_se_bb2 <- OmixM6A(se = glori_se[sample.int(nrow(glori_se), 500),7], method = "bbmix", subsample_bbmix = NULL)
#Run EM iteration 114 ...
#Current alpha fg:  0.4691825 ...
#Current beta fg:  1.258089 ...
#Current alpha bg:  3.685351 ...
#Current beta bg:  170.0848 ...
#Current fg proportion:  0.3332777 ...
#Current bg proportion:  0.6667223 ...
#Converged in 114 iterations.

#See fdr of fitted models, can be understood as BH adjusted p value.
assays(glori_se_bu)$prob_bg

# Visualize the mixture of beta binomial distribution functions
# Parameters
alpha_fg <- 0.4691825
beta_fg <- 1.258089
alpha_bg <- 3.685351
beta_bg <- 170.0848
prop_fg <- 0.3332777
prop_bg <- 1 - prop_fg

# Create a sequence of values
x_values <- seq(0, 1, length.out = 1000)

# Calculate the density of each beta distribution
density_fg <- dbeta(x_values, alpha_fg, beta_fg)
density_bg <- dbeta(x_values, alpha_bg, beta_bg)

# Plot
plot(x_values, density_fg, type = 'l', col = 'red', lwd = 2,
     main = 'Density of Beta Distributions',
     xlab = 'Value', ylab = 'Density', ylim = c(0, max(density_bg)))
lines(x_values, density_bg, type = 'l', col = 'blue', lwd = 2)

# Adding a legend
legend("topright", legend=c("Foreground", "Background"), col=c("red", "blue"), lty=1, cex=0.8)

plot_bmixFit(glori_se)
plot_bumixFit(glori_se)
plot_bbmixFit(glori_se)
