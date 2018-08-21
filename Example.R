# Download current version from Github
library(devtools)
install_github(repo="MCKnaus/dmlmt")
library(dmlmt)

# Get data
library(hdm)
data(pension)
Y = pension$tw; D = pension$p401
# Only main effects (toy example)
X = model.matrix(~ -1 + i2 + i3 + i4 + i5 + i6 + i7 + a2 + a3 + a4 + a5 +
                fsize + hs + smcol + col + marr + twoearn + db + pira + hown, data = pension)
## Consider also interactions if you have some time
# X = model.matrix(~ -1 + (i2 + i3 + i4 + i5 + i6 + i7 + a2 + a3 + a4 + a5 +
#                    fsize + hs + smcol + col + marr + twoearn + db + pira + hown)^2, data = pension)


###########################
#### Standard analysis ####
###########################

# Post-Lasso for nuisance parameters (might be slow)
stand_pl_bin <- dmlmt(X,D,Y)

# Lasso for nuisance parameters (faster)
stand_l_bin <- dmlmt(X,D,Y,pl=FALSE)

# Create multiple treatment by splitting control group
D_mult <- D
D_mult[runif(length(D))*(1-D)>0.5] <- 2
table(D_mult)

# Run the analysis with multiple treatments
stand_pl_mult <- dmlmt(X,D_mult,Y)
stand_l_mult <- dmlmt(X,D_mult,Y,pl=FALSE)


#################################
###     Extended Analysis     ###
#################################

# Consider different rules to select penalty term and calculate implied weights
se_rules <- c(-1,-.5,.5,1)

# Binary
ext_pl_bin <- dmlmt(X,D,Y,se_rule=se_rules,w=TRUE)

# Example how to plot the results
library(ggplot2)
df <- data.frame(SE_rule = factor(colnames(ext_pl_bin$SE_rule[[1]])
                                  ,levels = colnames(ext_pl_bin$SE_rule[[1]]))
                 ,coef = ext_pl_bin$SE_rule[[1]][1,],se = ext_pl_bin$SE_rule[[2]][1,])
j <- ggplot(df, aes(SE_rule, coef, ymin = coef-se, ymax = coef+se)) +
             geom_errorbar() + geom_point()

# Example how to check balancing with the package of your choice, e.g. cobalt
library(cobalt)
balance <- bal.tab(as.data.frame(X), treat = D,weights=ext_pl_bin$weights,method = "weighting",
                    s.d.denom = "pooled", disp.v.ratio = TRUE, disp.ks = TRUE, un = TRUE)
plot <- love.plot(balance,abs = TRUE, line=TRUE, var.order="unadjusted")

# Multiple
ext_pl_mult <- dmlmt(X,D_mult,Y,se_rule=se_rules,w=TRUE)


#################################
###     Generic ML input      ###
#################################

# Example with Random Forest instead of Lasso
library(grf)

# Initialize nuisance matrices
values <- sort(unique(D_mult))
ps_mat <- t_mat <- y_mat <- matrix(NA,length(Y),length(values))

# Get nuisance parameter predictions
for (tr in 1:length(values)){
  t_mat[,tr] <- as.numeric(D_mult == values[tr])
  rf_p <- regression_forest(X,t_mat[,tr])
  ps_mat[,tr] <- predict(rf_p, X)$predictions
  rf_y <- regression_forest(X[t_mat[,tr] == 1,],Y[t_mat[,tr] == 1])
  y_mat[,tr] <- predict(rf_y, X)$predictions
}

# Calculate generalized p-score and enforce common support
rf_gps <- gps_cs(ps_mat,t_mat)

# Potential outcomes
rf_PO <- PO_dmlmt(t_mat,Y,y_mat,rf_gps$p,cs_i=rf_gps$cs)
# ATE
rf_ATE <- TE_dmlmt(rf_PO$mu,rf_gps$cs)
