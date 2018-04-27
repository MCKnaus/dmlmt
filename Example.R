# Download current version from Github
library(devtools)
install_github(repo="MCKnaus/dmlmt")

# Get data
library(hdm)
data(pension)
Y = pension$tw; D = pension$p401
# Only main effects
X = model.matrix(~ -1 + i2 + i3 + i4 + i5 + i6 + i7 + a2 + a3 + a4 + a5 +
                fsize + hs + smcol + col + marr + twoearn + db + pira + hown, data = pension)
## Consider also interactions
# X = model.matrix(~ -1 + (i2 + i3 + i4 + i5 + i6 + i7 + a2 + a3 + a4 + a5 +
#                    fsize + hs + smcol + col + marr + twoearn + db + pira + hown)^2, data = pension)

# Create an indicator matrix with each column indicating one treatment
D_mat <- cbind(1-D,D)


###########################
#### Standard analysis ####
###########################

#### Analysis like in the paper
# Post-Lasso Logit for p(x) (slow)
library(dmlmt)
select_d <- post_lasso_cv(X,D,family = "binomial",output=FALSE)$names_pl
l_nm_d <- list(select_d[-1],select_d[-1])
gps <- gps_prep(X,D_mat,l_nm_d,cs=TRUE)

# Post-Lasso OLS for mu(x)
select_y0 <- post_lasso_cv(X[D==0,],Y[D==0],family = "gaussian",output=FALSE)$names_pl
select_y1 <- post_lasso_cv(X[D==1,],Y[D==1],family = "gaussian",output=FALSE)$names_pl
l_nm_y <- list(select_y0[-1], select_y1[-1])
y_mat <- y_prep(X,D_mat,Y,l_nm_y)

# Potential outcomes
PO <- PO_dmlmt(D_mat,Y,y_mat,gps$p,cs_i=gps$cs)
# ATE
ATE <- TE_dmlmt(PO$mu,gps$cs)



#### Alternative: Lasso instead of Post-Lasso for p(x) (substantially faster)
library(glmnet)
cvfit = cv.glmnet(X,D, family = "binomial")
ps_mat <- predict(cvfit, X, s = "lambda.min", type = "response")
ps_mat <- cbind(1-ps_mat,ps_mat)
gps <- gps_cs(ps_mat,D_mat)

select_y0 <- post_lasso_cv(X[D==0,],Y[D==0],family = "gaussian",output=FALSE)$names_pl
select_y1 <- post_lasso_cv(X[D==1,],Y[D==1],family = "gaussian",output=FALSE)$names_pl
l_nm_y <- list(select_y0[-1], select_y1[-1])
y_mat <- y_prep(X,D_mat,Y,l_nm_y)

# Potential outcomes
PO <- PO_dmlmt(D_mat,Y,y_mat,gps$p,cs_i=gps$cs)
# ATE
ATE <- TE_dmlmt(PO$mu,gps$cs)



###########################
#### Extended analysis ####
###########################

#### Different penalty term selection rules to check sensitivity
se_rules <- c(-1,-.5,.5,1)

select_d_se <- post_lasso_cv(X,D,family = "binomial",se_rule = se_rules,output=FALSE)
l_nm_d_se <- list(select_d_se$names_pl[-1],select_d_se$names_pl[-1])
gps <- gps_prep(X,D_mat,l_nm_d_se)

select_y0_se <- post_lasso_cv(X[D==0,],Y[D==0],family = "gaussian",se_rule = se_rules,output=FALSE)
select_y1_se <- post_lasso_cv(X[D==1,],Y[D==1],family = "gaussian",se_rule = se_rules,output=FALSE)
l_nm_y_se <- list(select_y0_se$names_pl[-1], select_y1_se$names_pl[-1])
y_mat <- y_prep(X,D_mat,Y,l_nm_y_se)

# Potential outcomes
PO <- PO_dmlmt(D_mat,Y,y_mat,gps$p,cs_i=gps$cs)
# ATE
ATE <- TE_dmlmt(PO$mu,gps$cs)

for (s in 1:length(select_d_se$names_Xse_pl)) {
  cat("\n\nSE rule used:",names(select_d_se$names_Xse_pl)[s],"\n")
  gps_se <- gps_prep(X,D_mat,list(select_d_se$names_Xse_pl[[s]][-1],select_d_se$names_Xse_pl[[s]][-1]),print=F)
  y_mat_se <- y_prep(X,D_mat,Y,list(select_y0_se$names_Xse_pl[[s]][-1], select_y1_se$names_Xse_pl[[s]][-1]))
  ### ATE
  PO <- PO_dmlmt(D_mat,Y,y_mat_se,gps_se$p,cs_i=gps_se$cs)
  ATE <- TE_dmlmt(PO$mu,gps_se$cs)
}


##### Getting the weights
w_ate <- PO_dmlmt_w(X,D_mat,Y,gps$p,l_nm_y,gps$cs)

# Also possible to enter previously calculated Nt x N weight matrix (e.g. from Random forest)
# w_dml <- rep(0,length(Y))
# for (i in 1:ncol(D_mat)) {
#   OLSw <- OLS_pred_w(add_intercept(X[,l_nm_y[[i]]]),D_mat[,i],Y,gps$cs)
#   w_dml <- w_dml + PO_dmlmt_w_gen(OLSw,D_mat[,i],Y,gps$p[,i],gps$cs)
# }

# Check balancing with the package of your choice, e.g. cobalt
library(cobalt)

balance <- bal.tab(as.data.frame(X), treat = D,weights=as.data.frame(rowSums(w_ate$w_dml)),method = "weighting",
                s.d.denom = "pooled", disp.v.ratio = TRUE, disp.ks = TRUE, un = TRUE)

plot <- love.plot(balance,abs = TRUE, line=TRUE, var.order="unadjusted")
