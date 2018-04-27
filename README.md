# dmlmt
Double Machine Learning for Multiple Treatments

This code implements the Double Machine Learning approach (Chernozhukov et al., 2018) 
for multiple treatments following Farrell (2015). 
With modifications for sensitivity analysis and balancing checks as described and applied in Knaus (2018). 
The cross-validated Post-Lasso is based on the [glmnet](https://github.com/cran/glmnet) package.

## Example

The following example shows for simplicity how the analysis works for a binary treatment variable. 
The data are taken from the [hdm](https://github.com/cran/hdm) package that is described in Chernozhukov, Hansen, & Spindler (2016).

```R
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
```

### Basic analysis
The following code shows how to estimate the basic average potential outcomes and treatment effect.
Following the analysis in Knaus (2018) all nuisance paramteres are estimated using cross-validated Post-Lasso.

```R
# Post-Lasso Logit for p(x) (slow)
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
```

One alternative is to use normal Lasso instead of Post-Lasso Logit.
This is usually much faster because it does not require to solve a full logistic model at each considered Lambda.

```R
require(glmnet)
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
```

### Extended analysis as proposed in Knaus (2018)
Run the analysis also for 1SE, 0.5SE, 0.5SE+ and 1SE+ rules to select the Lambda.

```R
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
```

Calculate the implied weights

```R
w_ate <- PO_dmlmt_w(X,D_mat,Y,gps$p,l_nm_y,gps$cs)
```

and use the package of you choice to assess balancing, e.g.

```R
require(cobalt)

balance <- bal.tab(as.data.frame(X), treat = D,weights=as.data.frame(rowSums(w_ate$w_dml)),method = "weighting",
                s.d.denom = "pooled", disp.v.ratio = TRUE, disp.ks = TRUE, un = TRUE)

plot <- love.plot(balance,abs = TRUE, line=TRUE, var.order="unadjusted")
```


## References

Chernozhukov, V., Chetverikov, D., Demirer, M., Duflo, E., Hansen, C., Newey, W., & Robins, J. (2018). Double/debiased machine learning for treatment and structural parameters. *The Econometrics Journal*, 21(1), C1-C68.

Chernozhukov, V., Hansen, C., & Spindler, M. (2016). 3High-Dimensional Metrics in R3. [arXiv:1603.01700](https://arxiv.org/abs/1603.01700)

Farrell, M. H. (2015). Robust inference on average treatment effects with possibly more covariates than observations. *Journal of Econometrics*, 189(1), 1-23.

Knaus, M. C. (2018). A Double Machine Learning Approach to Estimate the Effects of Musical Practice on
Student's Skills, in preparation.
