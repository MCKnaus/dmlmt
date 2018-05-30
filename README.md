# dmlmt
Double Machine Learning for Multiple Treatments

This code implements the Double Machine Learning approach (Chernozhukov et al., 2018) 
for multiple treatments following Farrell (2015). 
With modifications for sensitivity analysis and balancing checks as described and applied in Knaus (2018). 
The cross-validated Post-Lasso is based on the [glmnet](https://github.com/cran/glmnet) package.

## Example

The following example shows how the analysis works for a binary treatment variable. 
The data are taken from the [hdm](https://github.com/cran/hdm) package that is described in Chernozhukov, Hansen, & Spindler (2016).

```R
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
```

### Basic analysis
#### Binary treatment
The following code shows how to estimate the basic average potential outcomes and treatment effect.
Following the analysis in Knaus (2018) all nuisance paramteres are estimated using cross-validated Post-Lasso.

```R
stand_pl_bin <- dmlmt(X,D,Y)
```

One alternative is to use normal Lasso instead of Post-Lasso.
This is usually much faster because it does not require to solve a full logistic and OLS model at each considered Lambda.

```R
stand_l_bin <- dmlmt(X,D,Y,pl=FALSE)
```

#### Multiple treament
For expository purposes create a third treatment by randomly splitting the control group.

```R
D_mult <- D
D_mult[runif(length(D))*(1-D)>0.5] <- 2
table(D_mult)
```

Run the analysis with the multiple treatment.

```R
stand_pl_mult <- dmlmt(X,D_mult,Y)
stand_l_mult <- dmlmt(X,D_mult,Y,pl=FALSE)
```

### Extended analysis as proposed in Knaus (2018)
Run the analysis also for 1SE, 0.5SE, 0.5SE+ and 1SE+ rules to select the Lambda and calculate the weights for balancing checks.

```R
se_rules <- c(-1,-.5,.5,1)

# Binary
ext_pl_bin <- dmlmt(X,D,Y,se_rule=se_rules,w=TRUE)

# Example how to plot the results
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
```


## References

Chernozhukov, V., Chetverikov, D., Demirer, M., Duflo, E., Hansen, C., Newey, W., & Robins, J. (2018). Double/debiased machine learning for treatment and structural parameters. *The Econometrics Journal*, 21(1), C1-C68.

Chernozhukov, V., Hansen, C., & Spindler, M. (2016). *High-Dimensional Metrics in R*. [arXiv:1603.01700](https://arxiv.org/abs/1603.01700)

Farrell, M. H. (2015). Robust inference on average treatment effects with possibly more covariates than observations. *Journal of Econometrics*, 189(1), 1-23.

Knaus, M. C. (2018). A Double Machine Learning Approach to Estimate the Effects of Musical Practice on
Student's Skills, [arXiv](https://arxiv.org/abs/1805.10300)
