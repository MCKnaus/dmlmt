#' Double machine learning for binary and multiple treatments'
#'
#'
#' This function estimates treatment effects for binary and multiple treatments
#' using Double Machine Learning.
#'
#' @param x Matrix of covariates (N x p matrix)
#' @param t Vector of treament indicators. Will be ordered from 0 to T-1.
#' @param y Vector of outcomes
#' @param pl If TRUE Post-Lasso is used to estimate nuisance parameters, if FALSE Lasso
#' @param cs If TRUE, common support will be checked
#' @param q Quantile used for enforcing common support
#' @param cl Vector with cluster variables can be provided
#' @param print If TRUE, supporting information is printed
#' @param se_rule If not NULL, define, e.g., c(-1,1) to get 1SE and 1SE+ rule
#' @param w If TRUE, implied weights are calculated (only if pl=TRUE)
#' @param parallel If TRUE, cross-validation of \code{\link{post_lasso_cv}} parallelized
#' @param ... Pass \code{\link{glmnet}} and \code{\link{post_lasso_cv}} options
#'
#' @return \code{gps_prep} returns a list containing a N x # of treatments matrix
#'          with generalized propensity scores (\code{p}) and a boolean indicating common
#'          and a boolean indicating common support (\code{cs}).
#' @export

dmlmt <- function(x,t,y,pl=TRUE,cs=TRUE,q=1,cl=NULL,print=FALSE,se_rule = NULL,w=FALSE,parallel=FALSE,...) {

  # Checks
  if (!isTRUE(pl) & !is.null(se_rule)) stop("Different standard error rules only implemented for Post-Lasso")
  if (!is.null(se_rule) & !(is.numeric(se_rule))) stop("Please provide numeric SE rules (e.g. c(-1,1))")
  if (!isTRUE(pl) & isTRUE(parallel)) warning("Parallelization currently not implemented for Lasso"); parallel=FALSE

  # Important parameters
  n <- length(y)
  num_t <- length(table(t))

  # Check whether treatment is binary or multiple
  binary <- all(t %in% 0:1)

  SE_rule <- NULL
  if (!is.null(se_rule)) l_se_d <- vector("list",num_t)
  if (!is.null(se_rule)) l_se_y <- vector("list",num_t)

  if (binary) {
    cat("\n Binary treatment\n\n")
    t_mat <- cbind(1-t,t)
    if (isTRUE(pl)) {
      sel_d <- post_lasso_cv(x,t,family = "binomial",output=print,se_rule=se_rule,parallel=parallel,...)
      l_nm_d <- list(sel_d$names_pl[-1],sel_d$names_pl[-1])
      if (!is.null(se_rule)) l_se_d <- list(sel_d$names_Xse_pl,sel_d$names_Xse_pl)
      gps <- gps_prep(x,t_mat,l_nm_d,cs=cs,q=q,print=print)

      # Post-Lasso OLS for mu(x)
      sel_y0 <- post_lasso_cv(x[t==0,],y[t==0],output=print,se_rule=se_rule,parallel=parallel,...)
      sel_y1 <- post_lasso_cv(x[t==1,],y[t==1],output=print,se_rule=se_rule,parallel=parallel,...)
      l_nm_y <- list(sel_y0$names_pl[-1], sel_y1$names_pl[-1])
      if (!is.null(se_rule)) l_se_y <- list(sel_y0$names_Xse_pl,sel_y1$names_Xse_pl)
      y_mat <- y_prep(x,t_mat,y,l_nm_y)

    } else {
      cvfit_p <- cv.glmnet(x,t, family = "binomial",parallel=parallel,...)
      ps_mat <- predict(cvfit_p, x, s = "lambda.min", type = "response")
      ps_mat <- cbind(1-ps_mat,ps_mat)
      gps <- gps_cs(ps_mat,t_mat,cs=cs,q=q,print=print)

      y_mat <- matrix(NA,n,num_t)
      for (tr in 1:num_t) {
        cvfit_y <- cv.glmnet(x[t_mat[,tr]==1,],y[t_mat[,tr]==1],parallel=parallel,...)
        y_mat[,tr] <- predict(cvfit_y, x, s = "lambda.min")
      }
    }
  }  else { # end binary
    cat("\n Multiple treatment\n\n")
    t_mat <- matrix(NA,n,num_t)
    values <- as.numeric(names(table(t)))
    for (tr in 1:num_t){
      t_mat[,tr] <- as.numeric(t == values[tr])
    }

    if (isTRUE(pl)) {
      list_t <- vector("list",num_t)
      l_nm_y <- vector("list",num_t)

      for (tr in 1:num_t){
        sel_d <- post_lasso_cv(x,t_mat[,tr],family = "binomial",output=print,se_rule=se_rule,parallel=parallel,...)
        list_t[[tr]] <- sel_d$names_pl[-1]
        if (!is.null(se_rule)) l_se_d[[tr]] <- sel_d$names_Xse_pl
        sel_y <- post_lasso_cv(x[t_mat[,tr]==1,],y[t_mat[,tr]==1],output=print,se_rule=se_rule,parallel=parallel,...)
        l_nm_y[[tr]] <- sel_y$names_pl[-1]
        if (!is.null(se_rule)) l_se_y[[tr]] <- sel_y$names_Xse_pl
      }
      gps <- gps_prep(x,t_mat,list_t,cs=cs,q=q,print=print)
      y_mat <- y_prep(x,t_mat,y,l_nm_y)
    } else {
      ps_mat <- matrix(NA,n,num_t)
      y_mat <- matrix(NA,n,num_t)
      for (tr in 1:num_t){
        cvfit_p <- cv.glmnet(x,t_mat[,tr], family = "binomial",...)
        ps_mat[,tr] <- predict(cvfit_p, x, s = "lambda.min", type = "response")
        cvfit_y <- cv.glmnet(x[t_mat[,tr]==1,],y[t_mat[,tr]==1],parallel=parallel,...)
        y_mat[,tr] <- predict(cvfit_y, x, s = "lambda.min")
      }
      gps <- gps_cs(ps_mat,t_mat,cs=cs,q=q,print=print)
    }
  }

  # Potential outcomes
  PO <- PO_dmlmt(t_mat,y,y_mat,gps$p,cs_i=gps$cs,cl=cl)
  # ATE
  ATE <- TE_dmlmt(PO$mu,gps$cs,cl=cl)

  # Calculate the weights
  if (isTRUE(w)) weights <- rowSums(PO_dmlmt_w(x,t_mat,y,gps$p,l_nm_y,gps$cs)$w_dml)


  ### Different penalty term rules
  if (!is.null(se_rule)) {
    TE <- matrix(NA,nrow(ATE),length(se_rule))
    SE <- matrix(NA,nrow(ATE),length(se_rule))
    rownames(TE) <- rownames(ATE)
    rownames(SE) <- rownames(ATE)
    colnames(TE) <- names(sel_d$names_Xse_pl)
    colnames(SE) <- names(sel_d$names_Xse_pl)

    for (s in 1:length(se_rule)) {
      cat("\n\nSE rule used:",names(sel_d$names_Xse_pl)[s],"\n")
      list_se_d <- vector("list",num_t)
      list_se_y <- vector("list",num_t)
      for (tr in 1:num_t) {
        list_se_d[[tr]] <- l_se_d[[tr]][[s]][-1]
        list_se_y[[tr]] <- l_se_y[[tr]][[s]][-1]
      }
      se_gps <- gps_prep(x,t_mat,list_se_d,print=F)
      se_y_mat <- y_prep(x,t_mat,y,list_se_y)

      se_PO <- PO_dmlmt(t_mat,y,se_y_mat,se_gps$p,cs_i=se_gps$cs,cl=cl)
      se_ATE <- TE_dmlmt(se_PO$mu,se_gps$cs,cl=cl)
      TE[,s] <- se_ATE[,1,drop=F]
      SE[,s] <- se_ATE[,2,drop=F]
    }
    # Save treatment effects and standard errors
    TE <- cbind(TE[,se_rule<0,drop=F],ATE[,1,drop=F],TE[,se_rule>0,drop=F])
    colnames(TE)[sum(se_rule<0)+1] <- "Min"
    SE <- cbind(SE[,se_rule<0,drop=F],ATE[,2,drop=F],SE[,se_rule>0,drop=F])
    colnames(SE)[sum(se_rule<0)+1] <- "Min"
    SE_rule <- list(TE,SE)
  }

  ## Return results
  list("ATE" = ATE,"PO" = PO$results,"SE_rule"=SE_rule,"weights"=weights)
}



#' This function prepares the (generalized) propensity scores
#' and checks the common support after variable selection.
#'
#' @param x Matrix of covariates (N x p matrix)
#' @param t Matrix of binary treament indicators (N x # of treatments matrix)
#'          each column contains a binary indicator for each tratment
#' @param nm_list List with length = # of treatments containing list with the
#'          selected variable names for each treatment level
#' @param cs If TRUE, common support will be checked
#' @param q Quantile used for enforcing common support
#' @param print If TRUE, descriptives for p-scores and common support shown
#' @import psych
#'
#' @return \code{gps_prep} returns a list containing a N x # of treatments matrix
#'          with generalized propensity scores (\code{p}) and a boolean indicating common
#'          and a boolean indicating common support (\code{cs}).
#' @export

gps_prep <- function(x,t,nm_list,cs=TRUE,q=1,print=TRUE) {

  # Retrieve important info
  n <- nrow(x)
  num_t <- ncol(t)

  # Intialize matrices
  p_mat <- matrix(NA,n,num_t)
  colnames(p_mat) <- sprintf("Treatment %d",0:(num_t-1))
  minmax <- matrix(NA,2,num_t)
  cs_mat <- matrix(NA,n,num_t)

  for (i in 1:num_t) {
    xx <- add_intercept(x[,nm_list[[i]],drop=F])
    fit_val <- stats::glm.fit(xx,t[,i],family=binomial(link="logit"))
    if (!fit_val$converged) warning("Logit did not converge: GPS might be implausible, check!", call. = FALSE)
    p_mat[,i] <- fit_val$fitted.values

    for (j in 1:num_t) {
      minmax[1,j] <- stats::quantile(p_mat[t[,j]==1,i],1-q)
      minmax[2,j] <- stats::quantile(p_mat[t[,j]==1,i],q)
    }
    cs_mat[,i] <- (p_mat[,i] < max(minmax[1,]) | p_mat[,i] > min(minmax[2,]))
  }

  cs_ind <- rep(TRUE,n)

  if (isTRUE(cs)) cs_ind <- !apply(cs_mat,1,any)

  if (isTRUE(print)) {
    cat("\n\nPscores\n")
    print(psych::describe(p_mat))

    if (isTRUE(cs)) {
      cat("\nOff support\n", toString(sum(!cs_ind)))
      cat("\n\nPscores on support\n")
      print(psych::describe(p_mat[cs_ind,]))
    }
  }

  list("p"=p_mat,"cs"=cs_ind)
}



#' This function checks common support for a user provided matrix of
#' (generalized) propensity scores
#'
#' @param p_mat N x # of treatments matrix with (generalized) propensity scores
#' @param t Matrix of binary treament indicators (N x # of treatments matrix)
#'          each column contains a binary indicator for each tratment
#' @param q Quantile used for enforcing common support
#' @param print If TRUE, descriptives for p-scores and common support shown
#' @param cs If TRUE common support enforced, if FALSE only boolean vector added
#' @import psych
#' @return \code{gps_cs} returns a list containing a N x # of treatments matrix
#'          with generalized propensity scores (\code{p}) and a boolean indicating common
#'          and a boolean indicating common support (\code{cs}).
#' @export

gps_cs <- function(p_mat,t,q=1,print=TRUE,cs=TRUE) {
  # Retrieve important info
  n <- nrow(p_mat)
  num_t <- ncol(t)

  # Intialize matrices
  minmax <- matrix(NA,2,num_t)
  cs_mat <- matrix(NA,n,num_t)

  for (i in 1:num_t) {
    for (j in 1:num_t) {
      minmax[1,j] <- stats::quantile(p_mat[t[,j]==1,i],1-q)
      minmax[2,j] <- stats::quantile(p_mat[t[,j]==1,i],q)
    }
    cs_mat[,i] <- (p_mat[,i] < max(minmax[1,]) | p_mat[,i] > min(minmax[2,]))
  }

  cs_ind <- rep(TRUE,n)
  if (isTRUE(cs)) cs_ind <- !apply(cs_mat,1,any)

  if (isTRUE(print)) {
    cat("\n\nPscores\n")
    print(psych::describe(p_mat))

    if (isTRUE(cs)) {
      cat("\nOff support\n", toString(sum(!cs_ind)))
      cat("\n\nPscores on support\n")
      print(psych::describe(p_mat[cs_ind,]))
    }
  }

  list("p"=p_mat,"cs"=cs_ind)
}


#' This function prepares the fitted values of the outcome
#' after variable selection.
#'
#' @param x Matrix of covariates (N x p matrix)
#' @param t Matrix of binary treament indicators (N x # of treatments matrix)
#'          each column contains a binary indicator for each tratment
#' @param y Vector of outcomes
#' @param nm_list List with length = # of treatments containing list with the
#'          selected variable names for each treatment level
#'
#' @return \code{y_prep} returns a N x # of treatments matrix with the fitted values
#'          for each treatment level.
#' @export


y_prep <- function(x,t,y,nm_list) {

  # Retrieve important info
  n <- nrow(t)
  num_t <- ncol(t)
  y <- as.matrix(y,ncol=1)

  # Intialize matrix
  y_mat <- matrix(NA,n,num_t)

  for (i in 1:num_t) {
    xx <- add_intercept(x[,nm_list[[i]],drop=FALSE])
    coef <- glm.fit(xx[t[,i]==1,],y[t[,i]==1,])$coefficients
    # Check for NA coefficients
    if (any(is.na(coef)) == TRUE) {
      xx <- xx[,!is.na(coef)]
      coef <- coef[!is.na(coef)]
    }
    y_mat[,i] <- xx%*%coef

  }

  return(y_mat)
}


#' This function calculates the potential outcomes for all treatment levels.
#'
#' @param t Matrix of binary treament indicators (N x # of treatments matrix)
#'          each column contains a binary indicator for each tratment
#' @param y Vector of outcomes
#' @param y_mat N x # of treatments matrix with fitted outcome values
#' @param p_mat N x # of treatments matrix with (generalized) propensity scores
#' @param cs_i If not NULL, boolean vector to indicate that observation is on support
#' @param cl If not NULL, vector with cluster variables
#' @param print If TRUE, print results
#'
#' @return \code{PO_dmlmt} returns a list containing the \code{results} and a matrix
#'          with the individual values of the effcicient score (\code{mu}).
#' @export

PO_dmlmt <- function(t,y,y_mat,p_mat,cs_i=NULL,cl=NULL,print=TRUE) {

  # Retrieve important info
  n <- nrow(t)
  num_t <- ncol(t)
  if (is.null(cs_i)) rep(TRUE,n)

  # Initialize matrices
  w_ipw <- matrix(0,n,num_t)
  mu_mat <- matrix(NA,n,num_t)
  res <- matrix(NA,num_t,2)
  rownames(res) <- sprintf("Treatment %d",0:(num_t-1))
  colnames(res) <- c("PO","SE")

  for (i in 1:num_t) {
    w_ipw[cs_i,i] <- as.matrix(t[cs_i,i] / p_mat[cs_i,i],ncol=1)
    w_ipw[cs_i,i] <- norm_w_to_n(w_ipw[cs_i,i,drop=F])
  }

  for (i in 1:num_t) {
    # Potential outcome ES for individual
    mu_mat[cs_i,i] <- w_ipw[cs_i,i] * (y[cs_i] - y_mat[cs_i,i])  + y_mat[cs_i,i]
  }

  # Calculate Mean PO
  res[,1] <- colMeans(mu_mat,na.rm=TRUE)

  # Calculate SE for PO
  if (is.null(cl)) {
    for (i in 1:num_t) {
      res[i,2] <- sqrt(mean((mu_mat[cs_i,i] - mean(mu_mat[cs_i,i]))^2) / n)
    }
  } else {
    for (i in 1:num_t) {
      res[i,2] <- sqrt(sum(tapply(mu_mat[cs_i,i] - mean(mu_mat[cs_i,i]), cl[cs_i], sum)^2) / n^2)
    }
  }

  if (isTRUE(print)) {
    cat("\n\n Potential outcomes:\n")
    stats::printCoefmat(res)
  }

  list("results"=res,"mu"=mu_mat)
}


#' This function calculates the weights impled by Double Machine learning.
#'
#' @param x Matrix of covariates (N x p matrix)
#' @param t Matrix of binary treament indicators (N x # of treatments matrix)
#'          each column contains a binary indicator for each tratment
#' @param y Vector of outcomes
#' @param p_mat N x # of treatments matrix with (generalized) propensity scores
#' @param nm_list List with length = # of treatments containing list with the
#'          selected variable names in the outcome regression for each treatment level
#' @param cs_i If not NULL, boolean vector to indicate that observation is on support
#' @param cl If not NULL, vector with cluster variables
#' @import RandomFieldsUtils
#'
#' @return \code{PO_dmlmt_w} returns a list containing the Double Machine Learning Weights,
#'         the IPW weights, the OLS weights, and the adjustment weights.
#' @export

PO_dmlmt_w <- function(x,t,y,p_mat,nm_list,cs_i=NULL,cl=NULL) {

  # Retrieve important info
  n <- nrow(t)
  num_t <- ncol(t)
  y <- as.matrix(y,ncol=1)
  if (is.null(cs_i)) rep(TRUE,n)

  # Predict y's
  w_ipw <- matrix(0,n,num_t)
  w_ols <- matrix(0,n,num_t)
  w_adj <- matrix(0,n,num_t)

  for (i in 1:num_t) {
    xx <- add_intercept(x[,nm_list[[i]],drop=FALSE])
    # Remove multicollinear
    coef <- stats::glm.fit(xx[t[,i]==1,],y[t[,i]==1,])$coefficients
    # Check for NA coefficients
    if (any(is.na(coef)) == TRUE) {
      xx <- xx[,!is.na(coef),drop=F]
    }

    # IPW weights
    w_ipw[cs_i,i] <- as.matrix(t[cs_i,i] / p_mat[cs_i,i],ncol=1)
    w_ipw[cs_i,i] <- norm_w_to_n(w_ipw[cs_i,i,drop=F])

    # Get X'X for treated
    XtX <- crossprod(xx[t[,i]==1,])
    # Get X(X'X)-1 for treated
    XXtX <- xx[t[,i]==1,] %*% solvex(XtX)

    for (r in 1:n) {
      if (cs_i[r]==TRUE) {
        w_ol <- matrix(0,n,1)
        w_ol[t[,i]==1,] <- XXtX %*% xx[r,]
        w_ols[,i] <- w_ols[,i] + w_ol
        w_adj[,i] <- w_adj[,i] + w_ol * w_ipw[r,i]
      }
    }
  } # End t

  # Calculate weight matrix
  w_mat <- w_ipw + w_ols - w_adj

  list("w_dml"=w_mat,"w_ols"=w_ols,"w_ipw"=w_ipw,"w_adj"=w_adj)
}



#' This function calculates the OLS weights for the fitted values of
#' each observation
#'
#' @param x Matrix of covariates (N x p matrix)
#' @param t Vector of binary treament indicator
#' @param y Vector of outcomes
#' @param cs_i If not NULL, boolean vector to indicate that observation is on support
#'
#' @return Vector of weights
#' @example
#' # Also possible to enter previously calculated Nt x N weight matrix (e.g. from Random forest)
#' \dontrun{
#' w_dml <- rep(0,length(Y))
#' for (i in 1:ncol(D_mat)) {
#'   OLSw <- OLS_pred_w(add_intercept(X[,l_nm_y[[i]]]),D_mat[,i],Y,gps$cs)
#'   w_dml <- w_dml + PO_dmlmt_w_gen(OLSw,D_mat[,i],Y,gps$p[,i],gps$cs)
#' }
#' }
#'
#' @export


OLS_pred_w <- function(x,t,y,cs_i=NULL) {
  # Retrieve important info
  n <- length(t)
  n_t <- sum(t)
  y <- as.matrix(y,ncol=1)
  if (is.null(cs_i)) rep(TRUE,n)

  # Matrix to be filled
  w_mat <- matrix(0,n_t,n)

  # Remove multicollinear
  coef <- stats::glm.fit(x[t==1,],y[t==1,])$coefficients
  if (any(is.na(coef)) == TRUE) {
    x <- x[,!is.na(coef),drop=F]
  }

  # Get X'X for treated
  XtX <- crossprod(x[t==1,])
  # Get X(X'X)-1 for treated
  XXtX <- x[t==1,] %*% solvex(XtX)

  for (r in 1:n) {
    if (cs_i[r]==TRUE) {
      w_mat[,r] <- XXtX %*% x[r,]
    }
  }
  return(w_mat)
}


#' This function calculates the weights impled by Double Machine learning
#' if outcome weights are already calculated.
#'
#' @param yw_mat Matrix containing weights for each observation (Nt x N matrix)
#' @param t Vector of binary treament indicator
#' @param y Vector of outcomes
#' @param p Vector of propensity scores
#' @param cs_i If not NULL, boolean vector to indicate that observation is on support
#'
#' @return Vector of weights
#' @export

PO_dmlmt_w_gen <- function(yw_mat,t,y,p,cs_i=NULL) {
  # Retrieve important info
  n <- length(t)
  n_t <- sum(t)
  t <- as.matrix(t,ncol=1)
  if (is.null(cs_i)) cs_i <- rep(TRUE,n)

  w_ipw <- matrix(0,n,1)
  w_adj <- matrix(0,n,1)
  # IPW weights
  w_ipw[cs_i,1] <- as.matrix(t[cs_i,1] / p[cs_i],ncol=1)
  w_ipw[cs_i,1] <- norm_w_to_n(w_ipw[cs_i,1,drop=F])

  w_ols <- as.matrix(rep(0,n))
  w_ols[t==1,] <- rowSums(yw_mat)

  for (r in 1:n) {
    if (cs_i[r]==TRUE) {
      w_adj[t==1,1] <- w_adj[t==1,1] + yw_mat[,r] * w_ipw[r,1]
    }
  }

  # Calculate weight matrix
  w_mat <- w_ipw + rowSums(w_ols) - w_adj

  return(w_mat)
}



#' This function calculates the average treatment effects for all combinations.
#' Requires input from \code{PO_dmlmt}.
#'
#' @param mu Matrix mu with the individual values of the effcicient score from \code{PO_dmlmt}
#' @param cs_i If not NULL, boolean vector to indicate that observation is on support
#' @param cl If not NULL, vector with cluster variables
#' @param print If TRUE, print results
#'
#' @return \code{TE_dmlmt} returns the results average reatment effects.
#' @export

TE_dmlmt <- function(mu,cs_i=NULL,print=TRUE,cl=NULL) {

  # Retrieve important info
  n <- nrow(mu)
  num_t <- ncol(mu)
  nm_t <- sprintf("T%d",0:(num_t-1))
  if (is.null(cs_i)) cs_i <- rep(TRUE,n)

  # Initialize results
  res <- matrix(NA,sum(1:(num_t-1)),4)
  colnames(res) <- c("TE","SE","t","p")
  rownames(res) <- rep("Platzhalter",nrow(res))

  pos <- 1
  for (i in 1:(num_t-1)) {
    loc <- i+1
    for (j in loc:(num_t)) {
      eif_i <- mu[cs_i,j] - mu[cs_i,i]
      res[pos,1] <- mean(eif_i)
      if (is.null(cl)) {
        res[pos,2] <- sqrt(mean((eif_i-mean(eif_i))^2) / n)
      } else {
        res[pos,2] <-  sqrt(sum(tapply(eif_i - mean(eif_i), cl[cs_i], sum)^2) / n^2)
      }

      rownames(res)[pos] <- paste0(nm_t[j]," - ",nm_t[i])
      pos <- pos + 1
    }
  }

  # t-stat
  res[,3] <- res[,1] / res[,2]
  # p-value
  res[,4] <- 2 * stats::pt(abs(res[,3]),n,lower = FALSE )

  if (isTRUE(print)) {
    cat("\nAverage effects\n")
    stats::printCoefmat(res,has.Pvalue = TRUE)
    cat("\n# of obs on / off support:",toString(sum(cs_i))," / ",toString(sum(!cs_i)) ,"\n")
  }
  return(res)
}




#################################################################################
#################################################################################
################### Functions for ATET (not used and tested) ####################
#################################################################################
#################################################################################
#
# gps_prep_atet <- function(x,t,nm_list,pos_t=2,cs=TRUE,q=1) {
#
#   # Retrieve important info
#   n <- nrow(x)
#   num_t <- ncol(t)
#
#   # Intialize matrices
#   p_mat <- matrix(NA,n,num_t)
#   minmax <- matrix(NA,2,1)
#   cs_mat <- matrix(NA,n,num_t)
#
#   for (i in 1:num_t) {
#     xx <- add_intercept(x[,nm_list[[i]],drop=F])
#     fit_val <- glm.fit(xx,t[,i],family=binomial(link="logit"))
#     if (!fit_val$converged) warning("Logit did not converge: GPS might be implausible, check!", call. = FALSE)
#     p_mat[,i] <- fit_val$fitted.values
#   }
#
#   minmax[1,1] <- quantile(p_mat[t[,pos_t]==1,pos_t],1-q)
#   minmax[2,1] <- quantile(p_mat[t[,pos_t]==1,pos_t],q)
#
#   if (isTRUE(cs)) cs_ind <- !(p_mat[,pos_t] <= minmax[1,] | p_mat[,pos_t] >= minmax[2,])
#
#   cat("\n\nPscores\n")
#   print(psych::describe(p_mat))
#
#   if (isTRUE(cs)) {
#     cat("\nOff support\n", toString(sum(!cs_ind)))
#     cat("\n\nPscores on support\n")
#     print(psych::describe(p_mat[cs_ind,]))
#   }
#
#
#   list("p"=p_mat,"cs"=cs_ind)
# }
#
#
#
# #################################################################################
# #################################################################################
#
# PO_dmlmt_atet <- function(t,y,y_mat,p_mat,cs_i=NULL,pos_t=2,cl=NULL,print=TRUE) {
#
#   # Retrieve important info
#   n <- nrow(t)
#   num_t <- ncol(t)
#   if (is.null(cs_i)) rep(TRUE,n)
#
#   # Initialize matrices
#   w_ipw <- matrix(0,n,num_t)
#   mu_mat <- matrix(NA,n,num_t)
#   res <- matrix(NA,num_t,2)
#   rownames(res) <- sprintf("Treatment %d",0:(num_t-1))
#   rownames(res)[pos_t] <- paste0(rownames(res)[pos_t],"*")
#   colnames(res) <- c("PO","SE")
#
#   # Get pscore of treated
#   ps_t <- p_mat[cs_i,pos_t,drop=F]
#   p_t <- mean(t[cs_i,2])
#
#   for (i in 1:num_t) {
#     # Potential outcome ES for individual
#     mu_mat[cs_i,i] <- ((t[cs_i,pos_t]==1) * y_mat[cs_i,i]) / p_t +
#                       ((t[cs_i,i]==1) * ps_t * (y[cs_i] - y_mat[cs_i,i])) / (p_t * p_mat[cs_i,i])
#   }
#
#   # Calculate Mean PO
#   res[,1] <- colMeans(mu_mat,na.rm=TRUE)
#
#   # Calculate SE for PO
#   if (is.null(cl)) {
#     for (i in 1:num_t) {
#       res[i,2] <- sqrt(mean((mu_mat[cs_i,i] - mean(mu_mat[cs_i,i]))^2) / n)
#     }
#   } else {
#     for (i in 1:num_t) {
#       res[i,2] <- sqrt(sum(tapply(mu_mat[cs_i,i] - mean(mu_mat[cs_i,i]), cl[cs_i], sum)^2) / n^2)
#     }
#   }
#
#   if (isTRUE(print)) printCoefmat(res)
#
#   list("results"=res,"mu"=mu_mat)
# }
#
#
# #################################################################################
# #################################################################################
#
# ATET_dmlmt <- function(mu,cs_i,pos_t=2,print=TRUE,cl=NULL) {
#
#   num_t <- ncol(mu)
#   nm_t <- sprintf("T%d",0:(num_t-1))
#
#   # Initialize results
#   res <- matrix(NA,num_t-1,4)
#   colnames(res) <- c("TE","SE","t","p")
#   rownames(res) <- rep("Platzhalter",nrow(res))
#
#   pos <- 1
#   for (i in 1:(num_t)) {
#     if (i != pos_t) {
#       print("1")
#       eif_i <- mu[cs_i,pos_t] - mu[cs_i,i]
#       print("2")
#       res[pos,1] <- mean(eif_i)
#       print("3")
#       if (is.null(cl)) {
#         res[pos,2] <- sqrt(mean((eif_i-mean(eif_i))^2) / n)
#       } else {
#         res[pos,2] <-  sqrt(sum(tapply(eif_i - mean(eif_i), cl[cs_i], sum)^2) / n^2)
#       }
#       print("4")
#       rownames(res)[pos] <- paste0(nm_t[pos_t],"* - ",nm_t[i])
#       pos <- pos + 1
#     }
#   }
#
#   # t-stat
#   res[,3] <- res[,1] / res[,2]
#   # p-value
#   res[,4] <- 2 * pt(abs(res[,3]),n,lower = FALSE )
#
#   if (isTRUE(print)) {
#     cat("\nAverage effects\n")
#     printCoefmat(res,has.Pvalue = TRUE)
#     cat("\n# of obs on / off support",toString(sum(cs_i))," / ",toString(sum(!cs_i)) ,"\n")
#   }
#   return(res)
# }
