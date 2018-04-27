#' This function uses the \code{\link{glmnet}} package to estimate the coefficient paths
#' and cross-validates Lasso AND Post-Lasso
#'
#' @param x Matrix of covariates (number of observations times number of covariates matrix)
#' @param y vector of outcomes
#' @param w vector of weights
#' @param kf number of folds in k-fold CV
#' @param family Outcome family, currenlty "gaussian" and "binomial" supported
#' @param seed CV samples are drawn starting with a specified seed
#' @param output If TRUE, output and graphs are printed
#' @param d vector with binary treatment indicator, if not NULL,
#'          weights are normalized to sum to N in treated and control also in CV samples
#' @param parallel If TRUE, cross-validation parallelized
#' @param graph_file Option to define path to save graph
#' @param se_rule If not NULL, define, e.g., c(-1,1) to get 1SE and 1SE+ rule
#' @param ... Pass \code{\link{glmnet}} options
#' @import glmnet rms doParallel doSNOW
#'
#' @return List with the names of selected variables at cross-validated minima for Lasso and Post-Lasso
#' @export

post_lasso_cv <- function(x,y,w=NULL,kf = 10,family="gaussian", seed=NULL,output = TRUE,
                          d=NULL,parallel=FALSE,graph_file=NULL,se_rule = NULL,...) {

  # Create weights of ones if no weights are specified
  if (is.null(w)) {

    w <- as.matrix(rep(1,nrow(x)),nrow = nrow(x),ncol = 1)
  } else {
    w <- as.matrix(w,nrow = nrow(y),ncol = 1)
  }
  colnames(w) <- "w"

  # Normalize the weights either to N for all, or to N in treated and control group if d is specified
  if (is.null(d)) {
    w <- norm_w_to_n(w)
  } else {
    w <- norm_w_to_n(w,d)
  }

  # Lasso with full estimation sample
  lasso_full <- glmnet(x,y,weights = as.vector(w),family=family,...)
  if (output ==  TRUE) {
    print(lasso_full)
    if (!is.null(graph_file)) png(file = graph_file, units="in", width=7, height=10, res=1000,bg = "transparent")
    par(mfrow = c(2,1))
    plot(lasso_full, xvar = "lambda",label=TRUE)
  }

  coef_lasso_full <- coef(lasso_full)                   # Save coefficients to extract later the ones at the CV minima
  # print(coef_lasso_full)
  nm_coef_lasso_full <- rownames(coef_lasso_full)[-1]   # Save variable names that were used and kick out intercept
  lambda <- lasso_full$lambda                           # Save grid to use the same in cross-validation


  ###################################################
  ### Cross-validation with Lasso and Post Lasso ####
  ###################################################

  ## Figure out variables that were alway inactive to drop them in CV
  # Figure out which coefficients are inactive at each Lasso grid
  inact_coef_full <- (coef_lasso_full == 0)         # boolean matrix
  # Figure out number of active coefficients for graph
  num_act <- Matrix::colSums(1*!inact_coef_full)-1

  # Initialize matrices for MSE of Lasso and post-lasso for each grid point and CV fold
  cv_MSE_lasso <- matrix(nrow = kf,ncol = length(lambda))
  cv_MSE_post_lasso <- matrix(nrow = kf,ncol = length(lambda))

  # Get indicator for CV samples
  if (!is.null(seed)) {                                  # Use seed if specified
    set.seed(seed)
  }
  split = runif(nrow(x))
  cvgroup = as.numeric(cut(split,quantile(split, probs = seq(0,1,1/kf)),include.lowest = TRUE))  # groups for K-fold cv
  list <- 1:kf                            # Needed later in the loop to get the appropriate samples
  # table(cvgroup)

  ###########################################
  ########### CV NOT parallel ###############
  ###########################################
  if (parallel==FALSE) {

    ## Creating progress bar to get status of CV
    if (output ==  TRUE) {

      print(paste("Progress of",toString(kf),"fold cross-validation:"))
      progress.bar <- create_progress_bar("text")
      progress.bar$init(kf)
    }

    ## Start loop for cross-validation of Lasso and Post-Lasso

    for (i in 1:kf) {
      ## core CV program functionsxxx no. 9
      CV <- CV_core(x,y,w,d,cvgroup,list,i,lambda,family,...)

      ## Extract MSE of Lasso and Post-Lasso
      cv_MSE_lasso[i,] <- CV$MSE_lasso
      cv_MSE_post_lasso[i,] <- CV$MSE_post_lasso

      ## Update progress bar
      if (output ==  TRUE) progress.bar$step()

    }       # end loop over folds


    ###########################################
    ############### CV parallel ###############
    ###########################################

  }  else if (parallel == TRUE)  {

    ## Set up cluster of cores
    n_cores <- min(detectCores(),kf)

    if (output ==  TRUE) cat("\n CV parallelized using",toString(n_cores),"cores \n\n")
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)

    ## Parallelized foreach loop
    para_res <- foreach (i = 1:kf, .combine="rbind", .multicombine=TRUE, .packages=c("psych","glmnet","tcltk","RandomFieldsUtils","rms")
                         ,.export=c("CV_core","norm_w_to_n","fitted_values","add_intercept"),.inorder=FALSE) %dopar% {

                           if (output ==  TRUE) {
                             if(!exists("pb")) pb <- tkProgressBar(title="CV progress", min=1, max=kf)
                           }

                           ## core CV program functionsxxx no. 9
                           CV <- CV_core(x,y,w,d,cvgroup,list,i,lambda,family,...)

                           ## Extract MSE of Lasso and Post-Lasso
                           cv_MSE_lasso <- CV$MSE_lasso
                           cv_MSE_post_lasso <- CV$MSE_post_lasso

                           if (output ==  TRUE) setTkProgressBar(pb, i)

                           ## Matrices to be returned from cores
                           return(list(as.matrix(cv_MSE_lasso),as.matrix(cv_MSE_post_lasso)))
                         }     # end of parallel loop

    ## Shut down cluster
    stopCluster(cl)

    ## Extract results from list that is produced by foreach
    para_res <- as.matrix(do.call(cbind,para_res))
    cv_MSE_lasso <- t(para_res[,1:kf])
    cv_MSE_post_lasso <- t(para_res[,(kf+1):(2*kf)])
  }

  ## Calculate mean MSEs over all folds
  cv_MSE_lasso[is.na(cv_MSE_lasso)] <- max(cv_MSE_lasso) # can happen if glmnet does not go over the full grid
  cv_MSE_post_lasso[is.na(cv_MSE_post_lasso)] <- max(cv_MSE_post_lasso) # or Post-Lasso crashed
  mean_MSE_lasso <- colMeans(cv_MSE_lasso)
  mean_MSE_post_lasso <- colMeans(cv_MSE_post_lasso)
  mean_MSE_lasso[is.na(mean_MSE_lasso)] <- max(mean_MSE_lasso,na.rm=T) + 1e-7 # can happen if glmnet does not go over the full grid
  mean_MSE_post_lasso[is.na(mean_MSE_post_lasso)] <- max(mean_MSE_post_lasso,na.rm=T) + 1e-7 # or Post-Lasso crashed
  # Standard error of folds
  oneSE_lasso <- sqrt(apply(cv_MSE_lasso, 2, var)/kf)
  oneSE_post_lasso <- sqrt(apply(cv_MSE_post_lasso, 2, var)/kf)
  oneSE_lasso[is.na(oneSE_lasso)] <- 0          # otherwise the plot crashes
  oneSE_post_lasso[is.na(oneSE_post_lasso)] <- 0


  ## Get grid position of minimum MSE
  ind_min_l <- which.min(mean_MSE_lasso)
  ind_min_pl <- which.min(mean_MSE_post_lasso)
  ## Get names at minima
  names_l <- names(coef_lasso_full[,ind_min_l])[which(coef_lasso_full[,ind_min_l] != 0)]
  names_pl <- names(coef_lasso_full[,ind_min_pl])[which(coef_lasso_full[,ind_min_pl] != 0)]


  # Initialize lists for SE rules
  ind_Xse_l <- NULL
  nm_Xse_l <- NULL
  ind_Xse_pl <- NULL
  nm_Xse_pl <- NULL

  if (is.numeric(se_rule)) {
    ind_Xse_l <- rep(NA,length(se_rule))
    nm_Xse_l <- vector("list", length(se_rule))
    names(nm_Xse_l) <- paste("SE",se_rule)
    ind_Xse_pl <- rep(NA,length(se_rule))
    nm_Xse_pl <- vector("list", length(se_rule))
    names(nm_Xse_pl) <- paste("SE",se_rule)
    for (s in 1:length(se_rule)) {
      ind_Xse_l[s] <- find_Xse_ind(mean_MSE_lasso,ind_min_l,oneSE_lasso,se_rule[s])
      ind_Xse_pl[s] <- find_Xse_ind(mean_MSE_post_lasso,ind_min_pl,oneSE_post_lasso,se_rule[s])
      nm_Xse_l[[s]] <- names(coef_lasso_full[,ind_Xse_l[s]])[which(coef_lasso_full[,ind_Xse_l[s]] != 0)]
      nm_Xse_pl[[s]] <- names(coef_lasso_full[,ind_Xse_pl[s]])[which(coef_lasso_full[,ind_Xse_pl[s]] != 0)]
    }
  }

  ## Output
  if (output ==  TRUE) {

    # Calculate 1SE bands
    lasso_1se_up <- mean_MSE_lasso+oneSE_lasso
    lasso_1se_low <- mean_MSE_lasso-oneSE_lasso
    post_lasso_1se_up <- mean_MSE_post_lasso+oneSE_post_lasso
    post_lasso_1se_low <- mean_MSE_post_lasso-oneSE_post_lasso

    # Get ranges for the graph
    xrange <- range(log(lambda))
    yrange <- c(max(-1.7e+308,min(lasso_1se_low,post_lasso_1se_low)),
                min(1.7e+308,max(lasso_1se_up,post_lasso_1se_up)))

    # Plot mean lines
    ylab <- "Mean-squared Error"
    plot(xrange,yrange,type="n",xlab="Log Lambda",ylab=ylab)
    lines(log(lambda),mean_MSE_lasso,lwd=1.5,col="blue")
    lines(log(lambda),mean_MSE_post_lasso,lwd=1.5,col="red")

    # Plot upper and lower 1SE lines
    lines(log(lambda),lasso_1se_up,lty=2,lwd=1,col="blue")
    lines(log(lambda),lasso_1se_low,lty=2,lwd=1,col="blue")
    lines(log(lambda),post_lasso_1se_up,lty=2,lwd=1,col="red")
    lines(log(lambda),post_lasso_1se_low,lty=2,lwd=1,col="red")

    # Show location of minima and 1SE
    abline(v = log(lambda[ind_min_l]), lty = 1, col="blue")
    abline(v = log(lambda[ind_min_pl]), lty = 1, col="red")

    if (is.numeric(se_rule)) {
      for (s in 1:length(se_rule)) {
        abline(v = log(lambda[ind_Xse_l[s]]), lty = 3, col="blue")
        abline(v = log(lambda[ind_Xse_pl[s]]), lty = 3, col="red")
      }
    }

    # Print legend
    legend('top',c("Lasso MSE","Lasso MSE?1SE","Post-Lasso MSE","Post-Lasso MSE?1SE","# active coeff."), lty = c(1,2,1,2,1),
           col=c('blue','blue','red','red','forestgreen'),ncol=1,bty ="n",cex=0.7)

    # Open a new graph for number of coefficients to be written on existing
    par(new = TRUE)
    plot(log(lambda),num_act, axes=F, xlab=NA, ylab=NA, cex=1.2,type="l",col="forestgreen",lwd=1.5)
    axis(side = 4)
    mtext(side = 4, line = 3, "# active coefficients")

    if (!is.null(graph_file)) dev.off()

    # Comparison of minimum MSE
    cat("\n\n Minimum ",ylab," Lasso:",toString(min(mean_MSE_lasso,na.rm = TRUE)))
    cat("\n\n Minimum ",ylab," Post-Lasso:",toString(min(mean_MSE_post_lasso,na.rm = TRUE)))

    # Show names of active variables at respective minima
    cat("\n\n Active variables at minimum of Lasso: \n")
    print(names(coef_lasso_full[,ind_min_l])[which(coef_lasso_full[,ind_min_l] != 0)])

    cat("\n\n Active variables at minimum of Post-Lasso: \n")
    print(names(coef_lasso_full[,ind_min_pl])[which(coef_lasso_full[,ind_min_pl] != 0)])
  }

  ## Get Lasso coefficients at minimum of Lasso
  coef_min_l <- coef_lasso_full[,ind_min_l][which(coef_lasso_full[,ind_min_l] != 0)]
  ## Get Lasso coefficients at minimum of Post-Lasso
  coef_min_pl <- coef_lasso_full[,ind_min_pl][which(coef_lasso_full[,ind_min_pl] != 0)]

  ## Return names and coefficients
  list("names_l" = names_l,"names_pl" = names_pl,"names_Xse_l"=nm_Xse_l
                              ,"names_Xse_pl"=nm_Xse_pl
                              ,"coef_min_l" = coef_min_l,"coef_min_pl" = coef_min_pl)
}


#' This function extracts a subset of active variables (nm_act) of the
#' relavant variables from X'X and X'y to get out-of-sample predictions
#' for a matrix containing only the active variables
#'This speeds up the cross-validation for post-Lasso to a large extent
#'
#' @param XtX_all Crossproduct of all covariates
#' @param Xty_all Crossproduct of covariates and outcome
#' @param x_pred Covariates matrix of the PREDICTION sample
#' @param nm_act names of active variables
#'
#' @return fit_val fitted values in the prediction sample


fitted_values <- function (XtX_all,Xty_all,x_pred,nm_act) {
  # if (!require("RandomFieldsUtils",quietly=TRUE)) {
  #   stop("Package RandomFieldsUtils not installed (required here)!")
  # }

  # Extract relevant rows and columns
  XtX <- XtX_all[nm_act,nm_act]
  Xty <- Xty_all[nm_act,]

  # Calculate point estimates in estimation sample
  b <- tryCatch(solvex(XtX, Xty), error=function(e) NULL) # much faster than standard solve
  if (is.null(b)) {        # In case b has not been estimated
    # Option 1: replace by mean
    # print("X'X not invertible: Prediction replaced by mean.")
    # b <- solvex(XtX[1,1], Xty[1])
    # Fitted value only the means
    # fit_val <- as.matrix(x_pred[,1])%*%b
    # Option 2: Return NAs
    # fit_val <- as.matrix(rep(NA,nrow(x_pred)),ncol=1)
    # Option 3: Pass NULL value
    fit_val <- b
  } else {
    # Calculate fitted values in prediction sample
    fit_val <- x_pred%*%b
  }
  return(fit_val)
}


#### Function to normalize weights to N or to N in treated and controls separately

norm_w_to_n <- function(w,d=NULL) {
  ## Input:  w: vector of weights that should be normalized
  ##         d: treatment indicator
  ## Output: w: Normalized weights
  if (is.null(d)) {
    ## Normalize weights to sum up to N
    w <- w / sum(w)* nrow(w)
  } else {
    # Separate weights of treated and controls
    w1 <- w * d
    w0 <- w * (1-d)
    # Normalize weights to sum to N in both groups
    w1 <- w1 / sum(w1) * nrow(w)
    w0 <- w0 / sum(w0) * nrow(w)
    # Unify weights again
    w <- w1 + w0
  }
  return(w)
}


#' This function contains the core parts of the CV for Lasso and Post-Lasso

#' @param x covariate matrix to be used in CV
#' @param y: vector of outcomes
#' @param w: vector of weight
#' @param d: vector of treament indicators
#' @param cvgroup: categorical with k groups to identify folds
#' @param list: list 1:k
#' @param i: number of fold that is used for prediction
#' @param lambda: series of lambdas used
#' @param ... Pass \code{\link{glmnet}} options
#'
#' @return MSE_lasso / MSE_post_lasso: means squared errors for each lambda

CV_core <- function(x,y,w,d,cvgroup,list,i,lambda,family,...) {

  # Get estimation and prediction sample for this specific fold
  x_est_cv <- subset(x,cvgroup %in% list[-i])
  y_est_cv <- subset(y,cvgroup %in% list[-i])
  w_est_cv <- subset(w,cvgroup %in% list[-i])
  x_pred_cv <- subset(x,cvgroup %in% list[i])
  y_pred_cv <- subset(y,cvgroup %in% list[i])
  w_pred_cv <- subset(w,cvgroup %in% list[i])

  # Normalize the weights either to N for all or to N in treated and control group if d is specified
  if (is.null(d)) {
    w_est_cv <- norm_w_to_n(w_est_cv)
    w_pred_cv <- norm_w_to_n(w_pred_cv)
  } else {
    d_est_cv <- subset(d,cvgroup %in% list[-i])
    d_pred_cv <- subset(d,cvgroup %in% list[i])
    w_est_cv <- norm_w_to_n(w_est_cv,d_est_cv)
    w_pred_cv <- (norm_w_to_n(w_pred_cv,d_pred_cv) / 2)    # /2 to normalize to N again for MSE, not necessary
  }

  # Estimate Lasso for this fold using the grid of the full sample
  lasso_cv <- glmnet(x_est_cv, y_est_cv,lambda = lambda,weights=as.vector(w_est_cv),
                     family=family,...)
  coef_lasso_cv <- coef(lasso_cv)                                       # Save coefficients at each grid point


  # Predicted values with lasso coefficients for prediction sample at each grid
  fit_lasso <- predict(lasso_cv,x_pred_cv)
  if (family == "binomial") fit_lasso <- predict(lasso_cv,newx = x_pred_cv, type="response")

  if (ncol(fit_lasso) != length(lambda)) {
    fit_lasso <- cbind(fit_lasso,matrix(NA,nrow(fit_lasso),(length(lambda)-ncol(fit_lasso))))
  }

  ## Now calculate predicted values for post Lasso at each grid
  # Initialize first matrix for fitted values
  fit_post_lasso <- matrix(NA,nrow = nrow(fit_lasso),ncol = ncol(fit_lasso))

  # Figure out which coefficients are active at each Lasso grid
  act_coef <- (coef_lasso_cv != 0)         # boolean matrix

  ## Calculate full covariance matrix X'X and X'y once and select only the relevant in the loop below (much faster)
  # First, figure out variables that were active at least once and get only these
  act_once <- apply(act_coef,1,any)
  nm_all_act_coef <- rownames(act_coef)[act_once]
  if (length(nm_all_act_coef) == 1) {
    warning("No variables selected in one CV, even for lowest lambda, might reconsider choice of penalty term")
    fit_post_lasso <- fit_lasso
  } else {
    ## Create the covariate matrix to "manually" calculate the fitted values (much faster than the build in lm command)
    if (sum(abs(coef_lasso_cv[1,]))==0) {   # Indicates that no intercept was used
      x_all_act <- x_est_cv[,nm_all_act_coef]
    }
    else if (sum(abs(coef_lasso_cv[1,]))!=0) {    # Indicates that intercept was used
      x_all_act <- add_intercept(x_est_cv[,nm_all_act_coef[2:length(nm_all_act_coef)],drop=F])
      colnames(x_all_act)[1] <- "(Intercept)"
      # add intercept also to prediction sample
      x_pred_cv <- add_intercept(x_pred_cv[,nm_all_act_coef[2:length(nm_all_act_coef)],drop=F])
      colnames(x_pred_cv)[1] <- "(Intercept)"
      # print(nm_all_act_coef)
    }
    else {
      stop("Something strange happens with the intercepts at the Lasso path")
    }

    if (family == "gaussian") {
      # Add weights
      x_w <- apply(x_all_act,2,`*`,sqrt(w_est_cv))
      y_w <- y_est_cv * sqrt(w_est_cv)

      # Get X'X
      XtX_all <- crossprod(x_w)
      # Get X'y
      Xty_all <- crossprod(x_w,y_w)
    }


    # Initialize vector to save chosen variable names to figure out whether new variables were added from the last to the next grid
    nm_act_coef_prev = "random variable name, nobody should call a variable like this"

    ## Loop over the grid of Lambdas to get the Post-Lasso predictions
    for (j in 1:length(lambda)) {

      # Get names of active variables at this grid
      nm_act_coef <- rownames(act_coef)[act_coef[,j]]
      if (identical(nm_act_coef,character(0))==TRUE) {
        # print("No variable selected at this grid => no prediction")
        next
      }
      if (identical(nm_act_coef,nm_act_coef_prev) == TRUE & j>1) {
        # print("Same variables selected as in previous grid => Post-Lasso predictions remain unchanged")
        fit_post_lasso[,j] = fit_post_lasso[,(j-1)]
        next
      }
      else {     # Get prediction covariate matrix for that grid
        x_ols_pred <- as.matrix(x_pred_cv[,nm_act_coef])
      }

      # Get OLS Post-Lasso predictions for this grid point
      if (family == "gaussian") {
        fit <- fitted_values(XtX_all,Xty_all,x_ols_pred,nm_act_coef)
        if (is.null(fit) & j==1) {
          fit_post_lasso[,j] = rep(mean(y == 1),nrow(fit_post_lasso))
          next
        }
        if (is.null(fit) & j>1) {
          # cat("\n X'X not invertible at grid",toString(j),": Use last feasible coefficients")
          fit_post_lasso[,j] = fit_post_lasso[,(j-1)]
          next
        } else{
          fit_post_lasso[,j] <- fit
        }
      } else if (family == "binomial") {
        # Check if only intercept was selected and predict mean in this case
        if (length(nm_act_coef) == 1) {
          fit_post_lasso[,j] <- rep(mean(y_est_cv == 1),nrow(fit_post_lasso))
          next
        }
        # Otherwise run unrestricted logit
        # Get matrix
        x_logit <- x_est_cv[,nm_act_coef[-1]]
        # Check dimnesion and rank if more than one variable selected
        if (is.null(dim(x_logit))) {
          x_logit <- as.matrix(x_logit,ncol=1)
        } else {
          # Use last feasible predictions in case of collinearity
          if (qr(apply(x_logit,2,`*`,sqrt(w_est_cv)), tol = 1e-06)$rank < ncol(x_logit) & j==1) {
            fit_post_lasso[,j] = rep(mean(y_est_cv == 1),nrow(fit_post_lasso))
            next
          }
          if (qr(apply(x_logit,2,`*`,sqrt(w_est_cv)), tol = 1e-06)$rank < ncol(x_logit) & j>1) {
            fit_post_lasso[,j] = fit_post_lasso[,(j-1)]
            next
          }
        }

        # If ok run logit and get coefficients
        coef_logit <- lrm.fit(x_logit,y_est_cv,weights = w_est_cv,tol=1e-30)$coefficients

        if (is.null(coef_logit) & j==1) {
          fit_post_lasso[,j] = rep(mean(y_est_cv == 1),nrow(fit_post_lasso))
          next
        }
        if (is.null(coef_logit) & j>1) {
          fit_post_lasso[,j] = fit_post_lasso[,(j-1)]
          next
        }

        # Get linear index and prediction
        li <- x_ols_pred %*% as.matrix(coef_logit,ncol=1)
        fit_post_lasso[,j] <- exp(li) / (1 + exp(li))

      } else {
        stop(cat("\n Post-Lasso for this outcome family not implemented \n"))
      }


      # Update current active covariates for the check whether it changes in the next grid
      nm_act_coef_prev <- nm_act_coef

    }       # end loop over grids

  } # end if at least one var selected


  # Matrix with "real" outcome values for each grid
  y_rep <- matrix(rep(y_pred_cv,length(lambda)),nrow = nrow(fit_lasso),ncol = ncol(fit_lasso))

  # Get RMSE
  SE_lasso <- (y_rep - fit_lasso)^2
  SE_post_lasso <- (y_rep - fit_post_lasso)^2
  if (!is.null(w)) {                                     # Weight errors by "sampling weights"
    SE_lasso <- apply(SE_lasso,2,"*",w_pred_cv)
    SE_post_lasso <- apply(SE_post_lasso,2,"*",w_pred_cv)
  }
  MSE_lasso <- apply(SE_lasso,2,mean)
  MSE_post_lasso <- apply(SE_post_lasso,2,mean)

  list("MSE_lasso" = MSE_lasso,"MSE_post_lasso" = MSE_post_lasso)
  # if(family == "gaussian")  list("MSE_lasso" = MSE_lasso,"MSE_post_lasso" = MSE_post_lasso)
  # else list("MSE_lasso" = MSE_lasso,"MSE_post_lasso" = MSE_post_lasso,
  #           "BD_lasso"=BD_lasso,"BD_post_lasso"=BD_post_lasso)
}



#' Function finds the positionfor pre-specified SE rules

find_Xse_ind <- function(CV,ind_min,oneSE,factor) {
  cv_temp <- CV - (CV[ind_min] + abs(factor) * oneSE[ind_min])
  if (factor < 0) {
    for (i in ind_min:1) {
      ind <- i
      if (cv_temp[i] < 0) next
      else if (cv_temp[i] > 0) break
    }
  } else {
    for (i in ind_min:length(oneSE)) {
      ind <- i
      if (cv_temp[i] < 0) next
      else if (cv_temp[i] > 0) break
    }
  }
  return(ind)
}


#' Adds an intercept to a matrix

add_intercept <- function(mat) {
  if (is.null(dim(mat))) mat <- as.matrix(mat,ncol=1)
  mat <- cbind(rep(1,nrow(mat)),mat)
  colnames(mat)[1] <- "(Intercept)"
  return(mat)
}
