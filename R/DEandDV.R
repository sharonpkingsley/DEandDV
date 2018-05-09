calculateDE <- function(x, plot=FALSE){
  exprs<-x
  
  dropouts <- rowSums(exprs==0)
  dropout_proportion <- dropouts/ncol(exprs)
  row_sums <- rowSums(exprs)
  non_zeros <- ncol(exprs)-dropouts
  means <- row_sums/non_zeros
  
  
  model <- glm (dropout_proportion~means)
  
  p_x <- predict(model)
  d_mean <- (1-p_x) * means
  x_mean <- means
  m <- non_zeros
  n <- dropouts
  s_mean<- ((m*d_mean)+(n*x_mean))/(n+m)
  
  mean_nofix <- rowMeans(exprs)
  
  remove_outliers <- function(x, na.rm = TRUE, ...) {
    qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
    H <- 1.5 * IQR(x, na.rm = na.rm)
    y <- x
    y[x < (qnt[1] - H)] <- NA
    y[x > (qnt[2] + H)] <- NA
    y
  }
  
  s_mean_2<-remove_outliers(s_mean)
  mean_nofix_2<-remove_outliers(mean_nofix)
  
  if(plot==TRUE){
    plot(s_mean_2, mean_nofix_2)   
  }
  
  
  return(s_mean_2)
  
}


calculateDV <- function(x, plot=FALSE){
  exprs<-x
  
  dropouts <- rowSums(exprs==0)
  dropout_proportion <- dropouts/ncol(exprs)
  row_sums <- rowSums(exprs)
  non_zeros <- ncol(exprs)-dropouts
  means <- row_sums/non_zeros
  
  
  model <- glm (dropout_proportion~means)
  
  p_x <- predict(model)
  d_mean <- (1-p_x) * means
  x_mean <- means
  m <- non_zeros
  n <- dropouts
  s_mean<- ((m*d_mean)+(n*x_mean))/(n+m)
  
  d_sd <- (1-p_x) * means * means
  standard_deviations <- numeric()
  for(row in 1:nrow(exprs)){
    exprs1<-as.matrix(exprs[row,])
    exprs2<-as.matrix(exprs1[exprs1!=0])
    standard_deviations[row]<-sd(exprs2)
  }
  
  x_sd<-standard_deviations
  
  s_sd <- (((n*d_sd) +(m*x_sd))/(m+n)) + (((n*(d_mean-s_mean)*(d_mean-s_mean)) + (m*(x_mean-s_mean)*(x_mean-s_mean)))/(m+n))
  sd_nofix <- apply(exprs,1,sd)
  
  
  remove_outliers <- function(x, na.rm = TRUE, ...) {
    qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
    H <- 1.5 * IQR(x, na.rm = na.rm)
    y <- x
    y[x < (qnt[1] - H)] <- NA
    y[x > (qnt[2] + H)] <- NA
    y
  }
  
  s_sd_2<-remove_outliers(s_sd)
  sd_nofix_2<-remove_outliers(sd_nofix)
  
  
  
  
  if(plot==TRUE){
    plot(s_sd_2, sd_nofix_2)
  }
  
  
  return(s_mean_2)
  
}
