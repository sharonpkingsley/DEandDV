calculateDEandDV<- function(x, labels){
  input<-x

  # check labels matches matrix columns
  if(length(labels)!=ncol(x)) {
    stop("labels must be same size as number of samples");
  }
  # check valid labels
  check<-(labels%in%c(0,1))
  if(!all(check==TRUE)){
    stop("labels must contain values 0 or 1");
  }

  ##iniatialize matrices for the 2 groups of cells
  exprs1<-matrix(nrow=nrow(x), ncol = 0)
  exprs2<-matrix(nrow=nrow(x), ncol = 0)

  ##Populate matrices based on labels
  for(n in 1:length(labels)){
    if(labels[n]==0){
      exprs1<-cbind(exprs1,x[,n])
    }
    if(labels[n]==1){
      exprs2<-cbind(exprs2,x[,n])
    }

  }

  exprs<-exprs1


  #Normalization - Log counts per million
  cpm <- apply(exprs,2, function(x) (x/sum(x))*1000000)
  log.cpm <- log(cpm + 1)
  exprs<-log.cpm


  #Basic Values - No.of Dropout candidates, mean of expressed genes
  dropouts <- rowSums(exprs==0)
  dropout_proportion <- dropouts/ncol(exprs)
  row_sums <- rowSums(exprs)
  non_zeros <- ncol(exprs)-dropouts
  means <- numeric()
  for(row in 1:nrow(exprs)){
    if(non_zeros[row]==0){
      means[row] = 0;
    } else {
      means[row]<-row_sums[row]/non_zeros[row]
    }
  }

  #Model dropout proportions vs means
  model <- glm (dropout_proportion~means, na.action = na.exclude)
  #Predict probabily of dropout for each gene
  p_x <- predict(model)

  #Statistically determine the expected expression for each gene
  d_mean <- (1-p_x) * means
  x_mean <- means
  m <- non_zeros
  n <- dropouts
  s_mean<- ((m*d_mean)+(n*x_mean))/(n+m)

  #Statistically determine the variance for each gene
  d_sd <- (1-p_x) * means * means
  standard_deviations <- numeric()
  for(row in 1:nrow(exprs)){
    exprs_1<-as.matrix(exprs[row,])
    exprs_2<-as.matrix(exprs_1[exprs_1!=0])
    standard_deviations[row]<-sd(exprs_2)
    if(is.na(standard_deviations[row])){
      standard_deviations[row]=0
    }
  }
  x_sd<-standard_deviations
  s_sd <- (((n*d_sd) +(m*x_sd))/(m+n)) + (((n*(d_mean-s_mean)*(d_mean-s_mean)) + (m*(x_mean-s_mean)*(x_mean-s_mean)))/(m+n))

  #Remove outlier function
  remove_outliers <- function(x, na.rm = TRUE, ...) {
    qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
    H <- 1.5 * IQR(x, na.rm = na.rm)
    y <- x
    y[x < (qnt[1] - H)] <- NA
    y[x > (qnt[2] + H)] <- NA
    y
  }

  s_mean_2<-remove_outliers(s_mean)
  s_sd_2<-remove_outliers(s_sd)

  #Store mean, variance and number of samples for first group of cells
  m1<-s_mean
  s1<-s_sd
  n1<-ncol(exprs)



  #Repeat correction within second cell group
  exprs<-exprs2

  #Normalization - Log counts per million
  cpm <- apply(exprs,2, function(x) (x/sum(x))*1000000)
  log.cpm <- log(cpm + 1)
  exprs<-log.cpm

  #Basic Values - No.of Dropout candidates, mean of expressed genes
  dropouts <- rowSums(exprs==0)
  dropout_proportion <- dropouts/ncol(exprs)
  row_sums <- rowSums(exprs)
  non_zeros <- ncol(exprs)-dropouts
  means <- numeric()
  for(row in 1:nrow(exprs)){
    if(non_zeros[row]==0){
      means[row] = 0;
    } else {
      means[row]<-row_sums[row]/non_zeros[row]
    }
  }

  #Model dropout proportions vs means
  model <- glm (dropout_proportion~means, na.action = na.exclude)
  #Predict probabily of dropout for each gene
  p_x <- predict(model)

  #Statistically determine the expected expression for each gene
  d_mean <- (1-p_x) * means
  x_mean <- means
  m <- non_zeros
  n <- dropouts
  s_mean<- ((m*d_mean)+(n*x_mean))/(n+m)

  #Statistically determine the variance for each gene
  d_sd <- (1-p_x) * means * means
  standard_deviations <- numeric()
  for(row in 1:nrow(exprs)){
    exprs_1<-as.matrix(exprs[row,])
    exprs_2<-as.matrix(exprs_1[exprs_1!=0])
    standard_deviations[row]<-sd(exprs_2)
    if(is.na(standard_deviations[row])){
      standard_deviations[row]=0
    }
  }

  x_sd<-standard_deviations
  s_sd <- (((n*d_sd) +(m*x_sd))/(m+n)) + (((n*(d_mean-s_mean)*(d_mean-s_mean)) + (m*(x_mean-s_mean)*(x_mean-s_mean)))/(m+n))

  s_mean_2<-remove_outliers(s_mean)
  s_sd_2<-remove_outliers(s_sd)

  #Store mean, variance and number of samples for second group of cells
  m2<-s_mean
  s2<-s_sd
  n2<-ncol(exprs)

  #variance to standard deviation
  s1<-sqrt(s1)
  s2<-sqrt(s2)

  #T-test function
  t.test2 <- function(m1,m2,s1,s2,n1,n2)
  {
    se <- sqrt( (s1^2/n1) + (s2^2/n2) )
    # welch-satterthwaite df
    df <- ( (s1^2/n1 + s2^2/n2)^2 )/( (s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1) )

    t <- (m1-m2)/se
    dat <- c(m1-m2, se, t, 2*pt(-abs(t),df))
    names(dat) <- c("Difference of means", "Std Error", "t", "p-value")
    return(dat)
  }

  #initialize matrix to be returned
  toreturn<- matrix(nrow=nrow(x), ncol=7)
  colnames(toreturn)<- c("LFC", "t", "p", "p.adj", "F", "F.p", "F.p.adj")

  #Conduct DE analysis with a t-test
  #Conduct DV analysis with an f-test
  #Populate matrix with values
  for(row in 1:nrow(x)){
    tt<-t.test2(m1[row],m2[row],s1[row],s2[row], n1, n2)
    x<-rnorm(n1,m1[row],s1[row])
    y<-rnorm(n2,m2[row],s2[row])
    if((any(is.na(x)))||(any(is.na(y)))){

    }else{
      z<-var.test(x,y)
      toreturn[row,"F"]<-z$statistic
      toreturn[row,"F.p"]<-z$p.value
      toreturn[row,"F.p.adj"]<-p.adjust(z$p.value, method = "BH")
    }
    toreturn[row,"t"]<-tt["t"]
    toreturn[row,"p"]<-tt["p-value"]
    toreturn[row,"p.adj"]<-p.adjust(tt["p-value"], method = "BH")

    toreturn[row,"LFC"]<-m2[row]-m1[row]

  }

  rownames(toreturn)<-rownames(input)

  return(toreturn)

}
