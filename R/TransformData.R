#' Performs Box-cox Transformation on a vector
#'
#' Takes in a vector and transforms it using a Box-cox Power Trasformation
#' @param one.col Any numeric vector that needs to be transformed
#' @param returnLambda A boolean value which returns the optimal lambda value for a vector when se to TRUE
#' @return Returns the optimal lambda value for Box-cox transformation when "returnLambda" is set to TRUE and a Box-cox transformed vector otherwise.
#' @export
#' @usage boxcoxtransForecast(one.col, returnLambda = FALSE)
#' @example
#' a <- c(12,34,234,23,678, 768, 34, 34 ,78)
#' boxcoxtransForecast(a)
#'

boxcoxtransForecast <- function(one.col, returnLambda = FALSE)
{
  # require(forecast)
  if (returnLambda == TRUE)
    return(BoxCox.lambda(one.col))
  return(BoxCox(one.col, BoxCox.lambda(one.col)))
}

#' Performs Yeo-Johnson Transformation on a vector
#'
#' Takes in a vector and transforms it using a Yeo-Johnson Power Trasformation
#' @param one.col Any numeric vector that needs to be transformed
#' @param returnLambda A boolean value which returns the optimal lambda value for a vector when se to TRUE
#' @return Returns the optimal lambda value for Yeo-Johnson transformation when "returnLambda" is set to TRUE and a Yeo-Johnson transformed vector otherwise.
#' @export
#' @usage yeoJohnsonVGAM(one.col, returnLambda = FALSE)
#' @example
#' a <- c(12,34,234,23,678, 768, 34, 34 ,78)
#' yeoJohnsonVGAM(a)
#'

yeoJohnsonVGAM <- function(one.col, returnLambda = FALSE)
{
  # require(forecast)
  # require(VGAM)
  if (returnLambda == TRUE)
    return(BoxCox.lambda(one.col))
  return(yeo.johnson(one.col, BoxCox.lambda(one.col)))
}


#' Performs Lambert Transformation on a vector
#'
#' Takes in a vector and transforms it using a Lmabert-W Trasformation
#' @param one.col Any numeric vector that needs to be transformed
#' @param return_tau_mat A boolean value which returns the tau_mat vector when set to TRUE
#' @return Returns the tau_mat vector for Lambert-W transformation when "return_tau_mat" is set to TRUE and a Lambert-W transformed vector otherwise.
#' @export
#' @usage GaussianizeLambertW(one.col, return_tau_mat = FALSE)
#' @example
#' a <- c(12,34,234,23,678, 768, 34, 34 ,78)
#' GaussianizeLambertW(a, return_tau_mat = FALSE)
#'

GaussianizeLambertW <- function(one.col, return_tau_mat = FALSE)
{
  # require(LambertW)
  if ((sort(table(one.col), decreasing = T)[1]/length(one.col) * 100) < 50)
  {
    data.col <- data.frame(Gaussianize(one.col))
    tau_mat <- Gaussianize(one.col, return.tau.mat = TRUE)$tau.mat
  } else
  {
    data.col <- one.col
    tau_mat <- NULL
  }
  if (return_tau_mat == TRUE)
    return(tau_mat)
  return(data.col)
}


#' Transforms the dataset using the best possible transformation technique
#'
#' Takes in a data frame and performs the best possible transformation to each of the columns in the data frame
#' @param data Any data frame that has atlest one column with continuous data and that has to transformed into normal form
#' @param ChooseBestTrans Considers the best transformation based on one of the three values ("Shapiro P Value", "Pearson P Value", "Min skewness")
#' @description
#' VariableTransform initially pads the data inorder to eliminates all negative and zero values and then categorizes the data into normal, positive skewed and negative skewed based on the skewness score of each column of the given data frame
#'
#' \strong{Transformation for Positive skewed data:}
#'
#' Different transformation techiniques considered to transform positive skewed data are as follows:
#' \itemize{
#' \item Log Transformation
#' \item Square Root Transformation
#' \item Cube Root Transformation
#' \item Turkey Power Transformation
#' \item Box-Cox Power Transformation
#' \item Yeo-Johnson Transformation
#' \item Order-Norm Transformation
#' \item Lambert-W Transformation
#' }
#' Of all these techniques, Best technique is choosed based on the "ChooseBestTrans" argument provided.Default method will be Shapiro P-Value
#'
#' \strong{Transformation for Negative skewed data:}
#'
#' Different transformation techiniques considered to transform negative skewed data are as follows:
#' \itemize{
#' \item Square Transformation
#' \item Cube Transformation
#' \item Turkey Power Transformation
#' \item Box-Cox Power Transformation
#' \item Yeo-Johnson Transformation
#' \item Order-Norm Transformation
#' \item Lambert-W Transformation

#' }
#' Of all these techniques, Best technique is choosed based on the "ChooseBestTrans" argument provided.Default method will be Shapiro P-Value
#'
#' \strong{Scaling the dataset:}
#'
#' After merging all the datasets created (normally distributed data, positive skewd data, neggative skewed data), entire dataset is normalized and stored in an other dataframe.
#'
#' @return
#' Returns a list of 7 objects:
#' \describe{
#' \item{transformed_df}{Tranformed Dataset for all Continuous variables}
#' \item{scaled_df}{Scaled Dataset for all Continuous variables}
#' \item{original_dist}{Distribution of the dataset provided}
#' \item{neg_skew_trans}{Intermediate scores (for all the techniques) for Positive skewed data}
#' \item{pos_skew_trans}{Intermediate scores (for all the techniques) for Positive skewed data}
#' \item{trans_fit_model}{A list of model fit file, OrderNorm objects and tau_mat objects}
#' \item{complete_data}{Combination of transformed continuous data and categorical data}
#' }
#'
#' @export
#' @usage VariableTransform(data, ChooseBestTrans = "Shapiro P Value")

VariableTransform<-function(data, ChooseBestTrans = "Shapiro P Value")
{
  print("Transforming & Scaling data..")
  # require(rcompanion)
  # require(plyr)
  # require(data.table)
  # require(bestNormalize)

  data_org<-data
  distribution<-data_distribution(data,dv="NULL")
  distribution$is_dv<-NULL
  distribution<-data.table(distribution)[distribution=="Continous",]

  #Padding for Negetive & 0 Values. Also Storing the minvalue for a Reverse pad
  minvalues<-data.table(data)[,lapply(.SD,function(x){min(x)}),.SDcols=distribution$names]
  minvalues<-data.frame(t(minvalues))
  colnames(minvalues)<-"min_value"
  minvalues$names<-rownames(minvalues)
  data<-data.table(data)[,lapply(.SD,function(x){x+1-min(x)}),.SDcols=distribution$names]

  #Re-doing Distribution to be sure of any padding effects - Ideally should not have changed
  distribution<-data_distribution(data,dv="NULL")
  distribution_orginal<-distribution #Saving for return
  distribution$is_dv<-NULL
  distribution<-data.table(distribution)[distribution=="Continous",]

  #Catagorizing Variables
  normal_var<-as.character(data.table(distribution)[is_normal=="TRUE"|(skewness_res>0&skewness_res<=0.3)|(skewness_res<0&skewness_res>=-0.3),]$names)
  negetive_skewed<-as.character(data.table(distribution)[(skewness_res<=-0.3),]$names)
  positive_skewed<-as.character(data.table(distribution)[(skewness_res>=0.3),]$names)


  #Saving Normal df For later
  normal_df<-data.frame(data)[,normal_var]

  #Addressing Positive/Right Skew (log,sqrt,cube root,turkeypower, boxcoxpower)
  log_transformed<-data.table(data)[,lapply(.SD,function(x){log(x)}),.SDcols=positive_skewed]
  sqrt_transformed<-data.table(data)[,lapply(.SD,function(x){sqrt(x)}),.SDcols=positive_skewed]
  cube_root_transformed<-data.table(data)[,lapply(.SD,function(x){sign(x)*abs(x)^(1/3)}),.SDcols=positive_skewed]
  turkey_transformed<-data.table(data)[,lapply(.SD,function(x){transformTukey(x,plotit=FALSE,quiet = T)}),.SDcols=positive_skewed]
  turkey_lamda<-data.table(data)[,lapply(.SD,function(x){transformTukey(x,plotit=FALSE,returnLambda = T,quiet = T)}),.SDcols=positive_skewed]
  turkey_lamda<-t(turkey_lamda)
  turkey_lamda<-data.frame(turkey_lamda)
  turkey_lamda$names<-rownames(turkey_lamda)

  boxcox_transformed<-data.table(data)[,lapply(.SD,function(x){boxcoxtransForecast(x)}),.SDcols=positive_skewed]
  boxcox_lamda<-data.table(data)[,lapply(.SD,function(x){boxcoxtransForecast(x,returnLambda = TRUE)}),.SDcols=positive_skewed]
  boxcox_lamda<-t(boxcox_lamda)
  boxcox_lamda<-data.frame(boxcox_lamda)
  boxcox_lamda$names<-rownames(boxcox_lamda)

  yeo_johnson_transformed<-data.table(data)[,lapply(.SD,function(x){yeoJohnsonVGAM(x)}),.SDcols=positive_skewed]
  yeo_johnson_lamda<-data.table(data)[,lapply(.SD,function(x){yeoJohnsonVGAM(x,returnLambda = TRUE)}),.SDcols=positive_skewed]
  yeo_johnson_lamda<-t(yeo_johnson_lamda)
  yeo_johnson_lamda<-data.frame(yeo_johnson_lamda)
  yeo_johnson_lamda$names<-rownames(yeo_johnson_lamda)

  order_norm_transformed <- data.table(data)[,lapply(.SD,function(x){bestNormalize::orderNorm(x, warn = FALSE)$x.t}),.SDcols=positive_skewed]
  order_norm_objects <- list()
  for (i in 1:length(positive_skewed))
  {
    order_norm_objects[[positive_skewed[i]]] <- bestNormalize::orderNorm(data.frame(data)[,positive_skewed[i]], warn = FALSE)
  }

  lambertw_transformed<-data.table(data)[,lapply(.SD,function(x){GaussianizeLambertW(x)}),.SDcols=positive_skewed]
  colnames(lambertw_transformed) <- positive_skewed
  lambertw_tau_mat <- list()
  for (i in 1:length(positive_skewed))
  {
    lambertw_tau_mat[[positive_skewed[i]]] <- GaussianizeLambertW(data.frame(data)[,positive_skewed[i]],return_tau_mat = TRUE)
  }

  log_dist<-data.frame(data_distribution(log_transformed,dv="NULL"))
  sqrt_dist<-data.frame(data_distribution(sqrt_transformed,dv="NULL"))
  cube_root_dist<-data.frame(data_distribution(cube_root_transformed,dv="NULL"))
  turkey_dist<-data.frame(data_distribution(turkey_transformed,dv="NULL"))
  boxcox_dist<-data.frame(data_distribution(boxcox_transformed,dv="NULL"))
  yeo_johnson_dist<-data.frame(data_distribution(yeo_johnson_transformed,dv="NULL"))
  order_norm_dist<-data.frame(data_distribution(order_norm_transformed,dv="NULL"))
  lambertw_dist<-data.frame(data_distribution(lambertw_transformed,dv="NULL"))

  log_dist$is_dv<-NULL
  sqrt_dist$is_dv<-NULL
  cube_root_dist$is_dv<-NULL
  turkey_dist$is_dv<-NULL
  boxcox_dist$is_dv<-NULL
  yeo_johnson_dist$is_dv<-NULL
  order_norm_dist$is_dv<-NULL
  lambertw_dist$is_dv<-NULL

  colnames(log_dist)[3:ncol(log_dist)]<-paste0("log_trans_",colnames(log_dist)[3:ncol(log_dist)])
  colnames(sqrt_dist)[3:ncol(sqrt_dist)]<-paste0("sqrt_trans_",colnames(sqrt_dist)[3:ncol(sqrt_dist)])
  colnames(cube_root_dist)[3:ncol(cube_root_dist)]<-paste0("cube_root_trans_",colnames(cube_root_dist)[3:ncol(cube_root_dist)])
  colnames(turkey_dist)[3:ncol(turkey_dist)]<-paste0("turkey_trans_",colnames(turkey_dist)[3:ncol(turkey_dist)])
  colnames(boxcox_dist)[3:ncol(boxcox_dist)]<-paste0("boxcox_trans_",colnames(boxcox_dist)[3:ncol(boxcox_dist)])
  colnames(yeo_johnson_dist)[3:ncol(yeo_johnson_dist)]<-paste0("yeo_johnson_trans_",colnames(yeo_johnson_dist)[3:ncol(yeo_johnson_dist)])
  colnames(order_norm_dist)[3:ncol(order_norm_dist)]<-paste0("order_norm_trans_",colnames(order_norm_dist)[3:ncol(order_norm_dist)])
  colnames(lambertw_dist)[3:ncol(lambertw_dist)]<-paste0("lambertw_trans_",colnames(lambertw_dist)[3:ncol(lambertw_dist)])

  positive_skewed_transformation<-join_all(list(log_dist,sqrt_dist,cube_root_dist,turkey_dist,boxcox_dist,yeo_johnson_dist,order_norm_dist,lambertw_dist),type="left",by=c("names","distribution" ))

  if (ChooseBestTrans == "Pearson P Value")
  {
    # print("Pearson P value chosen as metric for best transformaion fit")
    method_chosen<- data.table(positive_skewed_transformation)[,lapply(.SD,function(x){(x)}),.SDcols=seq(7,length(positive_skewed_transformation),5)]
    colnames(method_chosen)<-c("log_transformed","sqrt_transformed","cube_root_transformed","turkey_transformed", "boxcox_transformed", "yeo_johnson_transformed", "order_norm_transformed", "lambertw_transformed")
    positive_skewed_transformation$chosen_method<-colnames(method_chosen)[apply(method_chosen,1,which.min)]
  } else if (ChooseBestTrans == "Min skewness")
  {
    # print("Minimum Skewness value chosen as metric for best transformaion fit")
    method_chosen<- data.table(positive_skewed_transformation)[,lapply(.SD,function(x){(x)}),.SDcols=seq(5,length(positive_skewed_transformation),5)]
    colnames(method_chosen)<-c("log_transformed","sqrt_transformed","cube_root_transformed","turkey_transformed", "boxcox_transformed", "yeo_johnson_transformed", "order_norm_transformed", "lambertw_transformed")
    positive_skewed_transformation$chosen_method<-colnames(method_chosen)[apply(method_chosen,1,which.min)]
  } else
  {
    print("'ChooseBestTrans' is not provided. Defaulted to 'Shapiro P value'")
    method_chosen<- data.table(positive_skewed_transformation)[,lapply(.SD,function(x){(x)}),.SDcols=seq(5,length(positive_skewed_transformation),5)]
    colnames(method_chosen)<-c("log_transformed","sqrt_transformed","cube_root_transformed","turkey_transformed", "boxcox_transformed", "yeo_johnson_transformed", "order_norm_transformed", "lambertw_transformed")
    positive_skewed_transformation$chosen_method<-colnames(method_chosen)[apply(method_chosen,1,which.max)]
  }

  positive_skewed_transformation<-left_join(positive_skewed_transformation,turkey_lamda,by="names")
  positive_skewed_transformation<-left_join(positive_skewed_transformation,boxcox_lamda,by="names")
  positive_skewed_transformation<-left_join(positive_skewed_transformation,yeo_johnson_lamda,by="names")
  positive_skewed_transformation<-left_join(positive_skewed_transformation,minvalues,by="names")
  positive_skewed_transformation$var<-ifelse(positive_skewed_transformation$chosen_method=="turkey_transformed",1,NA)
  positive_skewed_transformation$turkey_lamda<-positive_skewed_transformation$turkey_lamda*positive_skewed_transformation$var
  positive_skewed_transformation$var<-NULL
  positive_skewed_transformation$var_bc<-ifelse(positive_skewed_transformation$chosen_method=="boxcox_transformed",1,NA)
  positive_skewed_transformation$boxcox_lamda<-positive_skewed_transformation$boxcox_lamda*positive_skewed_transformation$var_bc
  positive_skewed_transformation$var_bc<-NULL
  positive_skewed_transformation$var_yj<-ifelse(positive_skewed_transformation$chosen_method=="yeo_johnson_transformed",1,NA)
  positive_skewed_transformation$yeo_johnson_lamda<-positive_skewed_transformation$yeo_johnson_lamda*positive_skewed_transformation$var_yj
  positive_skewed_transformation$var_yj<-NULL


  log_list<-data.table(positive_skewed_transformation)[chosen_method=="log_transformed",]$names
  sqrt_list<-data.table(positive_skewed_transformation)[chosen_method=="sqrt_transformed",]$names
  cube_root_list<-data.table(positive_skewed_transformation)[chosen_method=="cube_root_transformed",]$names
  turkey_list<-data.table(positive_skewed_transformation)[chosen_method=="turkey_transformed",]$names
  boxcox_list<-data.table(positive_skewed_transformation)[chosen_method=="boxcox_transformed",]$names
  yeo_johnson_list<-data.table(positive_skewed_transformation)[chosen_method=="yeo_johnson_transformed",]$names
  order_norm_list<-data.table(positive_skewed_transformation)[chosen_method=="order_norm_transformed",]$names
  lambertw_list<-data.table(positive_skewed_transformation)[chosen_method=="lambertw_transformed",]$names

  positive_df<-cbind(data.frame(log_transformed)[,log_list],
                     data.frame(sqrt_transformed)[,sqrt_list],
                     data.frame(cube_root_transformed)[,cube_root_list],
                     data.frame(turkey_transformed)[,turkey_list],
                     data.frame(boxcox_transformed)[,boxcox_list],
                     data.frame(yeo_johnson_transformed)[,yeo_johnson_list],
                     data.frame(order_norm_transformed)[,order_norm_list],
                     data.frame(lambertw_transformed)[,lambertw_list])
  colnames(positive_df)<-c(log_list,sqrt_list,cube_root_list,turkey_list, boxcox_list, yeo_johnson_list, order_norm_list, lambertw_list)

  rm(cube_root_dist,log_dist,turkey_dist,boxcox_dist,method_chosen,sqrt_dist,turkey_transformed,turkey_lamda,boxcox_transformed,boxcox_lamda,sqrt_transformed,cube_root_transformed,log_transformed,boxcox_list,cube_root_list,log_list,sqrt_list,turkey_list, yeo_johnson_dist, yeo_johnson_list, yeo_johnson_transformed, yeo_johnson_lamda, order_norm_dist, order_norm_list, order_norm_transformed, lambertw_dist, lambertw_list, lambertw_transformed)


  ####Note to self: Need to return positive_skewed_transformation & positive_df


  #Addressing Negetive/Left Skew (Square,Cube,turkeypower)
  sqr_transformed<-data.table(data)[,lapply(.SD,function(x){(x^2)}),.SDcols=negetive_skewed]
  cube_transformed<-data.table(data)[,lapply(.SD,function(x){sign(x)*abs(x)^(3)}),.SDcols=negetive_skewed]
  turkey_transformed<-data.table(data)[,lapply(.SD,function(x){transformTukey(x,plotit=FALSE,quiet = T)}),.SDcols=negetive_skewed]
  turkey_lamda<-data.table(data)[,lapply(.SD,function(x){transformTukey(x,plotit=FALSE,returnLambda = T,quiet = T)}),.SDcols=negetive_skewed]
  turkey_lamda<-t(turkey_lamda)
  turkey_lamda<-data.frame(turkey_lamda)
  turkey_lamda$names<-rownames(turkey_lamda)

  boxcox_transformed<-data.table(data)[,lapply(.SD,function(x){boxcoxtransForecast(x)}),.SDcols=negetive_skewed]
  boxcox_lamda<-data.table(data)[,lapply(.SD,function(x){boxcoxtransForecast(x,returnLambda = TRUE)}),.SDcols=negetive_skewed]
  boxcox_lamda<-t(boxcox_lamda)
  boxcox_lamda<-data.frame(boxcox_lamda)
  boxcox_lamda$names<-rownames(boxcox_lamda)

  yeo_johnson_transformed<-data.table(data)[,lapply(.SD,function(x){yeoJohnsonVGAM(x)}),.SDcols=negetive_skewed]
  yeo_johnson_lamda<-data.table(data)[,lapply(.SD,function(x){yeoJohnsonVGAM(x,returnLambda = TRUE)}),.SDcols=negetive_skewed]
  yeo_johnson_lamda<-t(yeo_johnson_lamda)
  yeo_johnson_lamda<-data.frame(yeo_johnson_lamda)
  yeo_johnson_lamda$names<-rownames(yeo_johnson_lamda)

  order_norm_transformed <- data.table(data)[,lapply(.SD,function(x){bestNormalize::orderNorm(x, warn = FALSE)$x.t}),.SDcols=negetive_skewed]
  for (i in 1:length(negetive_skewed))
  {
    order_norm_objects[[negetive_skewed[i]]] <- bestNormalize::orderNorm(data.frame(data)[,negetive_skewed[i]], warn = FALSE)
  }

  lambertw_transformed<-data.table(data)[,lapply(.SD,function(x){GaussianizeLambertW(x)}),.SDcols=negetive_skewed]
  colnames(lambertw_transformed) <- negetive_skewed
  for (i in 1:length(negetive_skewed))
  {
    lambertw_tau_mat[[negetive_skewed[i]]] <- GaussianizeLambertW(data.frame(data)[,negetive_skewed[i]],return_tau_mat = TRUE)
  }

  sqr_dist<-data.frame(data_distribution(sqr_transformed,dv="NULL"))
  cube_dist<-data.frame(data_distribution(cube_transformed,dv="NULL"))
  turkey_dist<-data.frame(data_distribution(turkey_transformed,dv="NULL"))
  boxcox_dist<-data.frame(data_distribution(boxcox_transformed,dv="NULL"))
  yeo_johnson_dist<-data.frame(data_distribution(yeo_johnson_transformed,dv="NULL"))
  order_norm_dist<-data.frame(data_distribution(order_norm_transformed,dv="NULL"))
  lambertw_dist<-data.frame(data_distribution(lambertw_transformed,dv="NULL"))

  sqr_dist$is_dv<-NULL
  cube_dist$is_dv<-NULL
  turkey_dist$is_dv<-NULL
  boxcox_dist$is_dv<-NULL
  yeo_johnson_dist$is_dv<-NULL
  order_norm_dist$is_dv<-NULL
  lambertw_dist$is_dv<-NULL

  colnames(sqr_dist)[3:ncol(sqr_dist)]<-paste0("sqr_trans_",colnames(sqr_dist)[3:ncol(sqr_dist)])
  colnames(cube_dist)[3:ncol(cube_dist)]<-paste0("cube_trans_",colnames(cube_dist)[3:ncol(cube_dist)])
  colnames(turkey_dist)[3:ncol(turkey_dist)]<-paste0("turkey_trans_",colnames(turkey_dist)[3:ncol(turkey_dist)])
  colnames(boxcox_dist)[3:ncol(boxcox_dist)]<-paste0("boxcox_trans_",colnames(boxcox_dist)[3:ncol(boxcox_dist)])
  colnames(yeo_johnson_dist)[3:ncol(yeo_johnson_dist)]<-paste0("yeo_johnson_trans_",colnames(yeo_johnson_dist)[3:ncol(yeo_johnson_dist)])
  colnames(order_norm_dist)[3:ncol(order_norm_dist)]<-paste0("order_norm_trans_",colnames(order_norm_dist)[3:ncol(order_norm_dist)])
  colnames(lambertw_dist)[3:ncol(lambertw_dist)]<-paste0("lambertw_trans_",colnames(lambertw_dist)[3:ncol(lambertw_dist)])

  negetive_skewed_transformation<-join_all(list(sqr_dist,cube_dist,turkey_dist,boxcox_dist,yeo_johnson_dist, order_norm_dist, lambertw_dist),type="left",by=c("names","distribution" ))

  if (ChooseBestTrans == "Pearson P Value")
  {
    # print("Pearson P value chosen as metric for best transformaion fit")
    method_chosen<- data.table(negetive_skewed_transformation)[,lapply(.SD,function(x){(x)}),.SDcols=seq(7,length(negetive_skewed_transformation),5)]
    colnames(method_chosen)<-c("sqr_transformed","cube_transformed","turkey_transformed","boxcox_transformed", "yeo_johnson_transformed", "order_norm_transformed", "lambertw_transformed")
    negetive_skewed_transformation$chosen_method<-colnames(method_chosen)[apply(method_chosen,1,which.min)]
  } else if (ChooseBestTrans == "Min skewness")
  {
    # print("Minimum Skewness value chosen as metric for best transformaion fit")
    method_chosen<- data.table(negetive_skewed_transformation)[,lapply(.SD,function(x){(x)}),.SDcols=seq(5,length(negetive_skewed_transformation),5)]
    colnames(method_chosen)<-c("sqr_transformed","cube_transformed","turkey_transformed","boxcox_transformed", "yeo_johnson_transformed", "order_norm_transformed", "lambertw_transformed")
    negetive_skewed_transformation$chosen_method<-colnames(method_chosen)[apply(method_chosen,1,which.min)]
  } else
  {
    # print("Shapiro P value chosen as metric for best transformaion fit")
    method_chosen<- data.table(negetive_skewed_transformation)[,lapply(.SD,function(x){(x)}),.SDcols=seq(5,length(negetive_skewed_transformation),5)]
    colnames(method_chosen)<-c("sqr_transformed","cube_transformed","turkey_transformed","boxcox_transformed", "yeo_johnson_transformed", "order_norm_transformed", "lambertw_transformed")
    negetive_skewed_transformation$chosen_method<-colnames(method_chosen)[apply(method_chosen,1,which.max)]
  }

  # method_chosen<- data.table(negetive_skewed_transformation)[,lapply(.SD,function(x){abs(x)}),.SDcols=c(4,9,14,19,24,29,34)]
  # colnames(method_chosen)<-c("sqr_transformed","cube_transformed","turkey_transformed","boxcox_transformed", "yeo_johnson_transformed", "order_norm_transformed", "lambertw_transformed")
  # negetive_skewed_transformation$chosen_method<-colnames(method_chosen)[apply(method_chosen,1,which.max)]
  negetive_skewed_transformation<-left_join(negetive_skewed_transformation,turkey_lamda,by="names")
  negetive_skewed_transformation<-left_join(negetive_skewed_transformation,boxcox_lamda,by="names")
  negetive_skewed_transformation<-left_join(negetive_skewed_transformation,yeo_johnson_lamda,by="names")
  negetive_skewed_transformation<-left_join(negetive_skewed_transformation,minvalues,by="names")
  negetive_skewed_transformation$var<-ifelse(negetive_skewed_transformation$chosen_method=="turkey_transformed",1,NA)
  negetive_skewed_transformation$turkey_lamda<-negetive_skewed_transformation$turkey_lamda*negetive_skewed_transformation$var
  negetive_skewed_transformation$var<-NULL
  negetive_skewed_transformation$var_bc<-ifelse(negetive_skewed_transformation$chosen_method=="boxcox_transformed",1,NA)
  negetive_skewed_transformation$boxcox_lamda<-negetive_skewed_transformation$boxcox_lamda*negetive_skewed_transformation$var_bc
  negetive_skewed_transformation$var_bc<-NULL
  negetive_skewed_transformation$var_yj<-ifelse(negetive_skewed_transformation$chosen_method=="yeo_johnson_transformed",1,NA)
  negetive_skewed_transformation$yeo_johnson_lamda<-negetive_skewed_transformation$yeo_johnson_lamda*negetive_skewed_transformation$var_yj
  negetive_skewed_transformation$var_yj<-NULL

  sqr_list<-data.table(negetive_skewed_transformation)[chosen_method=="sqr_transformed",]$names
  cube_list<-data.table(negetive_skewed_transformation)[chosen_method=="cube_transformed",]$names
  turkey_list<-data.table(negetive_skewed_transformation)[chosen_method=="turkey_transformed",]$names
  boxcox_list<-data.table(negetive_skewed_transformation)[chosen_method=="boxcox_transformed",]$names
  yeo_johnson_list<-data.table(negetive_skewed_transformation)[chosen_method=="yeo_johnson_transformed",]$names
  order_norm_list<-data.table(negetive_skewed_transformation)[chosen_method=="order_norm_transformed",]$names
  lambertw_list<-data.table(negetive_skewed_transformation)[chosen_method=="lambertw_transformed",]$names

  negetive_df<-cbind(data.frame(sqr_transformed)[,sqr_list],
                     data.frame(cube_transformed)[,cube_list],
                     data.frame(turkey_transformed)[,turkey_list],
                     data.frame(boxcox_transformed)[,boxcox_list],
                     data.frame(yeo_johnson_transformed)[,yeo_johnson_list],
                     data.frame(order_norm_transformed)[,order_norm_list],
                     data.frame(lambertw_transformed)[,lambertw_list])
  colnames(negetive_df) <- c(sqr_list, cube_list, turkey_list, boxcox_list, yeo_johnson_list, order_norm_list, lambertw_list)

  rm(cube_dist,turkey_dist,boxcox_dist,method_chosen,sqr_dist,turkey_transformed,boxcox_transformed,sqr_transformed,cube_transformed,minvalues,boxcox_list,cube_list,sqr_list,turkey_list, yeo_johnson_dist, yeo_johnson_lamda, yeo_johnson_list, yeo_johnson_transformed, order_norm_dist, order_norm_transformed, order_norm_list, lambertw_dist, lambertw_list, lambertw_transformed)


  ####Note to self: Need to return negetive_skewed_transformation & negetive_df


  transformed_df<-data.frame(cbind(positive_df,normal_df,negetive_df))
  complete_data<-data.frame(cbind(transformed_df,data.frame(data_org)[,(!colnames(data_org) %in% colnames(transformed_df))]))

  #Normalize Data
  Scaled_df<-data.table(transformed_df)[,lapply(.SD,function(x){((x-min(x))/(max(x)-min(x)))}),.SDcols=c(normal_var,negetive_skewed,positive_skewed)]

  # generating a model fit file
  pos_skew_tras <- positive_skewed_transformation[,c(1,(length(positive_skewed_transformation)-4):(length(positive_skewed_transformation)))]
  neg_skew_tras <- negetive_skewed_transformation[,c(1,(length(negetive_skewed_transformation)-4):(length(negetive_skewed_transformation)))]
  trans_fit <- rbind(pos_skew_tras,neg_skew_tras)

  ## Return ordernorm objects list
  return(list(transformed_df = transformed_df,
              scaled_df = Scaled_df,
              original_dist = distribution_orginal,
              neg_skew_trans = negetive_skewed_transformation,
              pos_skew_trans = positive_skewed_transformation,
              trans_fit_model = list(trans_fit = trans_fit,
                                     order_norm_obj = order_norm_objects,
                                     lambertw_tau_mat = lambertw_tau_mat),
              complete_data = complete_data))

}
