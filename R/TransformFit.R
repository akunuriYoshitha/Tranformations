#' Performs a Turkey power tranformation with predefined Lambda
#'
#' Takes in a vector and a lambda value to perform Turkey power transfomation
#' @param x A vector which neds to be tranformed
#' @param lambda A fixed lambda value for turkey power transformation
#' @return Returns log tranformattion if lambda equals zero and a power transormation otherwise
#' @export
#' @usage TurkeyFit(x, lambda)
#'
#' @example
#' a <- c(12,34,234,23,678, 768, 34, 34 ,78)
#' TurkeyFit(a, -5)
#'
TurkeyFit <- function(x, lambda)
{
  if (lambda > 0) {
    TRANS = x^lambda
  }
  if (lambda == 0) {
    TRANS = log(x)
  }
  if (lambda < 0) {
    TRANS = -1 * x^lambda
  }
  return(TRANS)
}

#' Fit function for test data transformations
#'
#' Takes in test data and fits the data based on the model file provided
#' @param data Test data to be tranformed
#' @param trans_fit_model A model file that captures the details of transformations done to train data. This is returned as a list element by VariableTransform function
#' @return A transformed dataset
#' @export
#' @usage TransformFit(data, trans_fit_model)

TransformFit <- function(data, trans_fit_model)
{
  # require(data.table)
  # require(dplyr)
  trans_fit <- trans_fit_model$trans_fit
  tau_mat_fit <- trans_fit_model$lambertw_tau_mat
  order_norm_obj <- trans_fit_model$order_norm_obj

  data_org <- data
  ## padding..
  padding_list <- as.character(data.table(trans_fit)$names)
  data <- data.frame(mapply(function(x,min_value){x+1-min_value}, data[,padding_list],data.table(trans_fit)$min_value))

  ## Transformations
  log_list<-as.character(data.table(trans_fit)[chosen_method=="log_transformed",]$names)
  sqrt_list<-as.character(data.table(trans_fit)[chosen_method=="sqrt_transformed",]$names)
  cu_rt_list<-as.character(data.table(trans_fit)[chosen_method=="cube_root_transformed",]$names)
  turkey_list<-as.character(data.table(trans_fit)[chosen_method=="turkey_transformed",]$names)
  boxcox_list<-as.character(data.table(trans_fit)[chosen_method=="boxcox_transformed",]$names)
  sqr_list<-as.character(data.table(trans_fit)[chosen_method=="sqr_transformed",]$names)
  cube_list<-as.character(data.table(trans_fit)[chosen_method=="cube_transformed",]$names)
  yeo_johnson_list <- as.character(data.table(trans_fit)[chosen_method=="yeo_johnson_transformed",]$names)
  lambertw_list <- as.character(data.table(trans_fit)[chosen_method=="lambertw_transformed",]$names)
  order_norm_list <- as.character(data.table(trans_fit)[chosen_method=="order_norm_transformed",]$names)


  log_trans<-data.table(data)[,lapply(.SD,function(x){log(x)}),.SDcols=log_list]
  sqrt_trans<-data.table(data)[,lapply(.SD,function(x){sqrt(x)}),.SDcols=sqrt_list]
  cube_root_trans<-data.table(data)[,lapply(.SD,function(x){sign(x)*abs(x)^(1/3)}),.SDcols=cu_rt_list]
  sqr_trans<-data.table(data)[,lapply(.SD,function(x){(x^2)}),.SDcols=sqr_list]
  cube_trans<-data.table(data)[,lapply(.SD,function(x){sign(x)*abs(x)^(3)}),.SDcols=cube_list]

  turkey_trans<- data.frame(mapply(TurkeyFit, data[,turkey_list],data.table(trans_fit)[chosen_method=="turkey_transformed",]$turkey_lamda))
  colnames(turkey_trans) <- turkey_list

  boxcox_trans<- data.frame(mapply(forecast::BoxCox, data[,boxcox_list],data.table(trans_fit)[chosen_method=="boxcox_transformed",]$boxcox_lamda))
  colnames(boxcox_trans) <- boxcox_list

  yeo_johnson_trans<- data.frame(mapply(yeo.johnson, data[,yeo_johnson_list],data.table(trans_fit)[chosen_method=="yeo_johnson_transformed",]$yeo_johnson_lamda))
  colnames(yeo_johnson_trans) <- yeo_johnson_list

  # lambertw_trans<- data.frame(mapply(Gaussianize, data[,lambertw_list], tau_mat = tau_mat_fit[lambertw_list]))
  for (i in 1:length(lambertw_list))
  {
    # print(i)
    if(!is.null(tau_mat_fit[lambertw_list[i]][[1]]))
    {
      if (i == 1)
      {
        lambertw_trans <- Gaussianize(data[,lambertw_list[i]],tau.mat = tau_mat_fit[lambertw_list[i]][[1]])
      }else
      {
        lambertw_trans <- cbind(lambertw_trans,Gaussianize(data[,lambertw_list[i]],tau.mat = tau_mat_fit[lambertw_list[i]][[1]]))
      }
    }else
    {
      if (i == 1)
      {
        lambertw_trans <- data[,lambertw_list[i]]
      }else
      {
        lambertw_trans <- cbind(lambertw_trans,data[,lambertw_list[i]])
      }
    }

  }
  lambertw_trans <- data.frame(lambertw_trans)
  colnames(lambertw_trans) <- lambertw_list

  order_norm_trans<- data.frame(mapply(predict, newdata = data[,order_norm_list], object = order_norm_obj[order_norm_list]))
  colnames(order_norm_trans) <- order_norm_list

  transformed_df <- bind_cols(lambertw_trans, order_norm_trans, turkey_trans, boxcox_trans, log_trans, sqrt_trans, cube_root_trans, sqr_trans, cube_trans, yeo_johnson_trans)
  complete_data<-data.frame(cbind(transformed_df,data.frame(data_org)[,(!colnames(data_org) %in% colnames(transformed_df))]))

  return(complete_data)
  # return(transformed_df)
}
