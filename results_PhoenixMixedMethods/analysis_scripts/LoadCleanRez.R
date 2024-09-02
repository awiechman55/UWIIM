##################Background Functions
####Load Results Data
#v and c are indices for the variable and city vectors
load_rez <- function(c, v1, v2=NULL, file_preamble="./Raw/", three_d=FALSE, shock_type = NULL)
{
  if(is.null(v2))
  {
    file_name <- paste(file_preamble,p_names$p[v1],"_",cities$City[c],".csv",sep="")
  }
  else
  {
    if(three_d)
    {
      file_name <- paste(file_preamble,p_names$p[v1],"_",p_names$p[v2],"_",cities$City[c],"_3d",sep="")
      if(shock_type == 2){
        file_name <- paste(file_name,"_highmag.csv",sep="")
      }else if(shock_type == 3){
        file_name <- paste(file_name,"_grad.csv",sep="")
      }else{
        file_name <- paste(file_name,".csv",sep="")
      }
      
    }else{
      file_name <- paste(file_preamble,p_names$p[v1],"_",p_names$p[v2],"_",cities$City[c],".csv",sep="")
    }
  }
  
  rez <- read.csv(file_name)
  
  return(rez)
}

####Fix Column Names 

fixColNames <- function(df,num_d,shock_type=NULL)
{
  if(num_d==1){
    colnames(df)[2:ncol(df)] <- metrics$Metric
  }else if(num_d==2){
    colnames(df)[3:ncol(df)] <- metrics$Metric
  }else if(shock_type ==3){
    colnames(df)[3:(ncol(df)-1)] <- metrics$Metric
    colnames(df)[ncol(df)] <- p_names$p[2]
  }else{
    colnames(df)[3:ncol(df)] <- metrics$Metric
  }
  
  return(df)
}
