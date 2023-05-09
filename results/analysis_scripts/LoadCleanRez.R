##################Background Functions
####Load Results Data
#v and c are indices for the variable and city vectors
load_rez <- function(c, v1, v2=NULL, file_preamble="./raw/", three_d=FALSE)
{
  if(is.null(v2))
  {
    file_name <- paste(file_preamble,p_names$p[v1],"_",cities$City[c],".csv",sep="")
  }
  else
  {
    if(three_d)
    {
      file_name <- paste(file_preamble,p_names$p[v1],"_",p_names$p[v2],"_",cities$City[c],"_3d.csv",sep="")
    }else{
      file_name <- paste(file_preamble,p_names$p[v1],"_",p_names$p[v2],"_",cities$City[c],".csv",sep="")
    }
  }
  
  rez <- read.csv(file_name)
  
  return(rez)
}

####Fix Column Names 

fixColNames <- function(df,num_d)
{
  if(num_d==1)
    colnames(df)[2:ncol(df)] <- metrics$Metric
  else
    colnames(df)[3:ncol(df)] <- metrics$Metric
  
  return(df)
}
