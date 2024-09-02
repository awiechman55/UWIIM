##############Calculate Sensitivity Metric 
############Background Function
####Calculate Nominal Performance Metrics####
calc_reliab_avg_nom <- function(c)
{
  rez_mag <- load_rez(c,1)
  return(tail(rez_mag$reliab_avg,1))
}

calc_reliab_total_nom <- function(c)
{
  rez_mag <- load_rez(c,1)
  return(tail(rez_mag$reliab_total,1))
}

calc_reliab_min_nom <- function(c)
{
  rez_mag <- load_rez(c,1)
  return(tail(rez_mag$reliab_min,1))
}

calc_fin_nom <- function(c)
{
  rez_mag <- load_rez(c,1)
  rez_mag <- fixColNames(rez_mag, 1)
  return(tail(rez_mag$Ending_PC_Revenue,1))
}

calc_d_nom <- function(c)
{
  rez_mag <- load_rez(c,1)
  rez_mag <- fixColNames(rez_mag, 1)
  return(tail(rez_mag$Ending_PC_Demand,1))
}

####Calculate Sensitivity Table for City, Variable 1, Variable 2#####
calc_SensMetric <- function(c,v1,v2=NULL)
{
  ##Load Raw Results Table
  df <- fixColNames(load_rez(c,v1,v2),ifelse(is.null(v2),1,2))
  
  ##Generate Nominal Performance Metrics
  reliab_avg_nom <- calc_reliab_avg_nom(c) 
  reliab_total_nom <- calc_reliab_total_nom(c)
  reliab_min_nom <- calc_reliab_min_nom(c)
  fin_nom <- calc_fin_nom(c)
  d_nom <- calc_d_nom(c)
  
  ##Calculate Sensitivity 
  if(v1==1)
  {
    df_sens <- df %>% mutate(Rel_Avg = ifelse(MagCAPShock==0,0,((Reliab_Avg-reliab_avg_nom)/reliab_avg_nom)/MagCAPShock),
                             #Rel_Total = ifelse(MagCAPShock==0,0,((Reliab_Total-reliab_total_nom)/reliab_total_nom)/MagCAPShock),
                             Rel_Min = ifelse(MagCAPShock==0,0,((Reliab_Min-reliab_min_nom)/reliab_min_nom)/MagCAPShock),
                             d_End = ifelse(MagCAPShock==0,0,((Ending_PC_Demand-d_nom)/d_nom)/MagCAPShock),
                             Rates_End = ifelse(MagCAPShock==0,0,(-(Ending_PC_Revenue-fin_nom)/fin_nom)/MagCAPShock))
  }
  else
  {
    def_CAP <- -0.284
    df_sens <- df %>% mutate(Rel_Avg = ((Reliab_Avg-reliab_avg_nom)/reliab_avg_nom)/def_CAP,
                             #Rel_Total = ((Reliab_Total-reliab_total_nom)/reliab_total_nom)/def_CAP,
                             Rel_Min = ((Reliab_Min-reliab_min_nom)/reliab_min_nom)/def_CAP,
                             d_End = ((Ending_PC_Demand-d_nom)/d_nom)/def_CAP,
                             Rates_End = (-(Ending_PC_Revenue-fin_nom)/fin_nom)/def_CAP)
  }
  
  ##Add City and Parameter Names Columns
  df_sens$City <- cities$City[c]
  df_sens$Param1_Name <- colnames(df_sens)[1]
  colnames(df_sens)[1] <- "Param1_Value"
  if(!is.null(v2))
  {
    df_sens$Param2_Name <- colnames(df_sens)[2]
    colnames(df_sens)[2] <- "Param2_Value"
    
    ##Extract Needed Columns for Sensitivity Table
    df_sens_final <- df_sens[,c("City","Param1_Name","Param1_Value","Param2_Name",
                                "Param2_Value", "Rel_Avg", "Rel_Min", 
                                "Rates_End","d_End")] 
  }else{
    ##Extract Needed Columns for Sensitivity Table
    df_sens_final <- df_sens[,c("City","Param1_Name","Param1_Value",
                                "Rel_Avg", "Rel_Min", "Rates_End","d_End")] 
  }
    
  return(df_sens_final)
}

####Calculate Sensitivity Table for Multiple Cities####
sens_Table <- function(cs,v1s,v2s=NULL)
{
  df_sens <- data.frame(City=character(),Param1_Name=character(),
                        Param1_Value=numeric(),Param2_Name=character(),
                        Param2_Value=numeric(),Rel_Avg=numeric(),
                        Rel_Min=numeric(),Rates_End=numeric(),d_End=numeric())
  
  for(c in 1:length(cs))
  {
    for(v1 in 1:length(v1s))
    {
      if(is.null(v2s))
      {
        df_sens_entry <- calc_SensMetric(cs[c],v1s[v1])
        df_sens <- rbind(df_sens,df_sens_entry)
      }else{
        for(v2 in 1:length(v2s))
        {
          df_sens_entry <- calc_SensMetric(cs[c],v1s[v1],v2s[v2])
          df_sens <- rbind(df_sens,df_sens_entry)
        }
      }
    }
  }
  
  return(df_sens)
}

####Calculate Default Average Sensitivities####
default_sens_avg <- function(cs, badShocks=FALSE, badCutoff_PS=-0.5, badCutoff_Q=-0.2)
{
  if(badShocks){
    sens_Table_list <- list()
    counter <- 1
    
    if(1 %in% cs){
      sens_Table_list[[counter]] <- filter(sens_Table(1,1), Param1_Value<=badCutoff_PS)
      counter <- counter + 1
    }
    if(2 %in% cs){
      sens_Table_list[[counter]] <- filter(sens_Table(2,1), Param1_Value<=badCutoff_PS)
      counter <- counter + 1
    }
    if(3 %in% cs){
      sens_Table_list[[counter]] <- filter(sens_Table(3,1), Param1_Value>=badCutoff_Q)
    }
    
    sens_Table_filt <- bind_rows(sens_Table_list)
    sens_avg <- sens_Table_filt %>% gather(Metric, Sens, 4:7) %>%
      group_by(City,Metric) %>% summarise(Default=mean(Sens))
  }else{
    sens_avg <- sens_Table(cs,1) %>% gather(Metric, Sens, 4:7) %>%
      group_by(City,Metric) %>% summarise(Default=mean(Sens))
  }
  
  return(sens_avg)
}


####Calculate Deviation of Average Sensitivity Table Given Sensitivity Table####
sens_avg_df <- function(cs, vs, badShocks=TRUE, badCutoff_PS=-0.5, badCutoff_Q=-0.2)
{
  ##Generate Default Table for Comparison
  default <- default_sens_avg(cs,badShocks=badShocks,badCutoff_PS=badCutoff_PS,
                              badCutoff_Q=badCutoff_Q)
  
  ##Collect Sensitivity Tables for all desired parameters
  df_sens <- sens_Table(cs,1,vs)
 
  ##Calculate Average Sensitivity Over All CAP Shortage Magnitudes
  if(badShocks)
  {
    
    sens_Table_list <- list()
    counter <- 1
    if(1 %in% cs){
      sens_Table_list[[counter]] <- filter(df_sens, City=="PHX" &
                                             Param1_Value<=badCutoff_PS)
      counter <- counter + 1
    }
    if(2 %in% cs){
      sens_Table_list[[counter]] <- filter(df_sens, City=="Sc" &
                                             Param1_Value<=badCutoff_PS)
      counter <- counter + 1
    }
    if(3 %in% cs){
      sens_Table_list[[counter]] <- filter(df_sens, City=="QC" &
                                             Param1_Value>=badCutoff_Q)
    }
      
    sens_Table_filt <- bind_rows(sens_Table_list)
    
    df_sens_avg <- sens_Table_filt %>% gather(Metric, Sens, 6:9) %>% 
      group_by(City,Metric,Param2_Name,Param2_Value) %>% 
      summarise(Avg_Sens=mean(Sens))
  }else{
    df_sens_avg <- df_sens %>% gather(Metric, Sens, 6:9) %>% 
      group_by(City,Metric,Param2_Name,Param2_Value) %>% 
      summarise(Avg_Sens=mean(Sens))
  }
  
  ##Calculate Deviation from Default
  df_sens_avg <- left_join(df_sens_avg, default, by=c("City","Metric"))
  df_sens_avg <- df_sens_avg %>% mutate(Dev=Avg_Sens-Default)
  
  return(df_sens_avg)
}





