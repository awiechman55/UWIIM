####Sensitivity Heat Map (Raw)####
plot_Sens_heat <- function(cs,v1,v2,metric=NA,noavg=TRUE,v1_def=NULL,v2_def=NULL,
                           cmin=NA,cmax=NA)
{
  df_heat <- sens_Table(cs,v1,v2)
  
  plt <- plot_heat(df_heat, cs, v1, v2, metric=metric,noavg=noavg,v1_def=v1_def,
                   v2_def=v2_def, cmin=cmin,cmax=cmax)
}

###Plot Heatmap given table
plot_heat <- function(df_heat, cs, v1, v2, metric=metric,noavg=noavg,v1_def=v1_def,
                      v2_def=v2_def, cmin=cmin,cmax=cmax,IF_thresh=FALSE,goal=NULL)
{
  df_heat_g <- df_heat %>% gather(Metric,Sens, 6:10)
  
  if(metric=="Shortage"){
    df_heat_g <- filter(df_heat_g,Metric=="Rel_Avg")
    df_heat_g <- df_heat_g %>% mutate(Sens = Sens*50)
  }else{
    df_heat_g <- filter(df_heat_g,Metric==metric)
  }
  
  #return(df_heat_g)
  
  #df_heat_g$Sens[which(df_heat_g$Sens>cmax)] <- NA
  #df_heat_g$Sens[which(df_heat_g$Sens<cmin)] <- NA
  
  if(IF_thresh){
    df_heat_g <- df_heat_g %>% mutate(Param2_Value= goal - Param2_Value)
  }
  
  if(p_names$p[v1] == df_heat$Param1_Name[1]){
    if(p_names$p[v2]==df_heat$Param2_Name[1]){
      plt <- filter(df_heat_g,!is.na(Sens)) %>% 
        ggplot(aes(y=Param1_Value,x=Param2_Value,fill=Sens)) + 
        xlab(p_labels$p[v2]) + ylab(p_labels$p[v1]) 
    }else{
      plt <- filter(df_heat_g,!is.na(Sens)) %>% 
        ggplot(aes(x=Param2_Value,y=ShockParam_Value,fill=Sens)) + 
        xlab(p_labels$p[v2]) + ylab(p_labels$p[v1]) 
    }
  }else{
    if(p_names$p[v2]==df_heat$Param1_Name[1]){
      plt <- filter(df_heat_g,!is.na(Sens)) %>% 
        ggplot(aes(x=Param1_Value,y=ShockParam_Value,fill=Sens)) + 
        xlab(p_labels$p[v2]) + ylab(p_labels$p[v1])
    }else{
      plt <- filter(df_heat_g,!is.na(Sens)) %>% 
        ggplot(aes(x=Param2_Value,y=ShockParam_Value,fill=Sens)) + 
        xlab(p_labels$p[v2]) + ylab(p_labels$p[v1])
    }
  }
  
  #return(df_heat_g)
  
  plt <- plt + geom_tile() + 
    theme_classic(base_size=18,base_family = "serif") 
  
  if(is.na(cmin)){
    plt <- plt + scale_fill_distiller(type="seq",direction=-1,palette="RdBu",
                                      limits=c(min(df_heat_g$Sens),max(df_heat_g$Sens)))
  }else{
    if(metric=="Rates_End"){
      plt <- plt + scale_fill_distiller(type="seq",direction=-1,palette="RdBu",
                                        limits=c(cmin,cmax),name="") +
        new_scale_fill() +
        geom_tile(data=filter(df_heat_g,Sens<cmin),aes(fill="#2166ac")) +
        scale_fill_manual(name="",labels=paste("<",cmin),values="#2166ac") +
        new_scale_fill() +
        geom_tile(data=filter(df_heat_g,Sens>cmax),aes(fill="#b2182b")) +
        scale_fill_manual(name="",labels=paste(">", cmax),values="#b2182b")
    }else{
      plt <- plt + scale_fill_distiller(type="seq",direction=-1,palette="RdBu",
                                        limits=c(cmin,cmax),name="") +
        new_scale_fill() +
        geom_tile(data=filter(df_heat_g,Sens>cmax),aes(fill="#b2182b")) +
        scale_fill_manual(name="",labels=paste(">",cmax),values="#b2182b") 
    }
  }
  
  #return(df_heat_g)
  
  if(v1 %in% c(7,8,9))
  {
    plt <- plt + scale_y_continuous(trans='log', breaks=c(4,22,110,220)) 
  }
  else if(v2 %in% c(7,8,9))
  {
    plt <- plt + scale_x_continuous(trans='log', breaks=c(4,22,110,220))
  }

  if(v1 == 1)
  {
    plt <- plt + scale_y_reverse()
  }
  if(IF_thresh){
    plt <- plt + scale_x_reverse()
  }
  if(!is.null(v1_def)){
    plt <- plt+geom_vline(xintercept=v1_def,lty="dashed")
  }
  if(!is.null(v2_def)){
    plt <- plt+geom_hline(yintercept=v2_def,lty="dashed")
  }
  
  return(plt)
}

##### Heatmap over average of high magnitude (given 3d data) #####
plot_Sens_heat_highmag <- function(cs, v1, v2, metric="Rel_Min",noavg=TRUE,v1_def=NULL,
                                   v2_def=NULL, cmin=NA, cmax=NA, badCutoff_PS=0)
{
  df_heat <- sens_avg_df_3d(cs,v1,v2,shock_type=1,avg_v=1,set_v_value=NULL,
                            badCutoff_PS=badCutoff_PS)
  
  df_heat_g <- filter(df_heat,Metric==metric)
  #df_heat_g$Sens[which(df_heat_g$Sens>1)] <- NA
  #df_heat_g$Sens[which(df_heat_g$Sens< -1)] <- NA
  
  if(v1 %in% c(7,8,9)){
    df_heat_g <- filter(df_heat_g, Param1_Value != 22)
  }
  
  
  plt <- df_heat_g %>% ggplot(aes(y=Param1_Value,x=Param2_Value,fill=Dev)) + 
    geom_tile() + 
    theme_classic(base_size=18,base_family = "serif") + 
    xlab(p_labels$p[v2]) + ylab(p_labels$p[v1]) 
  
  if(is.na(cmin)){
    plt <- plt + scale_fill_distiller(type="seq",direction=-1,palette="RdBu",
                                      limits=c(min(df_heat_g$Dev),max(df_heat_g$Dev)))
  }else{
    if(metric=="Rates_End"){
      plt <- plt + scale_fill_distiller(type="seq",direction=-1,palette="RdBu",
                                        limits=c(cmin,cmax),name="Avg. Dev.") +
        new_scale_fill() +
        geom_tile(data=filter(df_heat_g,Dev<cmin),aes(fill="#2166ac")) +
        scale_fill_manual(name="Avg. Dev",labels=paste("<",cmin),values="#2166ac")
    }else{
      plt <- plt + scale_fill_distiller(type="seq",direction=-1,palette="RdBu",
                                        limits=c(cmin,cmax),name="Avg. Dev.") +
        new_scale_fill() +
        geom_tile(data=filter(df_heat_g,is.na(Dev)),aes(fill="#b2182b")) +
        scale_fill_manual(name="Avg. Dev",labels=">1",values="#b2182b") 
    }
    
    
  }
  
  if(v1 %in% c(7,8,9))
  {
    plt <- plt + scale_y_continuous(trans='log', breaks=c(4,22,110,220)) 
  }
  else if(v2 %in% c(7,8,9))
  {
    plt <- plt + scale_x_continuous(trans='log', breaks=c(4,22,110,220))
  }
  
  if(v1 == 1)
  {
    plt <- plt + scale_y_reverse()
  }
  if(!is.null(v1_def)){
    plt <- plt+geom_vline(xintercept=v1_def,lty="dashed")
  }
  if(!is.null(v2_def)){
    plt <- plt+geom_hline(yintercept=v2_def,lty="dashed")
  }
  
  return(plt)
}

#####Plot Heatmap with 1 variable at constant (given 3d data)#####
plot_Sens_heat_1constant <- function(cs,file_v1, file_v2, shock_type, v_constant, 
                                     v_constant_value, metric=NA, noavg=TRUE,
                                     v1_def=NULL,v2_def=NULL, cmin=NA,cmax=NA,
                                     IF_thresh=FALSE,goal=NULL)
{
  df_heat_full <- sens_Table(cs, file_v1, file_v2, three_d = TRUE, 
                             shock_type = shock_type)
  
  #return(df_heat_full)
  
  if(file_v1 %in% c(8,9) && v_constant_value !=22){
    df_heat_filt <- filter(df_heat_full, Param1_Value != 22)
  }else if(file_v1 == 7 && v_constant_value !=110){
    df_heat_filt <- filter(df_heat_full, Param1_Value != 110)
  }else{
    if(file_v2 %in% c(8,9) && v_constant_value !=22){
      df_heat_filt <- filter(df_heat_full, Param2_Value != 22)
    }else if(file_v2 == 7 && v_constant_value !=110){
      df_heat_filt <- filter(df_heat_full, Param2_Value != 110)
    }else{
      df_heat_filt <- df_heat_full
    }
  }
  
  if(df_heat_full$Param1_Name[1]==v_constant){
    df_heat_filt <- filter(df_heat_filt, Param1_Value == v_constant_value)
    v2 = file_v2
    v1 = ifelse(shock_type==3,2,1)
  }else if(df_heat_full$Param2_Name[1] == v_constant){
    df_heat_filt <- filter(df_heat_filt, Param2_Value == v_constant_value)
    v2 = file_v1
    v1 = ifelse(shock_type==3,2,1)
  }else{
    df_heat_filt <- filter(df_heat_filt, ShockParam_Value == v_constant_value)
    v1 = file_v1
    v2 = file_v2
  }
  
  #return(df_heat_filt)
  
  
  plt <- plot_heat(df_heat_filt, cs, v1, v2, metric=metric,noavg=noavg,
                   v1_def=v1_def,v2_def=v2_def, cmin=cmin,cmax=cmax,
                   IF_thresh=IF_thresh,goal=goal)
  
  
  
  return(plt)
}

########Plot Sensitivity Deviation####
plot_Dev_heat <- function(cs,v1,v2,shock_type,metric,
                          avg_v=1,badCutoff_PS=0,cmin=NA,cmax=NA,
                          v1_def=NULL,v2_def=NULL,IF_thresh=FALSE,
                          goal=NULL){
  df <- sens_avg_df_3d(cs,v1,v2,shock_type,avg_v=avg_v,badCutoff_PS=badCutoff_PS)
  
  if(IF_thresh){
    df <- df %>% mutate(Param2_Value = goal - Param2_Value)
  }
  
  if(v1 %in% c(8,9)){
    df <- filter(df, Param1_Value != 22)
  }else if(v1 == 7){
    df <- filter(df, Param1_Value != 110)
  }else{
    if(v2 %in% c(8,9)){
      df <- filter(df, Param2_Value != 22)
    }else if(v2 == 7){
      df <- filter(df, Param2_Value != 110)
    }
  }
  
  
  
  if(metric=="Shortage"){
    df_heat <- filter(df,Metric=="Rel_Avg")
    df_heat <- df_heat %>% mutate(Dev=Dev*50)
  }else{
    df_heat <- filter(df,Metric==metric)
  }
  
  if(p_names$p[v1] == df_heat$Param1_Name[1]){
    if(p_names$p[v2]==df_heat$Param2_Name[1]){
      plt <- filter(df_heat,!is.na(Dev)) %>% 
        ggplot(aes(y=Param1_Value,x=Param2_Value,fill=Dev)) + 
        xlab(p_labels$p[v2]) + ylab(p_labels$p[v1]) 
    }else{
      plt <- filter(df_heat,!is.na(Dev)) %>% 
        ggplot(aes(x=Param2_Value,y=ShockParam_Value,fill=Dev)) + 
        xlab(p_labels$p[v2]) + ylab(p_labels$p[v1]) 
    }
  }else{
    if(p_names$p[v2]==df_heat$Param1_Name[1]){
      plt <- filter(df_heat,!is.na(Dev)) %>% 
        ggplot(aes(x=Param1_Value,y=ShockParam_Value,fill=Dev)) + 
        xlab(p_labels$p[v2]) + ylab(p_labels$p[v1])
    }else{
      plt <- filter(df_heat,!is.na(Dev)) %>% 
        ggplot(aes(x=Param2_Value,y=ShockParam_Value,fill=Dev)) + 
        xlab(p_labels$p[v2]) + ylab(p_labels$p[v1])
    }
  }
  
  plt <- plt + geom_tile() + 
    theme_classic(base_size=18,base_family = "serif") 
  
  if(is.na(cmin)){
    plt <- plt + scale_fill_distiller(type="seq",direction=-1,palette="RdBu",
                                      limits=c(min(df_heat$Dev),max(df_heat$Dev)))
  }else{
    if(metric=="Rates_End"){
      plt <- plt + scale_fill_distiller(type="seq",direction=-1,palette="RdBu",
                                        limits=c(cmin,cmax),name="") +
        new_scale_fill() +
        geom_tile(data=filter(df_heat,Dev<cmin),aes(fill="#2166ac")) +
        scale_fill_manual(name="",labels=paste("<",cmin),values="#2166ac") +
        #guides(fill=guide_legend(title=paste("Avg. Dev\n(Over Shocks)\nRate Burden"))) +
        new_scale_fill() +
        geom_tile(data=filter(df_heat,Dev>cmax),aes(fill="#b2182b")) +
        scale_fill_manual(name="",labels=paste(">", cmax),values="#b2182b")
    }else{
      plt <- plt + scale_fill_distiller(type="seq",direction=-1,palette="RdBu",
                                        limits=c(cmin,cmax),name="") +
        new_scale_fill() +
        geom_tile(data=filter(df_heat,Dev>cmax),aes(fill="#b2182b")) +
        scale_fill_manual(name="",labels=paste(">",cmax),values="#b2182b") 
    }
  }
  
  
  if(v1 %in% c(7,8,9))
  {
    plt <- plt + scale_y_continuous(trans='log', breaks=c(4,22,110,220)) 
  }
  else if(v2 %in% c(7,8,9))
  {
    plt <- plt + scale_x_continuous(trans='log', breaks=c(4,22,110,220))
  }
  
  if(v1 == 1)
  {
    plt <- plt + scale_y_reverse()
  }else if(IF_thresh){
    plt <- plt + scale_x_reverse()
  }
  
  if(!is.null(v1_def)){
    plt <- plt+geom_vline(xintercept=v1_def,lty="dashed")
  }
  if(!is.null(v2_def)){
    plt <- plt+geom_hline(yintercept=v2_def,lty="dashed")
  }
  
  return(plt)
}