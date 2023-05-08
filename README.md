# Urban Water Infrastructure Investment Model (UWIIM), PMA Version

This repository contains the source code for the model and sensitivity analysis configuration, raw outputs of the sensitivity analysis, R scripts used to analyze the outputs of sensitivity analysis, and the figures used in the published manuscript associated with this version of the Urban Water Infrastructure Investment Model (UWIIM). The model is a general coupled infrastructure dynamical systems model, written in the Julia programming language (version 1.8.4) as a discrete time dynamical system, that incorporates operational infrastructure and political-economic feedback processes governing urban water systems. The model has been parameterized to reflect the water resources, financial, and demand context of three cities in the Phoenix Metropolitan Area (PMA): Phoenix, Scottsdale, and Queen Creek. 

# Running the Model

The UWIIM can be run through the provided Jupyter notebook ("UWIIM_PMA.ipynb") or any Julia compiler with the "UWIIM_PMA.jl" file. Both files contain all background functions with heavily commented explanations that can also be found in the Supporting Information document (see "docs" folder).

After compiling the source code in a Julia compiler or the Jupyter notebook, the one-line function that runs the UWIIM is **run_UWIIM()**, which takes four possible inputs. 

## Model Inputs

1. *setup*: Setup function containing the parameters and initial conditions for the model run. This version offers four setup functions: **Default()**, **Phoenix()**, **Scottsdale()**, **QueenCreek()**. If the setup functions are called without any specified arguments, they will output the default parameter and initial conditions corresponding to that city as laid out in the Supporting Information document. However, a model user may change any of the parameter or initial conditions in the arguments of the setup functions as defined in Tables A1-4. 
2. *t_run*: number of years (time steps) to run the model. The default is 50.
3. *year_0*: initial year of the model run. This is primarily used for the time series plots. The default is 2010.
4. *units*: this specifies the water volume units to use for plotting. The default is acre-feet, "AF," which should be used **for all PMA cities.** Any other unit specification will create dimensional analysis errors in the parameters and initial conditions. However, in the Default setup, one can also use "gal" for this argument to use gallon-based water volume units. 

### Default Setup Changable Parameters and Initial Conditions

***Parameters***

| Full Parameter Name | Model Variable Name | Definition | Units |Default Value | Allowable Range |
| -------------- | ------ | ---------- | ----- |------------- | --------------- |
| Max Surface Streamflow Mean | $`\mu`$_s_max | Max mean surface inflow that the city can seek | Bgal/yr OR AFY | 200000 | $`[0,\infty)`$ |
| Max Ground Streamflow Mean | $`\mu`$_g_max | Max mean ground inflow (recharge) that the city can seek | Bgal/yr OR AFY | 0 | $`[0,\infty)`$ |
| Max Per-Capita Revenue | $`\pi`$_max | Max revenue that city can extract per citizen | Dollars/person/yr | 400 | $`[0, \infty)`$ |
| Streamflow Variation Scenario Type | $`\Delta`$C_v_scen | Variation Change Type (0 = no change, 1 = gradual, 2 = sudden) | Unitless | 0 | $`\{0,1,2\}`$ |
| Surface Streamflow Auto-correlation | $`\rho_s`$ | 1-year lagged auto-correlation in surface streamflow | Unitless | 0.2 | $`[0,1]`$ |
| Ground Inflow Auto-correlation | $`\rho_g`$ | 1-year lagged auto-correlation in ground inflow | Unitless | 0 | $`[0,1]`$ |
| Carrying Capacity | $`\kappa`$ | Population carrying capacity | Unitless | 1500000 | $`(0, \infty)`$ |
| Intrinsic Growth Rate | r | logistically fit, intrinsic growth rate of population | Unitless | 0.1 | $`[0,\infty)`$ |
| Max Delivery Efficiency | $`\eta`$_max | Max attainable delivery efficiency | Unitless | 2.0 | $`[0, 2]`$ |
| Max Surface Storage Capacity | $`\upsilon`$_bar_max_s | Max feasible surface storage capacity (multiple of inflow standard deviation) | Unitless | 12 | $`[0, \infty)`$ |
| Max Ground Storage Capacity | $`\upsilon`$_bar_max_g | Max feasible ground storage capacity (multiple of inflow standard deviation) | Unitless | 40 | $`[0, \infty)`$ |
| Max Legal Use of Ground Storage | a_gv | Max legally allowed use of ground storage (proportion) | Unitless | 1 | $`[0,1]`$ |
| Max Legal Use of Ground Inflow | a_gq | Max legally allowed use of ground inflow (proportion) | Unitless | 1 | $`[0,1]`$ |
| Max Legal Use of Surface Storage | a_sv | Max legally allowed use of surface storage (proportion) | Unitless | 1 | $`[0,1]`$ |
| Max Legal Use of Surface Inflow | a_sq | Max legally allowed use of surface inflow (proportion) | Unitless | 1 | $`[0,1]`$ |
| Min Per Capita Use | d_min | Min possible per-capita demand | (Bgal/yr OR AFY)/person | 0.04 | $`[0, \infty)`$ |
| Background Conservation Rate | $`\delta`$_dbar | Annual decay rate of long-term conservation measures (background) | Unitless | 0.0116 | $`[0,1]`$ |
| Surface Storage Decay Rate | $`\delta`$_v | Annual decay rate of surface storage capacity | Unitless | 0.001 | $`[0,1]`$ |
| Delivery Efficiency Decay Rate | $`\delta`$_$`\eta`$ | Annual decay rate of delivery efficiency | Unitless | 0.001 | $`[0,1]`$ |
| Surface Processing Decay Rate | $`\delta`$_w_s | Annual decay rate of surface processing capacity | Unitless | 0.001 | $`[0,1]`$ |
| Ground Processing Decay Rate | $`\delta`$_w_g | Annual decay rate of ground processing capacity | Unitless | 0.001 | $`[0,1]`$ |
| LT Demand Mgmt Implementation Time | $`\tau`$_d | Time to implement long-term demand management | yrs | 3 | $`[1,\infty)`$ |
| Delivery Efficiency Impl Time | $`\tau`$_$`\eta`$ | Time to implement delievery efficiency improvements | yrs | 4 | $[1,\infty)$ |
| Storage Capacity Impl Time | $`\tau`$_v | Time to implement surface storage capacity improvements | yrs | 5 | $`[1,\infty)`$ |
| Surface Processing Impl Time | $`\tau`$_w_s | Time to implement surface processing capacity improvements | yrs | 3 | $`[1,\infty)`$ |
| Ground Processing Impl Time | $`\tau`$_w_g | Time to implement ground processing capacity improvements | yrs | 3 | $`[1,\infty)`$ |
| Surface Flow Augmentation Impl Time | $`\tau`$_$`\mu`$_s | Time to implement surface flow augmentation improvements | yrs | 4 | $`[1,\infty)`$ |
| Ground Flow Augmentation Impl Time | $`\tau`$_$`\mu`$_g | Time to implement ground flow augmentation improvements | yrs | 4 | $`[1,\infty)`$ |
| Short-Term Goal Supply Sufficiency | $`\gamma`$_1 |  Goal short-term proportion between supply and demand | Unitless | 1 | $`[0,\infty)`$ |
| Long-Term Goal Supply Buffer | $`\gamma`$_2 | Goal long-term proportion between supply and demand | Unitless | 1.2 | $`[0,\infty)`$ |
| Minimum Debt Service Coverage Ratio | $`\gamma`$_3 | Minimum Debt Service Coverage Ratio Allowed in any Year | unitless | 2 | $`[0,\infty)`$ |
| Max Rate Increase | $`\psi`$_r | Max possible proportional increase in rates | Unitless | 0.06 | $`[0,\infty)`$ |
| Delivery Efficiency Priority | $`\beta`$_$`\eta`$ | Proportion of long-term investments directed to delivery efficiency | Unitless | 0.2 | $`[0,1]`$ |
| Storage Capacity Priority | $`\beta`$_v | Proportion of long-term investments directed to surface storage capacity | Unitless | 0.4 | $`[0,1]`$ |
| Surface Processing Priority | $`\beta`$_w_s | Proportion of long-term investments directed to surface processing capacity | Unitless | 0.3 | $`[0,1]`$ |
| Ground Processing Priority | $`\beta`$_w_g | Proportion of long-term investments directed to ground processing capacity | Unitless | 0 | $`[0,1]`$ |
| Surface Augmentation Priority | $`\beta`$_$`\mu`$_s | Proportion of long-term investments directed to surface flow augmentation | Unitless | 0 | $`[0,1]`$ |
| Ground Augmentation Priority | $`\beta`$_$`\mu`$_g | Proportion of long-term investments directed to ground flow augmentation | Unitless | 0 | $`[0,1]`$ |
| Short-Term Curtailment Sensitivity | $`\lambda`$_1 | Sensitivity (inst. friction component) in short-term investment | Unitless | 22 | $`[0, \infty)`$ |
| Long-Term Investment Sensitivity | $`\lambda`$_2 | Sensitivity (inst. friction component) in long-term investment | Unitless | 22 | $`[0, \infty)`$ |
| Rate-Setting Sensitivity | $`\lambda`$_3 | Sensitivity (inst. friction component) in rate-setting | Unitless | 22 | $`[0, \infty)`$ |
| ST Curtailment Activation Threshold | $`\epsilon`$_1 | Threshold for action (inst. friction component) in short-term investment | Unitless | 0 | $`[0, \infty)`$ |
| LT Investment Activation Threshold | $`\epsilon`$_2| Threshold for action (inst. friction component) in long-term investment | Unitless | 0 | $`[0, \infty)`$ |
| Rate-Setting Activation Threshold | $`\epsilon`$_3| Threshold for action (inst. friction component) in rate-setting | Unitless | 0 | $`[0, \infty)`$ |
| Operating Cost Function Coefficient | g_o | Coefficient in operating costs function (see equation) | Dollars/(persons*Vol) | 2.5 |  $`[0, \infty)`$ |
| LT Dem Mgmt Investment Coefficient | g_dbar | Coefficient in LT dem mgmt investment function | AFY/Dollars | 4.4 E-7 | $`[0,\infty)`$ |
| Del Eff Investment Coefficient | g_$`\eta`$ | Coefficient in delivery efficiency investment funciton | AFY/Dollars | 0.0033 | $`[0,\infty)`$ | 
| Storage Capacity Investment Coefficient | g_vbar | Coefficient in storage capacity investment function | AF/Dollars | 0.003 | $`[0, \infty)`$ | 
| SW Proc Capacity Investment Coefficient | g_ws | Coefficient in surface processing capacity investment function | AFY/Dollars | 0.00015 | $`[0,\infty)`$ |
| GW Proc Capacity Investment Coefficient | g_wg | Coefficient in ground processing capacity investment function | AFY/Dollars | 0.00015 | $`[0,\infty)`$ |
| SW Inflow Aug Investment Coefficient | g_$`\mu`$s | Coefficient in surface inflow augmentation investment function | AFY/Dollars | 0.0001 | $`[0,\infty)`$ |
| GW Inflow Aug Investment Coefficient | g_$`\mu`$g | Coefficient in ground inflow augmentation investment function | AFY/Dollars | 0.0001 | $`[0,\infty)`$ |
| Operating Cost Population Scale Factor | z_op | Population scale factor in operating cost function | unitless | 0.563 | $`[0,\infty)`$ |
| Operating Cost Demand Scale Factor | z_od | Demand scale factor in operating cost function | unitless | 0.831 | $`[0,\infty)`$ |
| LT Dem Mgmt Investment Scale Factor | z_dbar | Investment scale factor for long-term demand management | unitless | 1 | $`[0,\infty)`$ |
| Del Eff Investment Scale Factor | z_$`\eta`$ | Investment scale factor for delivery efficiency | unitless | 0.82 | $`[0,\infty)`$ |
| Storage Capacity Investment Scale Factor | z_vbar | Investment scale factor for storage capacity | unitless | 1 | $`[0,\infty)`$ |
| SW Proc Capacity Investment Scale Factor | z_ws | Investment scale factor for processing capacity | unitless | 1 | $`[0,\infty)`$ |
| GW Proc Capacity Investment Scale Factor | z_wg | Investment scale factor for processing capacity | unitless | 1 | $`[0,\infty)`$ |
| SW Inflow Aug Investment Scale Factor | z_$`\mu`$s | Investment scale factor for sw inflow augmentation | unitless | 1 | $`[0,\infty)`$ |
| GW Inflow Aug Investment Scale Factor | z_$`\mu`$g | Investment scale factor for gw inflow augmentation | unitless | 1 | $`[0,\infty)`$ |
| Case Specification | case | Indicate which case type to follow (0 = default, 1 = PHX) | Unitless | 0 |$`\{0,1\}`$ |
| Mean Surface Inflow Change Type | $`\Delta\mu`$_s_type | Type of change to surface mean inflow (0 = none, 1 = gradual, 2 = sudden) | Unitless | 0 | $`\{0,1,2\}`$ |
| Mean Surface Inflow Percent Change | $`\Delta\mu`$_s_pc | Percent change to surface mean inflow | Unitless | 0 | $`[0, \infty)`$ |
| Mean Sufrace Inflow Change Time | $`\Delta\mu`$_s_t | Time that sudden change occurs or time that gradual change will be complete | yrs | 0 | $`[0, \infty)`$ |
| Max Surface Processing | w_max_s | Max attainable surface processing (proportion of storage capacity + mean inflow) | Unitless | 1 | $`[0,1]`$ |
| Max Ground Processing | w_max_g | Max attainable ground processing (proportion of storage capacity + mean inflow) | Unitless | 0.5 | $`[0,1]`$ |
| Base Groundwater Use | $`\theta`$_g | Proportion of Demand that is usually served by groundwater | Unitless | 0.0 | $`[0,1]`$ |
| Projection Years | $`\tau`$_p | Years of Projection for Long-Term Investment Signal | Years | 5 | $`[0, \infty)`$ |
| ST Dem Mgmt Investment Effectiveness Coefficient | $`\alpha`$ | Coefficient for effectiveness in ST conservation investments | unitless | 0.5 | $`[0,\infty)`$ |
| Bond Life | $`\tau`$_b | Life of issued bonds | years | 15 | $`[0,\infty)`$ |
| Bond Interest Rate | i_b | Interest Rate of Issued Bonds | unitless | 0.04 | $`[0,1]`$ |
| Proportion of Rates from Fixed Charges | $`\beta`$_p | Proportion of the expected revenue to come from fixed per user charges | unitless | 0.5 | $`[0,1]`$ |
| Del Eff Investment Dollars Proportion | $`\phi`$__$`\eta`$ | Proportion of investment dollars to delivery efficiency | unitless | 0.6 | $`[0,1]`$ |
| Stor Capac Investment Dollars Proportion | $`\phi`$_v | Proportion of investment dollars to storage capacity | unitless | 0.3 | $`[0,1]`$ |
| Surface Proc Capac Investment Dollars Proportion | $`\phi`$_w_s | Proportion of investment dollars to surface processing capacity | unitless | 0.3 | $`[0,1]`$ |
| Ground Proc Capac Investment Dollars Proportion | $`\phi`$_w_g | Proportion of investment dollars to ground processing capacity | unitless | 0 | $`[0,1]`$ |
| Surface Inflow Investment Dollars Proportion | $`\phi`$_$`\mu`$_s | Proportion of investment dollars to surface inflow | unitless | 0 | $`[0,1]`$ |
| Ground Inflow Investment Dollars Proportion | $`\phi`$_$`\mu`$_g | Proportion of investment dollars to ground inflow | unitless | 0 | $`[0,1]`$ |

***Initial Conditions***
| Full Variable Name | Model Variable Name | Definition | Units |Default Value | Allowable Range |
| -------------- | ------ | ---------- | ----- |------------- | --------------- |
| Population Fill | p_0 | Proportion of population to carrying capacity | Unitless | 0.65 | $`[0, \infty)`$ |
| Actual Per-Capita Demand | $`\chi`$_0 | per-capita demand, accounting for short-term conservation (proportion of mean inflow) | Unitless | $8E-7$ | $`[0, \infty)`$ |
| Surface Storage Fill | $`\upsilon`$_s_0 | Fill proportion of surface storage | Unitless | 1 | $`[0, 1]`$ |
| Ground Storage Fill | $`\upsilon`$_g_0 | Fill proportion of ground storage | Unitless | 1 | $`[0, 1]`$ |
| Surface Storage Capacity | $`\upsilon`$_bar_s_0 | Surface storage capacity (multiple of inflow standard deviation) | Unitless | 4 | $`[0, \infty)`$ |
| Ground Storage Capacity | $`\upsilon`$_bar_g_0 | Ground storage capacity (multiple of inflow standard deviation) | Unitless | 1 | $`[0, \infty)`$ |
| Surface Processing Capacity | w_s_0 | Surface processing capacity (proportion of storage capacity & mean inflow) | Unitless | 1 | $`[0, \infty)`$ |
| Ground Processing Capacity | w_g_0 | Ground processing capacity (proportion of storage capacity & mean inflow) | Unitless | 0 | $`[0, \infty)`$ |
| Surface Inflow | q_s_0 | Surface inflow (proportion of mean) | Unitless | 1 | $`[0, \infty)`$ |
| Ground Inflow | q_g_0 | Ground inflow (proportion of mean) | Unitless | 0 | $`[0, \infty)`$ |
| Mean Surface Inflow | $`\mu`$_s_0 | Mean surface inflow | AFY or Bgal/yr | 200000 | $`[0, \infty)`$ |
| Mean Ground Inflow | $`\mu`$_g_0 | Mean ground inflow | AFY or Bgal/yr | 0 | $`[0, \infty)`$ |
| Surface Inflow Variation | C_v_s_0 | Surface flow coefficient of variation | Unitless | 0.1 | $`[0, \infty)`$ |
| Ground Inflow Variation | C_v_g_0 | Ground flow coefficient of variation | Unitless | 0.01 | $`[0, \infty)`$ |
| Delivery Efficiency | $`\eta`$_0 | Delivery efficiency (proportion of withdrawn delivered) | Unitless | 1 | $`[0, \infty)`$ |
| Base Per-Capita Demand | $`\chi`$bar_0 | Base per-capita demand, independent of ST conservation (proportion of mean inflow) | Unitless | $8E-7$ | $`[0, \infty)`$ |
| Per-Capita Revenue | f_0 | Proportion of per-capita revenue to max ($\pi$_max) | Unitless | 0.5 | $`[0,1]`$ | 
| Average Bond Investment | J_b_avg_0 | Average annual bond-sourced investment | $/yr | 69000000 | $`[0,\infty)`$ |

### PMA Cities Setup Changable Parameters and Initial Conditions

***Parameters***

| Full Parameter Name | Model Variable Name | Definition | Units |Default Value | Allowable Range |
| -------------- | ------ | ---------- | ----- |------------- | --------------- |
| Mean SRP Inflow | $`\mu`$_SRP | Mean inflow into PMA through SRP canals | AFY | 900000 | $`[0,\infty)`$ |
| Mean CAP Inflow | $`\mu`$_CAP | Mean inflow into PMA through CAP canals | AFY | 650491 | $`[0,\infty)`$ |
| Mean Groundwater Inflow | $`\mu`$_g_0 | Mean Groundwater Inflow into PMA | AFY | 690602 | $`[0,\infty)`$ |
| Max Per Capita Revenue | $`\pi`$_max | Max Per Capita Revenue | \$/yr | 1000 | $`[0,\infty)`$ |
| Carrying Capacity | $`\kappa`$ | Population carrying capacity | Unitless | 1500000 | $`(0, \infty)`$ |
| Intrinsic Growth Rate | r | logistically fit, intrinsic growth rate of population | Unitless | 0.1 | $(0,\infty)$ |
| Max Delivery Efficiency | $`\eta`$_max | Max attainable delivery efficiency | Unitless | 2.0 | $`[0, 2]`$ |
| Max Legal Use of Ground Storage | a_gv | Max legally allowed use of ground storage (proportion) | Unitless | 1 | $`[0,1]`$ |
| Max Legal Use of Ground Inflow | a_gq | Max legally allowed use of ground inflow (proportion) | Unitless | 1 | $`[0,1]`$ |
| Max Legal Use of Surface Storage | a_sv | Max legally allowed use of surface storage (proportion) | Unitless | 0 | $`[0,1]`$ |
| SRP Allocation | a_SRP | Proportion of SRP inflow allocated to city | Unitless | 0.30883 | $`[0,1]`$ | 
| CAP Allocation | a_CAP | Proportion of CAP inflow allocated to city | Unitless | 0.41168 | $`[0,1]`$ | 
| CAP Low Priority Allocation | a_CAP_low | Proportion of Low Priority CAP inflow allocated to city | Unitless | 0.5324 | $`[0,1]`$ |
| CAP High Priority Allocation | a_CAP_high | Proportion of High Priority CAP inflow allocated to city | Unitless | 0.3894 | $`[0,1]`$ |
| SRP NCS & Gatewater Allocation | a_SRP_NG | Proportion of SRP inflow allocation attributable to NCS and gatewater | Unitless | 0.06367 | $`[0,1]`$ |
| Min Per Capita Use | d_min | Min possible per-capita demand | (Bgal/yr OR AFY)/person | 0.04 | $`[0, \infty)`$ |
| Hard Infrastructure Decay | $`\delta`$ | Hard infrastructure decay value | Unitless | 0.05 | $`[0,1]`$ | 
| Background Conservation Rate | $`\delta`$_dbar | Annual decay rate of long-term conservation measures (background) | Unitless | 0.0003 | $`[0,1]`$ |
| LT Demand Mgmt Implementation Time | $`\tau`$_d | Time to implement long-term demand management | yrs | 1 | $`[1,\infty)`$ |
| Hard Infrastructure Implementation Time | $`\tau`$_i | Hard infrastructure implementation times | yrs | 3 | $`[1,\infty)`$ | 
| Long-Term Goal Supply Buffer | $`\gamma`$_2 | Goal long-term proportion between supply and demand | Unitless | 1.2 | $`[0,\infty)`$ |
| Minimum Debt Service Coverage Ratio | $`\gamma`$_3 | Minimum Debt Service Coverage Ratio Allowed in any Year | unitless | 2 | $`[0,\infty)`$ |
| Max Rate Increase | $`\psi`$_r | Max possible proportional increase in rates | Unitless | 0.06 | $`[0,\infty)`$ |
| Delivery Efficiency Priority | $`\beta`$_$`\eta`$ | Proportion of long-term investments directed to delivery efficiency | Unitless | 0.2 | $`[0,1]`$ |
| Surface Processing Priority | $`\beta`$_w_s | Proportion of long-term investments directed to surface processing capacity | Unitless | 0 | $`[0,1]`$ |
| Ground Processing Priority | $`\beta`$_w_g | Proportion of long-term investments directed to ground processing capacity | Unitless | 0.7 | $`[0,1]`$ |
| Short-Term Curtailment Sensitivity | $`\lambda`$_1 | Sensitivity (inst. friction component) in short-term investment | Unitless | 22 | $`[0, \infty)`$ |
| Long-Term Investment Sensitivity | $`\lambda`$_2 | Sensitivity (inst. friction component) in long-term investment | Unitless | 22 | $`[0, \infty)`$ |
| Rate-Setting Sensitivity | $`\lambda`$_3 | Sensitivity (inst. friction component) in rate-setting | Unitless | 22 | $`[0, \infty)`$ |
| ST Curtailment Activation Threshold | $`\epsilon`$_1 | Threshold for action (inst. friction component) in short-term investment | Unitless | 0 | $`[0, \infty)`$ |
| LT Investment Activation Threshold | $`\epsilon`$_2| Threshold for action (inst. friction component) in long-term investment | Unitless | 0 | $`[0, \infty)`$ |
| Rate-Setting Activation Threshold | $`\epsilon`$_3| Threshold for action (inst. friction component) in rate-setting | Unitless | 0 | $`[0, \infty)`$ |
| Operating Cost Function Coefficient | g_o | Coefficient in operating costs function (see equation) | Dollars/(persons*Vol) | 0.1435 |  $`[0, \infty)`$ |
| LT Dem Mgmt Investment Coefficient | g_dbar | Coefficient in LT dem mgmt investment function | AFY/Dollars | 5948 | $`[0,\infty)`$ |
| Operating Cost Population Scale Factor | z_op | Population scale factor in operating cost function | unitless | 0.5581 | $`[0,\infty)`$ |
| Operating Cost Demand Scale Factor | z_od | Demand scale factor in operating cost function | unitless | 1.0303 | $`[0,\infty)`$ |
| LT Dem Mgmt Investment Scale Factor | z_dbar | Investment scale factor for long-term demand management | unitless | 1 | $`[0,\infty)`$ |
| Del Eff Investment Scale Factor | z_$`\eta`$ | Investment scale factor for delivery efficiency | unitless | 1.01266 | $`[0,\infty)`$ |
| Proc Capacity Investment Scale Factor | z_w | Investment scale factor for processing capacity | unitless | 1.04019 | $`[0,\infty)`$ |
| Mean Surface Inflow Change Type | $`\Delta\mu`$_s_type | Type of change to surface mean inflow (0 = none, 1 = gradual, 2 = sudden) | Unitless | 2 | $`\{0,1,2\}`$ |
| Mean Surface Inflow Percent Change | $`\Delta\mu`$_s_pc | Percent change to surface mean inflow | Unitless | -0.284 | $`[0, \infty)`$ |
| Mean Sufrace Inflow Change Time | $`\Delta\mu`$_s_t | Time that sudden change occurs or time that gradual change will be complete | yrs | 14 | $`[0, \infty)`$ |
| Base Groundwater Use | $`\theta`$_g | Proportion of Demand that is usually served by groundwater | Unitless | 0.024 | $`[0,1]`$ |
| SRP Demand Proportion | $`\theta`$_1 | Proportion of Demand that is SRP eligible | Unitless | 0.5 | $`[0,1]`$ |
| Projection Years | $`\tau`$_p | Years of Projection for Long-Term Investment Signal | Years | 5 | $`[0, \infty)`$ |
| ST Dem Mgmt Investment Effectiveness Coefficient | $`\alpha`$ | Coefficient for effectiveness in ST conservation investments | unitless | 0.5 | $`[0,\infty)`$ |
| Bond Life | $`\tau`$_b | Life of issued bonds | years | 15 | $`[0,\infty)`$ |
| Bond Interest Rate | i_b | Interest Rate of Issued Bonds | unitless | 0.04 | $`[0,1]`$ |
| Proportion of Rates from Fixed Charges | $`\beta`$_p | Proportion of the expected revenue to come from fixed per user charges | unitless | 0.5 | $`[0,1]`$ |
| Del Eff Investment Dollars Proportion | $`\phi`$_$`\eta`$ | Proportion of investment dollars to delivery efficiency | unitless | 0.6605 | $`[0,1]`$ |
| Surface Proc Capac Investment Dollars Proportion | $`\phi`$_w_s | Proportion of investment dollars to surface processing capacity | unitless | 0.2964 | $`[0,1]`$ |
| Ground Proc Capac Investment Dollars Proportion | $`\phi`$_w_g | Proportion of investment dollars to ground processing capacity | unitless | 0.0331 | $`[0,1]`$ |
| Initial SRP Availability | A_SRP_0 | Initial Volume of Available SRP Water | AFY | 200275.18 | $`[0,\infty)`$ |

***Initial Conditions***

| Full Variable Name | Model Variable Name | Definition | Units |Default Value | Allowable Range |
| -------------- | ------ | ---------- | ----- |------------- | --------------- |
| Population Fill | p_0 | Proportion of population to carrying capacity | Unitless | 0.86466 | $[0, \infty)$ |
| Actual Per-Capita Demand | $`\chi`$_0 | per-capita demand, accounting for short-term conservation (proportion of mean inflow) | Unitless | 8.96274E-8 | $`[0, \infty)`$ |
| Ground Storage Fill | $`\upsilon`$_g_0 | Fill proportion of ground storage | Unitless | 0.53257 | $`[0, 1]`$ |
| Ground Processing Capacity | w_g_0 | Ground processing capacity (proportion of storage capacity & mean inflow) | Unitless | 0.005315 | $`[0, \infty)`$ |
| Surface Inflow Variation | C_v_s_0 | Surface flow coefficient of variation | Unitless | 0.001 | $[0, \infty)$ |
| Ground Inflow Variation | C_v_g_0 | Ground flow coefficient of variation | Unitless | 0.001 | $[0, \infty)$ |
| Delivery Efficiency | $`\eta`$_0 | Delivery efficiency (proportion of withdrawn delivered) | Unitless | 0.9742 | $[0, \infty)$ |
| Base Per-Capita Demand | $`\chi`$bar_0 | Base per-capita demand, independent of ST conservation (proportion of mean inflow) | Unitless | 8.96274E-8 | $`[0, \infty)`$ |
| Per-Capita Revenue | f_0 | Proportion of per-capita revenue to max ($\pi$_max) | Unitless | 0.23836 | $`[0,1]`$ | 
| Average Bond Investment | J_b_avg_0 | Average annual bond-sourced investment | $/yr | 69694375 | $`[0,\infty)`$ |

## Model Outputs

The run_UWIIM() function will generate three outputs in a general array. 

1. *Dataframe* of time series for state and auxilary variables over the model run

2. *Non-dimensional Time Series Plots* organized as a vector of 13 plots 

i. Shortage Plot (Supply, Demand, Shortage)

ii. Flows Plot (Inflows & Use) 

iii. Storage Plot (Fill Volume & Storage Capacity) 

iv. Error 

v. Attention

vi. Financial Flows (Revenue, Costs, Investments)

vii. Per-Capita Revenue 

viii. Delivery Efficiency 

ix. Per-Capita Demand (Base & Actual)

x. Processing Capacity

xi. Population

xii. Mean Inflow (including augmentation)

xiii. All plots aggregated into one figure

3. *Dimensional Time Series Plots* (using the units chosen in the function inputs) organized as a vector of 13 plots (same as above)

# Replicating the Sensitivity Analysis

To replicate the sensitivity analysis performed in the manuscript, use the "UWIIM_PMA_sensanalysis.ipynb" Jupyter notebook or the "UWIIM_PMA_sensanalysis.jl" Julia file in any Julia compiler. The script will output results files into the "results/raw" folder in both a Julia readable object and a CSV file for each parameter variation tested. 

Additionally, we provide the R scripts used to analyze the raw sensitivity analysis folders in the "results/analysis_scripts" folder. Refer to the "WRR_Figs" file for the script used to generate manuscript figures. All other scripts contain background functions to clean the data and generate the figures as specified in the "WRR_Figs" main script. The figures can be found in the "results/figs" folder. 
