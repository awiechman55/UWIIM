***Parameters (Default Case)***

| Full Parameter Name | Model Variable Name | Definition | Units |Default Value | Allowable Range |
| -------------- | ------ | ---------- | ----- |------------- | --------------- |
| Max Surface Streamflow Mean | $\mu$_s_max | Max mean surface inflow that the city can seek | Bgal/yr OR AFY | 200000 | $[0,\infty)$ |
| Max Ground Streamflow Mean | $\mu$_g_max | Max mean ground inflow (recharge) that the city can seek | Bgal/yr OR AFY | 0 | $[0,\infty)$ |
| Max Per-Capita Revenue | $\pi$_max | Max revenue that city can extract per citizen | Dollars/person/yr | 400 | $[0, \infty)$ |
| Streamflow Variation Scenario Type | $\Delta$C_v_scen | Variation Change Type (0 = no change, 1 = gradual, 2 = sudden) | Unitless | 0 | $\{0,1,2\}$ |
| Surface Streamflow Auto-correlation | $\rho_s$ | 1-year lagged auto-correlation in surface streamflow | Unitless | 0.2 | $[0,1]$ |
| Ground Inflow Auto-correlation | $\rho_g$ | 1-year lagged auto-correlation in ground inflow | Unitless | 0 | $[0,1]$ |
| Carrying Capacity | $\kappa$ | Population carrying capacity | Unitless | 1500000 | $(0, \infty)$ |
| Intrinsic Growth Rate | r | logistically fit, intrinsic growth rate of population | Unitless | 0.1 | $(0,\infty)$ |
| Max Delivery Efficiency | $\eta$_max | Max attainable delivery efficiency | Unitless | 2.0 | $[0, 2]$ |
| Max Surface Storage Capacity | $\bar{\upsilon}$_max_s | Max feasible surface storage capacity (multiple of inflow standard deviation) | Unitless | 12 | $[0, \infty)$ |
| Max Ground Storage Capacity | $\bar{\upsilon}$_max_g | Max feasible ground storage capacity (multiple of inflow standard deviation) | Unitless | 40 | $[0, \infty)$ |
| Max Legal Use of Ground Storage | a_gv | Max legally allowed use of ground storage (proportion) | Unitless | 1 | $[0,1]$ |
| Max Legal Use of Ground Inflow | a_gq | Max legally allowed use of ground inflow (proportion) | Unitless | 1 | $[0,1]$ |
| Max Legal Use of Surface Storage | a_sv | Max legally allowed use of surface storage (proportion) | Unitless | 1 | $[0,1]$ |
| Max Legal Use of Surface Inflow | a_sq | Max legally allowed use of surface inflow (proportion) | Unitless | 1 | $[0,1]$ |
| Min Per Capita Use | d_min | Min possible per-capita demand | (Bgal/yr OR AFY)/person | 0.04 | $[0, \infty)$ |
| Background Conservation Rate | $\delta$_dbar | Annual decay rate of long-term conservation measures (background) | Unitless | 0.0116 | $[0,1]$ |
| Surface Storage Decay Rate | $\delta$_v | Annual decay rate of surface storage capacity | Unitless | 0.001 | $[0,1]$ |
| Delivery Efficiency Decay Rate | $\delta$_$\eta$ | Annual decay rate of delivery efficiency | Unitless | 0.001 | $[0,1]$ |
| Surface Processing Decay Rate | $\delta$_w_s | Annual decay rate of surface processing capacity | Unitless | 0.001 | $[0,1]$ |
| Ground Processing Decay Rate | $\delta$_w_g | Annual decay rate of ground processing capacity | Unitless | 0.001 | $[0,1]$ |
| LT Demand Mgmt Implementation Time | $\tau$_d | Time to implement long-term demand management | yrs | 3 | $[1,\infty)$ |
| Delivery Efficiency Impl Time | $\tau$_$\eta$ | Time to implement delievery efficiency improvements | yrs | 4 | $[1,\infty)$ |
| Storage Capacity Impl Time | $\tau$_v | Time to implement surface storage capacity improvements | yrs | 5 | $[1,\infty)$ |
| Surface Processing Impl Time | $\tau$_w_s | Time to implement surface processing capacity improvements | yrs | 3 | $[1,\infty)$ |
| Ground Processing Impl Time | $\tau$_w_g | Time to implement ground processing capacity improvements | yrs | 3 | $[1,\infty)$ |
| Surface Flow Augmentation Impl Time | $\tau$_$\mu$_s | Time to implement surface flow augmentation improvements | yrs | 4 | $[1,\infty)$ |
| Ground Flow Augmentation Impl Time | $\tau$_$\mu$_g | Time to implement ground flow augmentation improvements | yrs | 4 | $[1,\infty)$ |
| Short-Term Goal Supply Sufficiency | $\gamma$_1 |  Goal short-term proportion between supply and demand | Unitless | 1 | $[0,\infty)$ |
| Long-Term Goal Supply Buffer | $\gamma$_2 | Goal long-term proportion between supply and demand | Unitless | 1.2 | $[0,\infty)$ |
| Minimum Debt Service Coverage Ratio | $\gamma$_3 | Minimum Debt Service Coverage Ratio Allowed in any Year | unitless | 2 | $[0,\infty)$ |
| Max Rate Increase | $\psi$_r | Max possible proportional increase in rates | Unitless | 0.06 | $[0,\infty)$ |
| Delivery Efficiency Priority | $\beta$_$\eta$ | Proportion of long-term investments directed to delivery efficiency | Unitless | 0.2 | $[0,1]$ |
| Storage Capacity Priority | $\beta$_v | Proportion of long-term investments directed to surface storage capacity | Unitless | 0.4 | $[0,1]$ |
| Surface Processing Priority | $\beta$_w_s | Proportion of long-term investments directed to surface processing capacity | Unitless | 0.3 | $[0,1]$ |
| Ground Processing Priority | $\beta$_w_g | Proportion of long-term investments directed to ground processing capacity | Unitless | 0 | $[0,1]$ |
| Surface Augmentation Priority | $\beta$_$\mu$_s | Proportion of long-term investments directed to surface flow augmentation | Unitless | 0 | $[0,1]$ |
| Ground Augmentation Priority | $\beta$_$\mu$_g | Proportion of long-term investments directed to ground flow augmentation | Unitless | 0 | $[0,1]$ |
| Short-Term Curtailment Sensitivity | $\lambda$_1 | Sensitivity (inst. friction component) in short-term investment | Unitless | 22 | $[0, \infty)$ |
| Long-Term Investment Sensitivity | $\lambda$_2 | Sensitivity (inst. friction component) in long-term investment | Unitless | 22 | $[0, \infty)$ |
| Rate-Setting Sensitivity | $\lambda$_3 | Sensitivity (inst. friction component) in rate-setting | Unitless | 22 | $[0, \infty)$ |
| ST Curtailment Activation Threshold | $\epsilon$_1 | Threshold for action (inst. friction component) in short-term investment | Unitless | 0 | $[0, \infty)$ |
| LT Investment Activation Threshold | $\epsilon$_2| Threshold for action (inst. friction component) in long-term investment | Unitless | 0 | $[0, \infty)$ |
| Rate-Setting Activation Threshold | $\epsilon$_3| Threshold for action (inst. friction component) in rate-setting | Unitless | 0 | $[0, \infty)$ |
| Operating Cost Function Coefficient | g_o | Coefficient in operating costs function (see equation) | Dollars/(persons*Vol) | 2.5 |  $[0, \infty)$ |
| LT Dem Mgmt Investment Coefficient | g_dbar | Coefficient in LT dem mgmt investment function | AFY/Dollars | 4.4 E-7 | $[0,\infty)$ |
| Del Eff Investment Coefficient | g_$\eta$ | Coefficient in delivery efficiency investment funciton | AFY/Dollars | 0.0033 | $[0,\infty)$ | 
| Storage Capacity Investment Coefficient | g_vbar | Coefficient in storage capacity investment function | AF/Dollars | 0.003 | $[0, \infty)$ | 
| SW Proc Capacity Investment Coefficient | g_ws | Coefficient in surface processing capacity investment function | AFY/Dollars | 0.00015 | $[0,\infty)$ |
| GW Proc Capacity Investment Coefficient | g_wg | Coefficient in ground processing capacity investment function | AFY/Dollars | 0.00015 | $[0,\infty)$ |
| SW Inflow Aug Investment Coefficient | g_$\mu$s | Coefficient in surface inflow augmentation investment function | AFY/Dollars | 0.0001 | $[0,\infty)$ |
| GW Inflow Aug Investment Coefficient | g_$\mu$g | Coefficient in ground inflow augmentation investment function | AFY/Dollars | 0.0001 | $[0,\infty)$ |
| Operating Cost Population Scale Factor | z_op | Population scale factor in operating cost function | unitless | 0.563 | $[0,\infty)$ |
| Operating Cost Demand Scale Factor | z_od | Demand scale factor in operating cost function | unitless | 0.831 | $[0,\infty)$ |
| LT Dem Mgmt Investment Scale Factor | z_dbar | Investment scale factor for long-term demand management | unitless | 1 | $[0,\infty)$ |
| Del Eff Investment Scale Factor | z_$\eta$ | Investment scale factor for delivery efficiency | unitless | 0.82 | $[0,\infty)$ |
| Storage Capacity Investment Scale Factor | z_vbar | Investment scale factor for storage capacity | unitless | 1 | $[0,\infty)$ |
| SW Proc Capacity Investment Scale Factor | z_ws | Investment scale factor for processing capacity | unitless | 1 | $[0,\infty)$ |
| GW Proc Capacity Investment Scale Factor | z_wg | Investment scale factor for processing capacity | unitless | 1 | $[0,\infty)$ |
| SW Inflow Aug Investment Scale Factor | z_$\mu$s | Investment scale factor for sw inflow augmentation | unitless | 1 | $[0,\infty)$ |
| GW Inflow Aug Investment Scale Factor | z_$\mu$g | Investment scale factor for gw inflow augmentation | unitless | 1 | $[0,\infty)$ |
| Case Specification | case | Indicate which case type to follow (0 = default, 1 = PHX) | Unitless | 0 |$\{0,1\}$ |
| Mean Surface Inflow Change Type | $\Delta\mu$_s_type | Type of change to surface mean inflow (0 = none, 1 = gradual, 2 = sudden) | Unitless | 0 | $\{0,1,2\}$ |
| Mean Surface Inflow Percent Change | $\Delta\mu$_s_pc | Percent change to surface mean inflow | Unitless | 0 | $[0, \infty)$ |
| Mean Sufrace Inflow Change Time | $\Delta\mu$_s_t | Time that sudden change occurs or time that gradual change will be complete | yrs | 0 | $[0, \infty)$ |
| Max Surface Processing | w_max_s | Max attainable surface processing (proportion of storage capacity + mean inflow) | Unitless | 1 | $[0,1]$ |
| Max Ground Processing | w_max_g | Max attainable ground processing (proportion of storage capacity + mean inflow) | Unitless | 0.5 | $[0,1]$ |
| Base Groundwater Use | $\theta$_g | Proportion of Demand that is usually served by groundwater | Unitless | 0.0 | $[0,1]$ |
| Projection Years | $\tau$_p | Years of Projection for Long-Term Investment Signal | Years | 5 | $[0, \infty)$ |
| ST Dem Mgmt Investment Effectiveness Coefficient | $\alpha$ | Coefficient for effectiveness in ST conservation investments | unitless | 0.5 | $[0,\infty)$ |
| Bond Life | $\tau$_b | Life of issued bonds | years | 15 | $[0,\infty)$ |
| Bond Interest Rate | i_b | Interest Rate of Issued Bonds | unitless | 0.04 | $[0,1]$ |
| Proportion of Rates from Fixed Charges | $\beta$_p | Proportion of the expected revenue to come from fixed per user charges | unitless | 0.5 | $[0,1]$ |
| Del Eff Investment Dollars Proportion | $\phi$_$\eta$ | Proportion of investment dollars to delivery efficiency | unitless | 0.6 | $[0,1]$ |
| Stor Capac Investment Dollars Proportion | $\phi$_v | Proportion of investment dollars to storage capacity | unitless | 0.3 | $[0,1]$ |
| Surface Proc Capac Investment Dollars Proportion | $\phi$_w_s | Proportion of investment dollars to surface processing capacity | unitless | 0.3 | $[0,1]$ |
| Ground Proc Capac Investment Dollars Proportion | $\phi$_w_g | Proportion of investment dollars to ground processing capacity | unitless | 0 | $[0,1]$ |
| Surface Inflow Investment Dollars Proportion | $\phi$_$\mu$_s | Proportion of investment dollars to surface inflow | unitless | 0 | $[0,1]$ |
| Ground Inflow Investment Dollars Proportion | $\phi$_$\mu$_g | Proportion of investment dollars to ground inflow | unitless | 0 | $[0,1]$ |