# Urban Water Infrastructure Investment Model (UWIIM), PMA Version

This repository contains the source code for the model and sensitivity analysis configuration, raw outputs of the sensitivity analysis, R scripts used to analyze the outputs of sensitivity analysis, and the figures used in the published manuscript associated with this version of the Urban Water Infrastructure Investment Model (UWIIM). The model is a general coupled infrastructure dynamical systems model, written in the Julia programming language (version 1.8.4) as a discrete time dynamical system, that incorporates operational infrastructure and political-economic feedback processes governing urban water systems. The model has been parameterized to reflect the water resources, financial, and demand context of three cities in the Phoenix Metropolitan Area (PMA): Phoenix, Scottsdale, and Queen Creek. 

# Running the Model

The UWIIM can be run through the provided Jupyter notebook ("UWIIM_PMA.ipynb") or any Julia compiler with the "UWIIM_PMA.jl" file. Both files contain all background functions with heavily commented explanations that can also be found in the Supporting Information document (see "docs" folder).

After compiling the source code in a Julia compiler or the Jupyter notebook, the one-line function that runs the UWIIM is **run_UWIIM()**, which takes four possible inputs. 

## Model Inputs

1. *setup*: Setup function containing the parameters and initial conditions for the model run. This version offers four setup functions: **Default()**, **Phoenix()**, **Scottsdale()**, **QueenCreek()**. If the setup functions are called without any specified arguments, they will output the default parameter and initial conditions corresponding to that city as laid out in the Supporting Information document. However, a model user may change any of the parameter or initial conditions in the arguments of the setup functions (see docs/variableTables). 
2. *t_run*: number of years (time steps) to run the model. The default is 50.
3. *year_0*: initial year of the model run. This is primarily used for the time series plots. The default is 2010.
4. *units*: this specifies the water volume units to use for plotting. The default is acre-feet, "AF," which should be used **for all PMA cities.** Any other unit specification will create dimensional analysis errors in the parameters and initial conditions. However, in the Default setup, one can also use "gal" for this argument to use gallon-based water volume units. 

## Model Outputs

The run_UWIIM() function will generate three outputs in a general array. 

1. *Dataframe* of time series for state and auxilary variables over the model run (see docs/variableTables)

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
