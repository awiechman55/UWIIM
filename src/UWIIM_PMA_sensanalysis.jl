#   Read Me
#   ≡≡≡≡≡≡≡≡≡

#   This is a notebook that contains the code used to conduct the sensitivity
#   analyzes described in the attached Supporting Information document for the
#   Phoenix Metropolitan Area cities.

#   1. Pre-Sensitivity Analysis Setup
#   ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

#   1.1 Sensitivity Analysis Specific Packages
#   ============================================

using SharedArrays #for storing sensitivity analysis results
using JLD2 #for saving objects 
using FlexiMaps #for log ranges

#   1.2 Load Model (UWIIM)
#   ========================

include("UWIIM_PMA.jl")

#   2. Sensitivity Analysis Functions
#   ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

#   2.1 Create Parameter Sets Given Range
#   =======================================

#   1 Variable
#   ––––––––––––

#This function creates an array where each row is a parameter vector that will be tested in the sensitivity analysis of 1 parameter. 
function createPSet_1var(default_p, p1_range, p1_index, p1_matrix, p1_matrix_index)
    p_var = p1_range
    num_traj = size(p_var)[1]
    pset = Array{Any}(nothing,num_traj,length(default_p))
    
    for i in 1:num_traj
        for j in 1:length(default_p)
            if j == p1_index
                if(p1_matrix)
                    pset[i,j] = reshape(p_var[i,:],1,length(p_var[i,:]))
                else
                    pset[i,j] = p_var[i]
                end
            else
                pset[i,j] = copy(default_p[j])
            end
        end
    end
    
    return pset
end;

#   2 Variables
#   –––––––––––––

#This function creates an array where each row is a parameter vector that will be tested in the sensitivity analysis of 2 parameters. 
function createPSet_2var(default_p,p1_range,p2_range,p1_index,p2_index, p1_matrix, p1_matrix_index, p2_matrix, p2_matrix_index, add_Δμ, Δμ_range)
    Δμ_def = copy(default_p[25])
    if(p1_matrix)
        p1_def = copy(default_p[p1_index])
        if(p2_matrix)
            p2_def = copy(default_p[p2_index])
            if(add_Δμ)
                p_var = [[i,j,k] for k in eachrow(Δμ_range), j in eachrow(p1_range), i in eachrow(p2_range)]
            else
                p_var = [[i,j] for j in eachrow(p1_range), i in eachrow(p2_range)]
            end
        else
            if(add_Δμ)
                p_var = [[i,j,k] for k in eachrow(Δμ_range), j in eachrow(p1_range), i in p2_range]
            else
                p_var = [[i,j] for j in eachrow(p1_range), i in p2_range]
            end
        end
    else
        if(p2_matrix)
            p2_def = copy(default_p[p2_index])
            if(add_Δμ)
                p_var = [[i,j,k] for k in eachrow(Δμ_range), j in p1_range, i in eachrow(p2_range)]
            else
                p_var = [[i,j] for j in p1_range, i in eachrow(p2_range)]
            end
        else
            if(add_Δμ)
                p_var = [[i,j,k] for k in eachrow(Δμ_range), j in p1_range, i in p2_range]
            else
                p_var = [[i, j] for j in p1_range, i in p2_range]
            end
        end
    end
    
    num_traj = length(p_var)
    
    pset = Array{Any}(nothing,num_traj,length(default_p))
    
    for i in 1:num_traj
        for j in 1:length(default_p)
            if j == p1_index
                if p1_matrix
                    temp = copy(p_var[i][2])
                    if j == p2_index
                        temp[p2_matrix_index] = p_var[i][1][p2_matrix_index]  
                    end   
                    
                    pset[i,j] = temp
                else
                    pset[i,j] = copy(p_var[i][2])
                end
            elseif j == p2_index
                if p2_matrix
                    if j != p1_index
                        pset[i,j] = copy(p_var[i][1])
                    end
                else
                    pset[i,j] = copy(p_var[i][1])
                end
            elseif (add_Δμ && j == 25)
                pset[i,j] = copy(p_var[i][3])
            else
                pset[i,j] = copy(default_p[j])
            end
        end
    end
    
    return pset
    
end;

#   2.2 Run Parameter Sets
#   ========================

#   Adjust g (if necessary)
#   –––––––––––––––––––––––––

#This function adjusts the dependent/auxilary parameters g and J_b to reflect the financial assumptions (see Supporting Information) given a new parameter setting 
function adjust_gandJ_b(p_i, x_0_i, QC)
    p_def = copy(p_i)
    x_0_def = copy(x_0_i)
    
    #Set Common Parameters
    μ_s_0 = p_def[1][1]
    μ_g_0 = p_def[1][2]
    μ_tot = μ_s_0 + μ_g_0 
    μ_SRP = p_def[1][3]
    μ_CAP = p_def[1][4]
    a_SRP = p_def[13][4]
    a_CAP = p_def[13][5]
    a_SRP_NCS = p_def[13][8]
    a_gq = p_def[13][2]
    a_gv = p_def[13][1]
    θ_1 = p_def[27][2]
    V_bar_g_0 = 200*a_gq*μ_g_0
    V_g_0 = x_0_def[4]*V_bar_g_0
    A_l_g_0 = a_gq*μ_g_0 + a_gv*max(0,V_g_0 - (100*a_gq*μ_g_0))
    κ = p_def[5]
    p_0 = x_0_def[1]
    χbar_0 = x_0_def[16]
    χ_0= x_0_def[2]
    f_0 = x_0_def[17]
    π_bar = p_def[2]
    g_o=p_def[22][1]
    z_op=p_def[23][1]
    z_od=p_def[23][2]
    τ_b = p_def[31]
    i_b = p_def[32]
    J_b_avg_0=x_0_def[18]*((1+15*0.04)/(1+τ_b*i_b)) #edit to reflect new bond life and interest rate
    η_0=x_0_def[15]
    w_g_0=x_0_def[8]
    w_s_0=x_0_def[7]
    δ_η=p_def[15][3]
    δ_w_s=p_def[15][5]
    δ_w_g=p_def[15][6]
    z_η=p_def[23][4]
    z_w=p_def[23][6]
    ϕ_η=p_def[35][1]
    ϕ_ws=p_def[35][3]
    ϕ_wg=p_def[35][4]
    
    if(QC==1)
        p_0_QC=74842/κ #2016
        χbar_0_QC = 16344.01/(μ_tot*p_0_QC*κ) #2016
        χ_0_QC=χbar_0_QC
        f_0_QC=23734654/(p_0_QC*κ*π_bar) #2016
        J_b_avg_0_QC=3492917.5*((1+15*0.04)/(1+τ_b*i_b)) #2016 and account for new bond life and interest rate
        A_SRP_0_QC=1
        
        g_new=set_g(μ_s_0,μ_g_0,μ_tot,μ_SRP,μ_CAP,a_SRP,a_CAP,a_SRP_NCS,a_gq,a_gv,θ_1,V_g_0, A_l_g_0, V_bar_g_0,0,
            κ,p_0_QC,χbar_0_QC,χ_0_QC,f_0_QC,π_bar,g_o,z_op,z_od,J_b_avg_0_QC,τ_b,i_b,η_0,w_g_0,w_s_0,δ_η,δ_w_s,
            δ_w_g,z_η, z_w, ϕ_η, (ϕ_ws+ϕ_wg), A_SRP_0_QC)
    else
        if(p_def[27][2]==0.5) #If Phoenix
            A_SRP_0 = 200275.18
        else #If Scottsdale
            A_SRP_0 = 13632.80
        end
        
        g_new=set_g(μ_s_0,μ_g_0,μ_tot,μ_SRP,μ_CAP,a_SRP,a_CAP,a_SRP_NCS,a_gq,a_gv,θ_1,V_g_0, A_l_g_0, V_bar_g_0,0,
            κ,p_0,χbar_0,χ_0,f_0,π_bar,g_o,z_op,z_od,J_b_avg_0,τ_b,i_b,η_0,w_g_0,w_s_0,δ_η,δ_w_s,
            δ_w_g,z_η, z_w, ϕ_η, (ϕ_ws+ϕ_wg), A_SRP_0)
    end
    
    return [g_new J_b_avg_0]
end;

#   Results Array to Dataframe
#   ––––––––––––––––––––––––––––

#This function converts the results array into a dataframe
function rezDF(rez_array,p1_name,p2_name)
    if(p2_name == "nothing") #if it is just a 1 parameter variation sensitivity analysis
        rez = DataFrame(p1_name=rez_array[:,2,1], reliab_total=vec(mean(rez_array[:,4,:],dims=2)), reliab_avg=vec(mean(rez_array[:,5,:],dims=2)),reliab_min=vec(mean(rez_array[:,6,:],dims=2)), 
            d_end = vec(mean(rez_array[:,7,:], dims=2)), π_end = vec(mean(rez_array[:,8,:],dims=2)),Δμ_pc = vec(mean(rez_array[:,9,:],dims=2)))
    else #if it is a 2 parameter variation sensitivity analysis
        rez = DataFrame(p1_name=rez_array[:,2,1], p2_name = rez_array[:,3,1], reliab_total=vec(mean(rez_array[:,4,:],dims=2)), reliab_avg=vec(mean(rez_array[:,5,:],dims=2)),
            reliab_min=vec(mean(rez_array[:,6,:],dims=2)), d_end = vec(mean(rez_array[:,7,:], dims=2)), π_end = vec(mean(rez_array[:,8,:],dims=2)),Δμ_pc = vec(mean(rez_array[:,9,:],dims=2)))
        p2_symbol=Symbol(p2_name)
        rename!(rez,:p2_name => p2_symbol)
    end
    
    p1_symbol=Symbol(p1_name)
    rename!(rez,:p1_name => p1_symbol)
    
    return rez
end;

#   Gather Results to Array
#   –––––––––––––––––––––––––

#This function conducts the sensitivity analysis on the given parameter set by running the parameter set for the given amount of time and recording the summary statistics (Supporting Information)
function runPSet(pset, p1_name, p2_name, p1_index, p2_index, p1_matrix, p1_matrix_index, p2_matrix, p2_matrix_index, x_0, num_t, num_it, QC, f=f)
    num_traj = size(pset)[1]
    num_p = size(pset)[2]
    
    #Create empty results tuple to store results of each trajectory
    rez_array = SharedArray{Float64}(num_traj,9,num_it)
     
    #Generate and Store Trajectories
    for i in 1:num_traj
        ####Set up Parameters
        p_i = copy(pset[i,:])
        
        ####Set Initial Conditions Given Parameters
        x_0_i = copy(x_0[1:18])
        #Add Planned Investment Memory
        for k in 1:length(p_i[16])
            if(p_i[16][k]>1)
                for i in 1:p_i[16][k]-1
                    append!(x_0_i,0.0)
                end
            end
        end
        
        # Fix g and J_b_avg if Relevant (if varying τ_b, i_b, δ)
        p_newgJb = [15, 31, 32]
        if(p1_index in p_newgJb || p2_index in p_newgJb)
            gJb_new = adjust_gandJ_b(p_i, x_0_i, QC)
            p_i[22][3] = gJb_new[1]
            p_i[22][5] = gJb_new[2]
            p_i[22][6] = gJb_new[3]
        end
        
        ######Define System
        DS = DiscreteDynamicalSystem(f, x_0_i, p_i)
        
        ######Run Model for Defined System
        for j in 1:num_it #loop for multiple iterations of same p set 
            #Generate Trajectory and Store It
            tr_j = trajectory(DS, num_t)
            
            ##Note the Iteration Number
            rez_array[i,1,j] = j
            
            ##Store the Varied Parameters Used
            #p1
            if(p1_matrix)
                rez_array[i,2,j] = copy(p_i[p1_index][p1_matrix_index])
            else
                rez_array[i,2,j] =copy(p_i[p1_index])
            end
            
            #p2
            if(p2_index==nothing)
                
            elseif(p2_index == 0)
                rez_array[i,3,j] = 0
            else
                if(p2_matrix)
                    rez_array[i,3,j] = copy(p_i[p2_index][p2_matrix_index])
                else
                    rez_array[i,3,j] = copy(p_i[p2_index])
                end
            end
            
            #Δμ
            rez_array[i,9,j] = copy(p_i[25][2])
            
            ##Record State Variables for each year in iteration
            D_bar_ij = zeros(num_t+1) #Create empty matrix to record total baseline demand
            S_ij = zeros(num_t+1) #Create empty matrix to record supply
            reliab_ij = zeros(num_t+1)
            d_ij = zeros(num_t+1) #Create empty matrix to record PC demand 
            π_ij = zeros(num_t+1) #Create empty matrix to record PC rates
            
            for t in 1:(num_t+1) #Loop Over All Time Measurements in a Trajectory
                D_bar_ij[t] = D_bar(tr_j[t],p_i,t-1) # Record total baseline demand at t
                S_ij[t] = S(tr_j[t],p_i,t-1) #Record total supply at t
                reliab_ij[t] = min(1,S_ij[t]/D_bar_ij[t])
                d_ij[t] = (D_bar_ij[t]/P(tr_j[t],p_i,t-1))*892.15
                π_ij[t] = tr_j[t][17]*p_i[2]
            end
            
            #####Record Reliability Metrics
            reliab_total = sum(reliab_ij)
            rez_array[i,4,j] = reliab_total
            reliab_avg = reliab_total/(num_t+1) #includes initial condition
            rez_array[i,5,j] = reliab_avg
            reliab_min = minimum(reliab_ij)
            rez_array[i,6,j] = reliab_min
             
            #####Record Normalized Demand, Resiliency, Storage Capacity Metrics
            d_end = d_ij[end]
            rez_array[i,7,j] = d_end
            
            #####Record Ending PC Rate
            π_end = π_ij[end]
            rez_array[i,8,j] = π_end
        end 
    end
    
    
    ##Convert Results into a Data Frame
    rez = rezDF(rez_array, p1_name, p2_name)
    
    return rez
end;

#   2.3 Run Sensitivity Analysis and Save Results
#   ===============================================

function runSA(p1_index,p1_matrix,p1_matrixindex,p1_saindex,p2_index,p2_matrix,p2_matrixindex,p2_saindex,cities_name,cities_x0,cities_p,file_preamble,num_t,num_it,add_Δμ)
    ##run the sensitivity analysis on each city 
    for c in 1:length(cities_name)
        c_p = copy(cities_p[c])
        c_x0 = copy(cities_x0[c])
        c_name = cities_name[c]
        
        print(c_name*" ")
    
        ##different protocols for 1d vs 2d sensitivity analysis
        if(p2_index==nothing)
            ##Create pset
            pset = createPSet_1var(c_p,p_ranges[p1_saindex],p1_index,p1_matrix,p1_matrixindex)
            
            ##run iteractions and generate results
            if(c_name == "QC")
                rez=runPSet(pset, p_names[p1_saindex],"nothing", p1_index, nothing, p1_matrix, p1_matrixindex, nothing, nothing, c_x0, num_t, num_it,1)
            else
                rez=runPSet(pset, p_names[p1_saindex],"nothing", p1_index, nothing, p1_matrix, p1_matrixindex, nothing, nothing, c_x0, num_t, num_it,0)
            end
            
            ##Record file name for saving 
            file_name = file_preamble*p_names[p1_saindex]*"_"*c_name
            file_name_rez = file_name*"_rez"
            file_name_csv = file_name*".csv"
        else
            ##Create pset
            pset = createPSet_2var(c_p,p_ranges[p1_saindex],p_ranges[p2_saindex],p1_index,p2_index,p1_matrix,p1_matrixindex,p2_matrix,p2_matrixindex,add_Δμ,copy(p_ranges[1]))
            
            ##run iteractions and generate results
            if(c_name == "QC")
                rez=runPSet(pset, p_names[p1_saindex],p_names[p2_saindex], p1_index, p2_index, p1_matrix, p1_matrixindex, p2_matrix, p2_matrixindex, c_x0, num_t, num_it,1)
            else
                rez=runPSet(pset, p_names[p1_saindex],p_names[p2_saindex], p1_index, p2_index, p1_matrix, p1_matrixindex, p2_matrix, p2_matrixindex, c_x0, num_t, num_it,0)
            end
            
            ##Record file name for saving 
            file_name = file_preamble*p_names[p1_saindex]*"_"*p_names[p2_saindex]*"_"*c_name
            if(add_Δμ)
                file_name = file_name*"_3d"
            end
            
            file_name_rez = file_name*"_rez"
            file_name_csv = file_name*".csv"
        end
        
        ##save results    
        save_object(file_name_rez,rez)
        CSV.write(file_name_csv,rez)
    end
end;

#   3. Define Sensitivity Analysis Variable Ranges (PMA Case)
#   ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

#   PMA City Defaults
#   ===================

PHX_setup = copy(Phoenix())
PHX_p = copy(PHX_setup[1])
PHX_x_0 = copy(Float64.(PHX_setup[2]))

Sc_setup = copy(Scottsdale())
Sc_p = copy(Sc_setup[1])
Sc_x_0 = copy(Float64.(Sc_setup[2]))

QC_setup = copy(QueenCreek())
QC_p = copy(QC_setup[1])
QC_x_0 = copy(Float64.(QC_setup[2]));

PMA_names = ["PHX", "Sc", "QC"]
PMA_p = [PHX_p, Sc_p, QC_p]
PMA_x0 = [PHX_x_0, Sc_x_0, QC_x_0];

#   Parameter Ranges (Same Across Cities)
#   =======================================

#   Combine Ranges into Vectors
#   –––––––––––––––––––––––––––––

num_p = length(PHX_p)
p_ranges = Array{Any,1}(nothing,num_p)
p_names = Array{Any,1}(nothing,num_p)
p_ranges_spec = Array{Any,1}(nothing,num_p);

#   Suddent CAP Change Magnitude (\Delta\mu)
#   ––––––––––––––––––––––––––––––––––––––––––

Δμ_s_range = collect(-1.0:0.01:0.0)

Δμ_range = zeros(length(Δμ_s_range), length(PHX_p[25])) 
for i in 1:length(Δμ_s_range)
    Δμ_range[i,:] = copy(PHX_p[25])
    Δμ_range[i,2] = Δμ_s_range[i]
end;

p_ranges[1] = Δμ_range
p_names[1] = "MagCAPShock"
p_ranges_spec[1] = Δμ_s_range;

#   Projection Year Range (\tau^p)
#   ––––––––––––––––––––––––––––––––

τ_p_range = collect(0:1:10)

p_ranges[6] = τ_p_range
p_names[6] = "ProjYears"
p_ranges_spec[6] = τ_p_range;

#   Institutional Friction Sensitivity (\lambda)
#   ––––––––––––––––––––––––––––––––––––––––––––––

λ_1 = copy(PHX_p[20][1])
λ_2 = copy(PHX_p[20][2])
λ_3 = copy(PHX_p[20][3])

λ_1_range = collect(maprange(log,4,220,length=50))
λ_2_range = collect(maprange(log,4,220,length=50))
λ_3_range = collect(maprange(log,4,220,length=50))

λ_range_1 = zeros(length(λ_1_range),3)
λ_range_2 = zeros(length(λ_2_range),3)
λ_range_3 = zeros(length(λ_3_range),3)

for i in 1:length(λ_1_range)
    λ_range_1[i,:] = [λ_1_range[i] λ_2 λ_3]
    λ_range_2[i,:] = [λ_2 λ_2_range[i] λ_3]
    λ_range_3[i,:] = [λ_3 λ_2 λ_3_range[i]]
end

p_ranges[7] = λ_range_1
p_names[7] = "IF_sens_s"
p_ranges_spec[7] = λ_1_range;
p_ranges[8] = λ_range_2
p_names[8] = "IF_sens_l"
p_ranges_spec[8] = λ_2_range;
p_ranges[9] = λ_range_3
p_names[9] = "IF_sens_r"
p_ranges_spec[9] = λ_3_range;

#   Institutional Friction Activation Threshold (\epsilon)
#   ––––––––––––––––––––––––––––––––––––––––––––––––––––––––

ϵ_1 = copy(PHX_p[21][1])
ϵ_2 = copy(PHX_p[21][2])
ϵ_3 = copy(PHX_p[21][3])

ϵ_1_range = collect(0.0:0.01:0.5)
ϵ_2_range = collect(0.0:0.01:0.5)
ϵ_3_range = collect(0.0:0.02:1)

ϵ_range_1 = zeros(length(ϵ_1_range),3)
ϵ_range_2 = zeros(length(ϵ_2_range),3)
ϵ_range_3 = zeros(length(ϵ_3_range),3)

for i in 1:length(ϵ_1_range)
    ϵ_range_1[i,:] = [ϵ_1_range[i] ϵ_2 ϵ_3]
    ϵ_range_2[i,:] = [ϵ_1 ϵ_2_range[i] ϵ_3]
end

for i in 1:length(ϵ_3_range)
    ϵ_range_3[i,:] = [ϵ_1 ϵ_2 ϵ_3_range[i]]
end

p_ranges[10] = ϵ_range_1
p_names[10] = "IF_thresh_s"
p_ranges_spec[10] = ϵ_1_range;
p_ranges[11] = ϵ_range_2
p_names[11] = "IF_thresh_l"
p_ranges_spec[11] = ϵ_2_range;
p_ranges[12] = ϵ_range_3
p_names[12] = "IF_thresh_r"
p_ranges_spec[12] = ϵ_3_range;

#   Goal Safety Factor and Debt Service Coverage Ratio (\gamma_2 & \gamma_3)
#   ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

γ_1 = copy(PHX_p[17][1])
γ_2 = copy(PHX_p[17][2])
γ_3 = copy(PHX_p[17][3])
γ_2_range = collect(1:0.01:1.3);
γ_3_range = collect(1.5:0.01:2.5)

γ_range_2 = zeros(length(γ_2_range),3)
γ_range_3 = zeros(length(γ_3_range),3)

for i in 1:length(γ_2_range)
    γ_range_2[i,:] = [γ_1 γ_2_range[i] γ_3]
end

for i in 1:length(γ_3_range)
    γ_range_3[i,:] = [γ_1 γ_2 γ_3_range[i]]
end

p_ranges[13] = γ_range_2
p_names[13] = "Goal_l"
p_ranges_spec[13] = γ_2_range;

p_ranges[21] = γ_range_3
p_names[21] = "DSCR"
p_ranges_spec[21] = γ_3_range;

#   Max Rate Increase Proportion (\psi_r)
#   –––––––––––––––––––––––––––––––––––––––

ψ_s = copy(PHX_p[18][1])
ψ_r = copy(PHX_p[18][2])
ψ_r_range = collect(0.01:0.01:0.5)
ψ_range_r = zeros(length(ψ_r_range),2)

for i in 1:length(ψ_r_range)
    ψ_range_r[i,:] = [ψ_s ψ_r_range[i]]
end

p_ranges[15] = ψ_range_r
p_names[15] = "Choice_r"
p_ranges_spec[15] = ψ_r_range;

#   Background Conservation Rate (\delta_{\bar{d}})
#   –––––––––––––––––––––––––––––––––––––––––––––––––

δ_dbar_range = collect(0:0.00005:0.001);
δ_d = copy(PHX_p[15][1])
δ_v = copy(PHX_p[15][3])
δ_η = copy(PHX_p[15][4])
δ_w_s = copy(PHX_p[15][5])
δ_w_g = copy(PHX_p[15][6])
δ_range_dbar = zeros(length(δ_dbar_range),6)

for i in 1:length(δ_dbar_range)
    δ_range_dbar[i,:]=[δ_d δ_dbar_range[i] δ_v δ_η δ_w_s δ_w_g]
end;

p_ranges[16] = δ_range_dbar
p_names[16] = "Decay_dbar"
p_ranges_spec[16] = δ_dbar_range;

#   Hard Infrastructure Decay Rate (\delta)
#   –––––––––––––––––––––––––––––––––––––––––

δ_range_hard = collect(0.025:0.005:0.075);

δ_vbar = copy(PHX_p[15][4])
δ_dbar = copy(PHX_p[15][2])
δ_d = copy(PHX_p[15][1])

δ_range = zeros(length(δ_range_hard),6)

for i in 1:length(δ_range_hard)
    δ_range[i,:] = [δ_d δ_dbar δ_range_hard[i] δ_vbar δ_range_hard[i] δ_range_hard[i]]
end
    
p_ranges[17] = δ_range
p_names[17] = "Decay_hard"
p_ranges_spec[17] = δ_range_hard;

#   Time to Implementation (\tau^i)
#   –––––––––––––––––––––––––––––––––

τ_i_range = collect(1:1:5);
τ_i_def = copy(PHX_p[16])

τ_range_i = zeros(length(τ_i_range),7)

for i in 1:length(τ_i_range)
    τ_range_i[i,:] = [τ_i_def[1] τ_i_range[i] τ_i_def[3] τ_i_range[i] τ_i_range[i] τ_i_def[6] τ_i_def[7]]
end

p_ranges[18] = τ_range_i
p_names[18] = "Impl_Time"
p_ranges_spec[18] = τ_i_range;

#   Proportion of Rates from Fixed Costs (\beta_p)
#   ––––––––––––––––––––––––––––––––––––––––––––––––

β_p_range = collect(0:0.01:1)

p_ranges[4] = β_p_range
p_names[4] = "FixedPropRates"
p_ranges_spec[4] = β_p_range;

#   Short-Term Conservation Sensitivity (\alpha)
#   ––––––––––––––––––––––––––––––––––––––––––––––

α_range = collect(0.01:0.01:0.95);

p_ranges[19] = α_range
p_names[19] = "Conserv_Sens"
p_ranges_spec[19] = α_range;

#   Bond Length (\tau_b)
#   ––––––––––––––––––––––

τ_b_range = collect(5:1:30)

p_ranges[20] = τ_b_range
p_names[20] = "Bond_Time"
p_ranges_spec[20] = τ_b_range;

#   Interest Rate (i_b)
#   –––––––––––––––––––––

i_b_range = collect(0:0.01:0.1)

p_ranges[22] = i_b_range
p_names[22] = "IntRate"
p_ranges_spec[22] = i_b_range;

#   4. Run Sensitivity Analysis Experiments
#   ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

#   List of Sensitivity Analysis Runs
#   ===================================

###### Globaal Parameters
file_preamble = "./results/raw/"

num_t = 50
num_it = 1

###### 1D Experiments - 30 minutes
#1: Magnitude of Sudden CAP Shock
#2: Max Rate Increase 
#3: Long-Term Safety Buffer Goal 
#4: IF_sens_s
#5: IF_sens_l 
#6: IF_sens_r
#7: IF_thresh_s
#8: IF_thresh_l 
#9: IF_thresh_r 
#10: Background Conservation
#11: Hard Infrastructure Decay
#12: Implementation Time 
#13: Conservation Sensitivity 
#14: Goal Debt Service Ratio
#15: Bond Life 
#16: Interest Rate
#17: Projection Years
#18: Fixed Proportion of Rates

runs_1d_index = [25,18,17,20,20,20,21,21,21,15,15,16,30,17,31,32,29,36]
runs_1d_matrix = [true,true,true,true,true,true,true,true,true,true,true,true,false,true,false,false,false,false]
runs_1d_matrixindex = [2,2,2,1,2,3,1,2,3,2,3,2,nothing,3,nothing,nothing,nothing,nothing]
runs_1d_saindex = [1,15,13,7,8,9,10,11,12,16,17,18,19,21,20,22,6,4]

#for i in 1:length(runs_1d_index)
#    print(p_names[runs_1d_saindex[i]]*" run time: ")
#    @time runSA(runs_1d_index[i],runs_1d_matrix[i],runs_1d_matrixindex[i],runs_1d_saindex[i],nothing,nothing,
#        nothing,nothing,PMA_names,PMA_x0,PMA_p,file_preamble,num_t,num_it,false)
#end

###### 2D Experiments
#Primary Sensitivity Analysis
##1: Magnitude of CAP Shock vs. IF ST Invest Sensitivity 
##2: Magnitude of CAP Shock vs. IF LT Invest Sensitivity
##3: Magnitude of CAP Shock vs. IF Rate Sensitivity 
##4: Magnitude of CAP Shock vs. IF ST Thresh
##5: Magnitude of CAP Shock vs. IF LT Thresh 
##6: Magnitude of CAP Shock vs. IF Rate Thresh 
##7: Institutional Friction (ST) 
##8: Institutional Friction (Long-Term) 
##9: Institutional Friction (Rates) 

#Additional Sensitivity Analysis
##10: Magnitude of CAP Shock vs. Infrastructure Decay
##11: Magnitude of CAP Shock vs. Projection Years
##12: Magnitude of CAP Shock vs. Implementation Years
##13: Magnitude of CAP Shock vs. Bond Life
##14: Magnitude of CAP Shock vs. Interest Rate
##15: Magnitude of CAP Shock vs. Max Rate Increase
##16: Magnitude of CAP Shock vs. Fixed Proportion of Rates
##17: Magnitude of CAP Shock vs. DSCR
##18: Magnitude of CAP Shock vs. Conservation Sensitivity
##19: Magnitude of CAP Shock vs. Background Demand Decay
##20: Magnitude of CAP Shock vs. Long-Term Safety Buffer Goal

runs_2d_index_1 = [25,25,25,25,25,25,20,20,20,25,25,25,25,25,25,25,25,25,25,25]
runs_2d_matrix_1 = [true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true]
runs_2d_matrixindex_1 = [2,2,2,2,2,2,1,2,3,2,2,2,2,2,2,2,2,2,2,2]
runs_2d_saindex_1 = [1,1,1,1,1,1,7,8,9,1,1,1,1,1,1,1,1,1,1,1]
runs_2d_index_2 = [20,20,20,21,21,21,21,21,21,15,29,16,31,32,18,36,17,30,15,17]
runs_2d_matrix_2 = [true,true,true,true,true,true,true,true,true,true,false,true,false,false,true,false,true,false,true,true]
runs_2d_matrixindex_2 = [1,2,3,1,2,3,1,2,3,3,nothing,2,nothing,nothing,2,nothing,3,nothing,2,2]
runs_2d_saindex_2 = [7,8,9,10,11,12,10,11,12,17,6,18,20,22,15,4,21,19,16,13];

for i in 20:length(runs_2d_index_1)
    print(p_names[runs_2d_saindex_1[i]]*"-"*p_names[runs_2d_saindex_2[i]]*" run time: ")
    @time runSA(runs_2d_index_1[i],runs_2d_matrix_1[i],runs_2d_matrixindex_1[i],runs_2d_saindex_1[i],runs_2d_index_2[i],
        runs_2d_matrix_2[i],runs_2d_matrixindex_2[i],runs_2d_saindex_2[i],PMA_names,PMA_x0,PMA_p,file_preamble,num_t,num_it,false)
end;

#   Report Meta Information on Sensitivity Analysis
#   =================================================

#   Load Results
#   ––––––––––––––

function load_rez(c,v,n_dim;file_preamble="./results/raw/")
    if(n_dim==1)
        file_name = file_preamble*p_names[v]*"_"*PMA_names[c]*"_rez"
    else
        file_name = file_preamble*p_names[v[1]]*"_"*p_names[v[2]]*"_"*PMA_names[c]*"_rez"
    end
    
    rez = load(file_name,"single_stored_object")
    
    return rez
end;

#   Save Global Tables
#   ––––––––––––––––––––

#This function records meta information on the sensitivity analysis, including city names, varied parameters names, and metrics recorded
function write_p_info(num_t,num_it;p_names=p_names,p_ranges_spec=p_ranges_spec,file_preamble = "./results/")
    p_names_file_name = file_preamble*"p_names.csv"
    cities_file_name = file_preamble*"cities.csv"
    metrics_file_name = file_preamble*"metrics.csv"
    
    p_names_cut=[]
    for p in 1:length(p_names)
        if(p_names[p]==nothing)
            push!(p_names_cut,"N/A")
        else
            push!(p_names_cut,p_names[p])
        end
    end
    
    p_names_df = DataFrame(p=p_names_cut)
    cities_df=DataFrame(City=PMA_names)
    
    
    test=load_rez(1,1,1)
    metrics = names(test)[2:end]
    metrics_df = DataFrame(Metric=metrics)
    
    CSV.write(p_names_file_name, p_names_df)
    CSV.write(cities_file_name, cities_df)
    CSV.write(metrics_file_name, metrics_df)
end;

write_p_info(50,1)