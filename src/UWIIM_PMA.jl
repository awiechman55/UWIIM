#   Read Me
#   ≡≡≡≡≡≡≡≡≡

#   This published version of the Urban Water Infrastructure Investment Model
#   (UWIIM) contains the code (Julia) used to run the model for a default city
#   (general city with a representative surface and groundwater soure) or the
#   three Phoenix Metropolitan Area (PMA) cities. The UWIIM is a discrete time
#   dynamical systems model with annual time steps.
# 
#   Section 1 contains the dynamical system definition. Each function is
#   accompanied by markdown equation defintions. For additional explanation, see
#   the Supporting Information from the attached manuscript. All state variables
#   and their dynamic equations are aggregated together in the "Equations of
#   Motion" in section 1.4.
# 
#   Section 2 contains information for how the model is parameterized and
#   initialized, including the setup for the general default city, and PMA
#   city-specific setups. Each setup is defined in terms of a function that
#   takes as its input, alterations to the originally intended parameters and
#   initial conditions (as specified in the Supporting Information document).
# 
#   Section 3 contains information for how to run the model, including the
#   single function that allows one to run a model instance with a certain
#   defined setup (see Section 2) and receive an output for the given model
#   instance. This output contains a (1) dataframe of all state and multiple
#   auxiliary variables over time, (2) multiple time series plots of the
#   variables in their reduced dimensional form, and (3) multiple time series
#   plots of the variables in their common dimensional form (i.e. AFY)
# 
#   Section 4 contains example outputs for the default PMA setups. A user can
#   alter the parameter and initial condition settings in the setup function to
#   compare outputs.
# 
#   One notational difference between the code and the manuscript is the
#   reliance on non-dimensional forms of the state variables in the Equations of
#   Motion. The non-dimensional forms allow us to relate all variables to the
#   units of mean inflow (e.g., AFY), carrying capacity, and maximum per-capita
#   annual revenue ($/(person*yr)). Those non-dimensional forms are the
#   following:
# 
#   Demand: \chi_t = \frac{d_t}{\mu_t}
# 
#   Population: p_t = \frac{P_t}{\kappa_t}
# 
#   Storage Volume: \upsilon_{i,t} = \frac{V_{i,t}}{\bar{V}_{i,t}}
# 
#   Storage Capacity: \bar{\upsilon}_{i,t} =
#   \frac{\bar{V}_{i,t}}{c^v_{i,t}\mu_{i,t}}
# 
#   Processing Capacity: w_{i,t} = \frac{A^w_{i,t}}{\bar{V}_{i,t} + \mu_{i,t}}
# 
#   Streamflow: q_{i,t} = \frac{Q_{i,t}}{\mu_{i,t}}
# 
#   Per-Capita Annual Revenue: f_t = \frac{\hat{\pi}_t}{\bar{\pi}_t}

#   Packages Used
#   ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

using Pkg
using Plots; pyplot()
using DynamicalSystems
using Distributions
using DataFrames
using ForwardDiff
using CSV

#   1. Model Definition
#   ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

#   1.1 Definitions of Auxilary Operational System Variables
#   ==========================================================

#   Mean Inflow
#   –––––––––––––

# :\
# 
#   mut = \mu^st + \mu^g_t $

μ(x,p,t) = x[11]+x[12];

#   Streamflow, Inflow
#   ––––––––––––––––––––

#   Total Inflow (Q)
#   ------------------

# :q_t
# 
#   = \frac{Qt}{\mut} $

function q(x,p,t)
    return Q(x,p,t)/μ(x,p,t)
end;

# :Q_t
# 
#   = q^st\mu^st + q^gt\mu^gt $

function Q(x,p,t)
    return Q_s(x,p,t) + Q_g(x,p,t)
end;

# :Q
# 
#   ^st = q^st \mu^s_t $

function Q_s(x,p,t)
    return x[9]*x[11]
end;

# :Q
# 
#   ^gt = q^gt \mu^g_t $

function Q_g(x,p,t)
    return x[10]*x[12]
end;

#   Available Inflows (Q^a)
#   -------------------------

#   Surface Available Inflows (Q^{a,s})
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅

# :Q
# 
#   ^{a,s} = \begin{cases} Q^{a,SRP} + Q^{a,CAP} & \text{if} \quad cases = 1 \
#   a^{s,q}Q^s_t & \text{otherwise} \end{cases} $

function Q_a_s(x,p,t) 
    if(p[24]==1)
        return Q_a_SRP(x,p,t) + Q_a_CAP(x,p,t)
    else
        return p[13][4]*Q_s(x,p,t)
    end
end;

#   CAP and SRP Available Inflows (Q^{a,CAP} & Q^{a,SRP})
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅

# :\
# 
#   tilde{Q}^{CAP,L}t = max(0,70022-(\bar{\mu}2 - (Q^st - \mu1))) $
# 
# :\
# 
#   tilde{Q}^{CAP,H}t = (Q^st - \mu1) - \tilde{Q}^{CAP,L}t $
# 
# :\
# 
#   tilde{Q}^{CAP}t = a^{q,2H}\tilde{Q}^{CAP,H}t + a^{q,2L}\tilde{Q}^{CAP,L}_t $

function Q_a_CAP(x,p,t)
    Q_CAP = Q_s(x,p,t) - p[1][3]
    short = p[1][4] - Q_CAP
    NIA_avail = max(0,70022-short)
    high_avail = Q_CAP - NIA_avail
    
    return p[13][6]*NIA_avail + p[13][7]*high_avail
end;

# :\
# 
#   tilde{Q}^{SRP}t = min(\mu1a^{q,1}, \frac{Dt}{\etat} \xi1 + \mu1 a^{q,NG}) $

function Q_a_SRP(x,p,t)
    max_prop_use = (D_bar(x,p,t)/x[15])*p[27][2]
    NG = p[13][8]*p[1][3]
    max_allocation = p[1][3]*p[13][4]
    
    return min(max_allocation, max_prop_use + NG)  
end;

#   Groundwater Available Inflows
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅

# :Q
# 
#   ^{a,g} = a^{g,q}Q^g_t $

Q_a_g(x,p,t) = p[13][2]*Q_g(x,p,t);

#   Total Available Inflows
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅

# :\
# 
#   tilde{Q}t = a^q Qt $

Q_a(x,p,t) = Q_a_s(x,p,t) + Q_a_g(x,p,t);

# :q
# 
#   ^at = \frac{\tilde{Q}t}{\mu_t} $

q_a(x,p,t) = Q_a(x,p,t)/μ(x,p,t);

#   Population
#   ––––––––––––

#   Total Population (P)
#   ----------------------
# 
# :P_t
# 
#   = p_t \kappa $

P(x,p,t) = x[1]*p[5];

#   Projected Population (P^{proj})
#   ---------------------------------
# 
#   *For Long-Term Investment Controller if Projection Setting is Turned On
# 
# :P
# 
#   ^{proj}t = \frac{Kt}{\frac{Kt-Pt}{Pt}exp(-r\taup) + 1} $

function P_proj(x,p,t) #project population out τ_p years
    P_t = P(x,p,t) #current population 
    κ = copy(p[5]) #Assume same carrying capacity .
    
    return κ/(((κ-P_t)/P_t)*exp(-p[6]*p[29]) + 1)
end;

function P_proj_1(x,p,t) #only project population out 1 year
    P_t = P(x,p,t) #current population 
    κ = copy(p[5]) #Assume same carrying capacity 
    
    return κ/(((κ-P_t)/P_t)*exp(-p[6]) + 1)
end;

#   Demand
#   ––––––––

#   Pre-ST Conservation Demand (D)
#   --------------------------------

# :D_t
# 
#   = dt Pt = \chit \mut P_t $
# 
# :\
# 
#   bar{D}t = \bar{d}t Pt = \bar{\chi}t \mut Pt $

D(x,p,t) = x[2]*P(x,p,t)*μ(x,p,t);

D_bar(x,p,t) = x[16]*P(x,p,t)*μ(x,p,t);

#   Post-ST Conservation Demand (D^{ST})
#   --------------------------------------

#   
# d^{ST}_t = (\chi_t - H^d_t)\mu_t
# $
# 
# *Note, $H^d_t
# 
#   is given in units of \chi (per capita demand/AFY mean inflow)
# 
# :D
# 
#   ^{ST}t = d^{ST}tP_t $

d_ST(x,p,t) = (x[2]-u_1(x,p,t))*μ(x,p,t);

D_ST(x,p,t) = d_ST(x,p,t)*P(x,p,t);

#   Projected Demand (D^{proj})
#   -----------------------------

# :\
# 
#   hat{D}t = \hat{P}t \bar{d}t
#   \left(1-\delta^{\bar{d}}\right)^{\tau^{\bar{d}}p} = \hat{P}t \bar{\chi}t
#   \mut \left(1-\delta^{\bar{d}}\right)^{\tau^{\bar{d}}p} $

D_proj(x,p,t) = P_proj(x,p,t)*x[16]*μ(x,p,t)*(1-p[15][2])^p[29]; #projects demand for τ_p years

D_proj_1(x,p,t) = P_proj_1(x,p,t)*x[16]*μ(x,p,t)*(1-p[15][2]); #only projects demand for 1 year

#   Storage
#   –––––––––

#   Storage Capacity (\bar{V})
#   ----------------------------

# :\
# 
#   bar{V}t = \bar{\upsilon}t c^v{t} \mut $

Vbar_s(x,p,t) = x[5]*x[13]*x[11];

Vbar_g(x,p,t) = x[6]*x[14]*x[12];

#   Volume (V)
#   ------------

# :V_t
# 
#   = \upsilont \bar{V}t $

V_s(x,p,t) = x[3]*Vbar_s(x,p,t);

V_g(x,p,t) = x[4]*Vbar_g(x,p,t);

#   Projected Volume (\upsilon_{proj})
#   ------------------------------------

#   Project Surface Water Volume (\upsilon^{proj,s})
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅
# 
#   Assume no change in surface storage volume over projection period. Note, the
#   PMA cases have near zero surface storage capacity.

function υ_proj_s(x,p,t)
    return x[3]
end;

#   Project Groundwater Volume (\upsilon^{proj,g}) & Delivery Efficiency
#  (\eta^{proj})
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅
# 
#   We project groundwater volume together with delivery efficiency for
#   computational efficiency because the delivery efficiency in each year
#   updates based on available water, which is a function of groundwater storage
#   volume.
# 
#   The projection process follows the steps outlined in the Supporting
#   Information

function υ_and_η_proj_g(x,p,t)
    #Note Current State
    A_l_s_now = A_l_s(x,p,t)
    Inflow = Q_a_g(x,p,t) + Q_b(x,p,t)
    
    #Define variables to hold updated outflows, availability, and delivery efficiency
    A_next = A(x,p,t)
    O_g_next = O_g(x,p,t)
    O_g_prev = copy(O_g_next)
    O_proj = 0
    η_now = copy(x[15])
    η_prev = copy(x[15])
    η_next = copy(x[15])
    new_SW = 0
    new_SW_prev = 0
    O_s_add_next = 0
    O_s_add_prev = 0
    A_w_s_now = A_w_s(x,p,t)
    A_w_s_prev = copy(A_w_s_now)
    A_w_s_next = copy(A_w_s_now)
    A_prev=0
    
    if(p[29]>0)
        #Number of investments that can be implemented in the projection period
        possibleInvests_η = floor(Int,min(p[29],p[16][2]-1)) 
        possibleInvests_w_s = floor(Int,min(p[29],p[16][4]-1)) 
        
        #Loop for Each Projection Year and Calculate Projected Groundwater Use
        for y in 1:p[29]
            #0: Add Groundwater Use to Sum
            O_proj += O_g_next
            
            #1: New Surface Water Processing Availability
            A_w_s_prev = copy(A_w_s_next)
            if(possibleInvests_w_s>0)
                A_w_s_next += x[plan_index_k(x,p,t)[4]+(y-1)]
                possibleInvests_w_s -= 1
            end
            
            #2: New Surface Water Availability
            new_SW_prev = copy(new_SW)
            O_s_add_prev = copy(O_s_add_next)
            if(A_l_s_now>=A_w_s_next)
                new_SW = A_w_s_next - A_w_s_prev
            else
                new_SW=0
            end
            O_s_add_next += new_SW
            
            #3: New Total Availability
            A_prev = copy(A_next)
            A_next += Inflow - O_g_next + new_SW
            
            #4: New Delivery Efficiency
            η_prev = copy(η_next)
            
            if(possibleInvests_η>0)
                u_impl_η=x[plan_index_k(x,p,t)[2]+(y-1)]
                η_next = η_prev + (u_impl_η/A_prev) 
                possibleInvests_η -= 1
                
                η_next = ifelse(η_next>p[11],p[11],η_next)
            else
                η_next = copy(η_prev)
            end
            
            #5: New Groundwater Use
            O_g_prev = copy(O_g_next)
            O_g_next = (η_prev/η_next)*(O_g_prev + O_s_add_prev) - O_s_add_next
        end
        
        υ_proj_t=max(100*p[13][2]*p[1][2],(V_g(x,p,t) + Inflow*p[29] - O_proj))/Vbar_g(x,p,t)
        η_proj_t=η_next
        return [υ_proj_t η_proj_t]
    else
        return [x[4] x[15]]
    end
end;

#   Available Water
#   –––––––––––––––––

#   Total Available (A)
#   ---------------------

# :A_t
# 
#   = A^st + A^gt $

A(x,p,t) = A_s(x,p,t) + A_g(x,p,t);

#   Available by Types of Sources
#   -------------------------------

#   Surface Water Availability
# 
# :A
# 
#   ^st = min(A^{l,s}t + A^{T,s}_t) $

A_s(x,p,t) = min(A_l_s(x,p,t), A_w_s(x,p,t));

#   Groundwater Availability
# 
# :A
# 
#   ^gt = min( A^{l,g}t + A^{T,g}_t) $

A_g(x,p,t) = min(A_l_g(x,p,t), A_w_g(x,p,t));

#   Legally Available Volume for Withdrawal (A^l)
#   -----------------------------------------------

#   How much of the physically present water is the system legally entitled to
#   use. Determined by a proportion of physically available volume (proportional
#   allocation).
# 
#   
# A^{l,s}_t = \tilde{Q}^s_t + a^{s,v}V^s_t
# $
# 
# $
# A^{l,g}_t = a^{g,v}V^g_t + a^{s,q}Q^g_t
# $
# 
# For Phoenix case, available $V^g_t
# 
#   is determined after subtracting the 100 year safe-yield allowance
#   (100a^{q,g}Q^g_t) from V^g_t to determine the available surplus credits. If
#   \beta^c=1 (cheating safe-yield), the city can use all of V^g_t.

function A_l_s(x,p,t)
    Q_a = Q_a_s(x,p,t)
    
    return Q_a + p[13][3]*V_s(x,p,t)
end;

function A_l_g(x,p,t) 
    #Available Groundwater From Stored Water
    if(p[24]==1) #Phoenix Metro Area Case
        A_l_g_v = p[13][1]*max(0,V_g(x,p,t)-p[13][2]*p[1][2]*100) ##only count available V if above 100-yr SY amount 
    else
        A_l_g_v = p[13][1]*V_g(x,p,t)
    end
    
    #Available Groundwater From Annual Inflow
    A_l_g_q = p[13][2]*Q_g(x,p,t);
    
    return A_l_g_v + A_l_g_q
end;

#   Technically Available Volume for Withdrawal (A^w)
#   ---------------------------------------------------

#   The technological limit is steered by processing capacity (i.e., treatment,
#   pumping, etc.). This specifies the limit to what the system can process from
#   its physical inflow.
# 
# :A
# 
#   ^wt = wt (\bar{V} + \mu_t) $

A_w_s(x,p,t) = x[7]*(Vbar_s(x,p,t) + x[11]);

A_w_g(x,p,t) = x[8]*(Vbar_g(x,p,t) + x[12]);

#   Maximum Processing Capacity (\bar{w})
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅

# :\
# 
#   bar{w}{i,t} = \frac{A^l{i,t}}{\mu{i,t} + \bar{V}{i,t}} $

function w_max_s(x,p,t)
    W_max_s = A_l_s(x,p,t)
        
    return W_max_s/(x[11] + Vbar_s(x,p,t))
end;

function w_max_g(x,p,t)
    A_l_g_t = A_l_g(x,p,t)
    w_max_g = (A_l_g_t)/(Vbar_g(x,p,t)+x[12]) #(taken out) 1.6 accounts for the desire to have surplus capacity for peak intrannual demands
    
    return w_max_g
end;

#   Projected Available Volume (A^{proj})
#   ---------------------------------------

#   Projected Infrastructure States
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅

##this function identifies the index of the state variable vector (x) associated with planned investments for each infrastructure type k
function plan_index_k(x,p,t)
    index_k = zeros(7)
    counter = 19
    
    for k in 1:7
        if(p[16][k]>1)
            index_k[k] = copy(counter)
            counter+=p[16][k]-1
        end
    end
    
    return floor.(Int,index_k)
end;

#   Projected Delivery Efficiency (\eta^{proj})
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅

function η_proj(x,p,t)
    return υ_and_η_proj_g(x,p,t)[2]
end;

#   Projected Storage Capacity (\bar{V}^{proj})
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅

function Vbar_proj_s(x,p,t)
    Vbar_t = Vbar_s(x,p,t)
    
    ##Determine the number of planned investments that will be implemented by the end of the projection period
    if(p[29]==0) #no projection (τ_p = 0)
        return Vbar_t
    elseif(p[16][3]==1) #immediate implementation (τ_i = 1)
        return Vbar_t
    else
        possibleInvests = floor(Int,min(p[29],p[16][3]-1)) #calculates the number of possible stored investments that need to be taken into account in the projection
        if(possibleInvests==1)#if only implement 1 planned investment
            return Vbar_t + x[plan_index_k(x,p,t)[3]]
        else
            return Vbar_t + sum(x[plan_index_k(x,p,t)[3]:(plan_index_k(x,p,t)[3]+possibleInvests-1)])
        end
    end
end;

function υbar_proj_s(x,p,t) #convert the projected surface storage capacity to υ non-dimensional units
    return Vbar_proj_s(x,p,t)/(x[13]*x[11])
end;

#   Projected Processing Capacity (w^{proj})
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅

function w_proj_s(x,p,t)
    Vbar_t = Vbar_s(x,p,t)
    Vbar_proj_t = Vbar_proj_s(x,p,t)
    μ_t = copy(x[11])
    μ_proj_t = μ_t
    A_w_t = (x[7]*(Vbar_t+μ_t))
    
    ##Determine the number of planned investments that will be implemented by the end of the projection period
    if(p[29]==0) #no projection (τ_p = 0)
        return x[7]
    elseif(p[16][4]==1) #immediate implementation (τ_i = 1)
        return A_w_t/(Vbar_proj_t + μ_proj_t)
    else
        possibleInvests = floor(Int,min(p[29],p[16][4]-1)) #calculates the number of possible stored investments that need to be taken into account in the projection
        if(possibleInvests==1)#only implement 1 planned investment
            return (A_w_t + x[plan_index_k(x,p,t)[4]])/(Vbar_proj_t + μ_proj_t)
        else
            return (A_w_t + sum(x[plan_index_k(x,p,t)[4]:(plan_index_k(x,p,t)[4]+possibleInvests-1)]))/(Vbar_proj_t + μ_proj_t)
        end
    end
end;

function w_proj_g(x,p,t)
    Vbar_t = Vbar_g(x,p,t)
    Vbar_proj_t = Vbar_t
    μ_t = copy(x[12])
    μ_proj_t = μ_t
    A_w_t = (x[8]*(Vbar_t+μ_t))
    
    ##Determine the number of planned investments that will be implemented by the end of the projection period
    if(p[29]==0) #no projection (τ_p = 0)
        return x[8]
    elseif(p[16][5]==1) #immediate implementation (τ_i = 1)
        return A_w_t/(Vbar_proj_t + μ_proj_t)
    else
        possibleInvests = floor(Int,min(p[29],p[16][5]-1))  #calculates the number of possible stored investments that need to be taken into account in the projection
        if(possibleInvests==1)#only implement 1 planned investment
            return (A_w_t + x[plan_index_k(x,p,t)[5]])/(Vbar_proj_t + μ_proj_t)
        else
            return (A_w_t + sum(x[plan_index_k(x,p,t)[5]:(plan_index_k(x,p,t)[5]+possibleInvests-1)]))/(Vbar_proj_t + μ_proj_t)
        end
    end
end;

function w_max_g_proj(x,p,t)
    A_l_g_t = A_proj_l_g(x,p,t)
    w_max_g = (A_l_g_t)/(Vbar_g(x,p,t)+x[12]) 
    
    return w_max_g
end;

#   Total Projected Available Volume (A^{proj})
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅

#   Same format as A_t but need to project \upsilon. Assume available purchased
#   water remains the same.
# 
# :A
# 
#   ^{proj}t = min(A^{l,proj}t, A^{w,proj}_t) $

A_proj(x,p,t) = A_proj_s(x,p,t) + A_proj_g(x,p,t);

A_proj_s(x,p,t) = min(A_proj_l_s(x,p,t), A_proj_w_s(x,p,t));

A_proj_g(x,p,t) = min(A_proj_l_g(x,p,t), A_proj_w_g(x,p,t));

#   Projected Legally Available Volume (A^{proj,l})
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅

#   Follows the same algorithm as the actual legal availability assuming that
#   planned investments are implemented, the projected groundwater volume, and
#   the same mean inflow

function A_proj_l_s(x,p,t) 
    #Available Surface Water from Storage
    A_proj_l_sv = p[13][3]*υ_proj_s(x,p,t)*Vbar_proj_s(x,p,t)
    
    #Available Surface Water from Inflow
    if(p[24]==1) #Phoenix Metro Area Cases
        #CAP
        Q_CAP = x[11] - p[1][3]
        short = p[1][4] - Q_CAP
        NIA_avail = max(0,70022-short)
        high_avail = Q_CAP - NIA_avail
        CAP_avail = p[13][6]*NIA_avail + p[13][7]*high_avail
        
        #SRP
        max_prop_use = (D_proj(x,p,t)/η_proj(x,p,t))*p[27][2]
        NG = p[13][8]*p[1][3]
        max_allocation = p[1][3]*p[13][4]
        
        SRP_avail = min(max_allocation, max_prop_use + NG) 
        
        return A_proj_l_sv + SRP_avail + CAP_avail
    else
        return  A_proj_l_sv + p[13][4]*x[11]
    end
end; 

function A_proj_l_g(x,p,t) 
    ###Stored Groundwater Availability
    if(p[34]==0)
        if(p[24]==1)
            A_l_g_v = p[13][1]*max(0,υ_and_η_proj_g(x,p,t)[1]*Vbar_g(x,p,t)-p[13][2]*p[1][2]*100) ##only count available V if above 100-yr SY amount 
         else
            A_l_g_v = p[13][1]*υ_and_η_proj_g(x,p,t)[1]*Vbar_g(x,p,t)
        end
    else
        A_l_g_v = υ_and_η_proj_g(x,p,t)[1]*Vbar_g(x,p,t)
    end
    
    ###Groundwater Inflow Availability
    if(p[24]==1)
        A_l_g_q = p[13][2]*x[12]
    else
        A_l_g_q = p[13][2]*x[12];
    end
    
    return A_l_g_v + A_l_g_q
end; 

#   Projected Technically Available Volume (A^{proj,T})
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅

# :A
# 
#   ^{T,proj}t = \hat{w}t(\hat{\bar{V}}t + \mut) $

A_proj_w_s(x,p,t) = w_proj_s(x,p,t)*(Vbar_proj_s(x,p,t)+x[11]);

A_proj_w_g(x,p,t) = w_proj_g(x,p,t)*(Vbar_g(x,p,t)+x[12]);

#   Outflows, Releases
#   ––––––––––––––––––––

#   Total Use/Outflows (O)
#   ------------------------

#   Outflows either are used to satisfy demand (O^d_t) or release flood water
#   (O^f_t)
# 
# :O_t
# 
#   = O^dt + O^ft $

O(x,p,t) = O_d(x,p,t) + O_f(x,p,t);

#   Demand-Related Outflows (O^d)
#   -------------------------------

# :O
# 
#   ^dt = O^st + O^g_t $

O_d(x,p,t) = O_s(x,p,t) + O_g(x,p,t);

#   Surface Water Demand-Related Outflows (O^s)
#   ---------------------------------------------

# :O
# 
#   ^st = \begin{cases} O^1t + O^2t & \text{if} \quad cases = 1 \ min(A^st,
#   \frac{\tilde{D}t}{\etat}(1-\theta^g) & \text{otherwise} \end{cases} $

function O_s(x,p,t)
    if(p[24]==1) #Phoenix Metro Area Case
        return O_1(x,p,t) + O_2(x,p,t) #sum CAP and SRP use
    else 
        need = (D_ST(x,p,t)/x[15])*(1-p[27][1])
        return min(A_s(x,p,t), need)
    end
end;

#   SRP Use (O^1)
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅

# :O
# 
#   ^1t = min(\mu1a^{SRP}, A^{T,s}t, \frac{\tilde{D}t}{\etat}\xi^{SRP} +
#   \mu1a^{NCS}, \frac{\tilde{D}t}{\etat}(1-\xi^{g})) $

function O_1(x,p,t) 
    #calculate need after considering annual groundwater use
    need = (D_ST(x,p,t)/x[15])*(1-p[27][1])
    
    #calculate available water considering legal and technical constraints
    tech_avail = A_w_s(x,p,t)
    allocation = p[1][3]*p[13][4]
    dem_op = (D_ST(x,p,t)/x[15])*p[27][2] + p[13][8]*p[1][3] #on-project demand (demand that can use basic SRP rights) + annual NCS & gatewater 
    
    #return min of need and available
    return min(allocation, tech_avail, dem_op, need)
end;

#   CAP Use (O^2)
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅

# :O
# 
#   ^2t = min(\tilde{Q}^{CAP}t, \frac{\tilde{D}t}{\etat} - O^1t, A^{T,s}t) $

function O_2(x,p,t)
    #calculate need
    need =(D_ST(x,p,t)/x[15])*(1-p[27][1])
    need_left = need - O_1(x,p,t)
    
    #calculate available water considering legal and technical constraints
    legal_avail = Q_a_CAP(x,p,t)
    tech_avail = A_w_s(x,p,t)
    
    return min(legal_avail, tech_avail, need_left)
end;

#   Banked Water (Q^b)
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅

# :Q
# 
#   ^b = \begin{cases} Q^{a,CAP}t - O^2t & \text{if} \quad cases = 1 \ 0 &
#   \text{otherwise} \end{cases} $

function Q_b(x,p,t) 
    if(p[24]==1) #if PMA case
        return Q_a_CAP(x,p,t) - O_2(x,p,t);
    else
        return 0
    end
end; 

#   Groundwater Demand-Related Outflows (O^g_t)
#   ---------------------------------------------

# :O
# 
#   ^gt = min\left(\frac{\tilde{D}t}{\etat} - O^st, A^g_t \right) $

function O_g(x,p,t)
    need = D_ST(x,p,t)/x[15]
    need_left = need - O_s(x,p,t) #after considering surface water use
    
    return min(need_left, A_g(x,p,t))
end;

#   Flood Release Outflows (O^f_t)
#   --------------------------------

# :excess
# 
#   = (V^st + a^{s,q} Q^st - (O^{s}t + Q^bt)) - \bar{V}_t $
# 
# :O
# 
#   ^f_t = \begin{cases} excess & \text{if} \quad excess > 0 \ 0 &
#   \text{otherwise} \end{cases} $

function O_f(x,p,t)
    excess = (V_s(x,p,t) + Q_a_s(x,p,t) - (O_s(x,p,t)+Q_b(x,p,t))) - Vbar_s(x,p,t) 
    
    return max(excess, 0.0)
end;

#   Supply
#   ––––––––

#   Current Supply
#   ----------------

# :S_t
# 
#   = \etat At $

S(x,p,t) = x[15]*A(x,p,t);

#   Projected Supply
#   ------------------

# :\
# 
#   hat{S}t = \hat{\eta}t \hat{A}_t $

S_proj(x,p,t) = η_proj(x,p,t)*A_proj(x,p,t);

#   Finances
#   ––––––––––

#   Revenue
#   ---------

#   Actual Revenue Generated (R)
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅

# :R_t
# 
#   = Pt\hat{\pi}t(\beta^{(\pi)}p +
#   (1-\beta^{(\pi)}p)\frac{\tilde{d}t}{\bar{\chi}t\mu_t}) $

R(x,p,t) = P(x,p,t)*p[2]*x[17]*(p[36]+((1-p[36])*(d_ST(x,p,t)/(x[16]*μ(x,p,t)))));

#   Projected Revenue for Next Year (R^{proj})
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅

# :R
# 
#   ^{proj}t = P^{proj,1}t f_t \bar{\pi} $

R_proj(x,p,t) = P_proj_1(x,p,t)*x[17]*p[2];

#   Operating Costs
#   -----------------

#   Actual Operating Costs (C^o)
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅

# :C
# 
#   ^ot = go Pt^{z{op}}\tilde{D}t^{z{od}} $

C_o(x,p,t) = p[22][1]*(P(x,p,t)^p[23][1])*(D_bar(x,p,t)^p[23][2]);

#   Projected Operating Costs (C^{o,proj})
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅

# :C
# 
#   ^{o,proj}t = go
#   \left(P^{proj,1}t\right)^{z{op}}\left(\tilde{D}^{proj,1}t\right)^{z{od}} $

C_o_proj(x,p,t) = p[22][1]*(P_proj_1(x,p,t)^p[23][1])*(D_proj_1(x,p,t)^p[23][2]);

#   Debt Service
#   --------------

#   Actual Debt Service (C^d)
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅

#   C^d_t = \tilde{J}^b_{t}(1+\tau_bi_b)

C_d(x,p,t) = x[18]*(1+p[31]*p[32]);

#   Needed Debt Service in Next Year (C^{d,need})
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅

# :C
# 
#   ^{d,need}t = \tilde{J}^bt(1+\taubib - \frac{1}{\taub} - ib) +
#   \hat{J}^bt(\frac{1}{\taub} + i_b) $

function C_d_need(x,p,t) #needed debt service for the next year
    return x[18]*(1+p[31]*p[32]-(1/p[31])-p[32]) + J_b(x,p,t)*((1/p[31])+p[32])
end;

#   1.2 Definition of Signal & Error for Each Action Situation
#   ============================================================

#   Signal - Safety Factor for Short-Term Curtailment (M_1)
#   –––––––––––––––––––––––––––––––––––––––––––––––––––––––––

# :M_
# 
#   {1,t} = \frac{St}{Dt} $

M_1(x,p,t) = S(x,p,t)/D(x,p,t);

#   Signal - Projected Safety Fator for Investment (M_2)
#   ––––––––––––––––––––––––––––––––––––––––––––––––––––––

# :M_
# 
#   {2,t} = \frac{S^{proj}t}{D^{proj}t} $

M_2(x,p,t) = S_proj(x,p,t)/D_proj(x,p,t);

#   Signal - Debt Service Capacity Ratio for Rate-Settting (M_3)
#   ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

# :M_
# 
#   {3,t} = \frac{R^{proj}t-C^{o,proj}t}{C^{d,need}_t} $

M_3(x,p,t) = (R_proj(x,p,t) - C_o_proj(x,p,t))/C_d_need(x,p,t);

#   Error - Short-Term (e_1)
#   ––––––––––––––––––––––––––

# :e_
# 
#   {1,t} = \gamma1-M{1,t} $

e_1(x,p,t) = p[17][1] - M_1(x,p,t);

#   Error - Projected for Investment (e_2)
#   ––––––––––––––––––––––––––––––––––––––––

# :e_
# 
#   {2,t} = \gamma2 - M{2,t} $

e_2(x,p,t) = p[17][2] - M_2(x,p,t);

#   Error - Rate-Setting (e_3)
#   ––––––––––––––––––––––––––––

# :e_
# 
#   {3,t} = \gamma3 - M{3,t} $

function e_3(x,p,t)
    return p[17][3] - M_3(x,p,t)
end;

#   1.3 Definition of Controller Steps for Each Action Situation
#   ==============================================================

#   1.3.1 Attention (Y)
#   –––––––––––––––––––––

# :Y_
# 
#   {j,t} = \frac{1}{1 + \mathrm{exp}(-\lambdaj(e{j,t} - \epsilon_j))} $

function Y_1(x,p,t)
    e_1t = e_1(x,p,t)
    
    return 1/(1+exp(-p[20][1]*(e_1t - p[21][1])))
end;

function Y_2(x,p,t)
    e_2t = e_2(x,p,t)
    
    return 1/(1+exp(-p[20][2]*(e_2t - p[21][2])))
end;

function Y_3(x,p,t)
    e_3t = e_3(x,p,t)
    
    return 1/(1+exp(-p[20][3]*(e_3t - p[21][3])))
end;

#   1.3.2 Controller Response: Short-term Curtailment Action Situation (u_1)
#   ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

# :u_
# 
#   {1,t} = \alpha \chit Y{1,t}(1-\frac{\chi{min}}{\chit}) $

function u_1(x,p,t)
    χ_min = p[14]/μ(x,p,t) #convert d_min to χ units
    
    u_1_t = x[2]*Y_1(x,p,t)*p[30]*(1-(χ_min/x[2]))
    
    if (u_1_t > (x[2] - χ_min)) #if investment would push the per-capita demand below the minimum demand, lower the investment
        u_1_t = x[2] - χ_min
    end
    
    return u_1_t
end;

#   1.3.3 Controller Response: Investment Action Situations (u_{2,k})
#   –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

#   1.3.3.1 Maintenance Investment Needs (u^{m,need} & J^{m,need})
#   ----------------------------------------------------------------

#   Total Maintenance Investment Need (in dollars) (J^{m,need})
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅

function J_k_m_need(x,p,t)
    J_k_m = zeros(7)
    
    J_k_m[2] = J_m_η(x,p,t)
    J_k_m[3] = J_m_v(x,p,t)
    J_k_m[4] = J_m_w_s(x,p,t)
    J_k_m[5] = J_m_w_g(x,p,t)
    
    return J_k_m
end;

J_m_need(x,p,t) = sum(J_k_m_need(x,p,t));

#   Total Maintenance Investment Need (in volumetric capacity) (u^{m,need})
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅

u_m_need(x,p,t) = u_m_η(x,p,t) + u_m_v(x,p,t) + u_m_w_s(x,p,t) + u_m_w_g(x,p,t);

#   Maintenance Investment Need for Each Infrastructure Type (u^{m,need}_k &
#  J^{m,need}_k)
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅

#   
# u^{m,need}_k = \left(\tilde{I}_{k,t} - I_{k,t}(1-\delta_k)\right)A_t
# $
# 
# where $\tilde{I}_{k,t}
# 
#   is the goal maintenance state of I_k (if I_{k,t} is above maximum capacity,
#   it is allowed to decay to the maximum, \bar{I}. Otherwise, the goal is the
#   current state)
# 
# :J
# 
#   ^m{k,t} = \begin{cases} g\eta \left(\etat u^{m,need}{\eta,t}\right)^{z\eta}
#   & \text{if delivery efficiency} \ gk \left(u^{m,need}{k,t}\right)^{zk} &
#   \text{otherwise} \end{cases} $

#   Delivery Efficiency

function u_m_η(x,p,t)
    #Note Current State
    x_t = ForwardDiff.value(x[15])
    A_t = ForwardDiff.value(A(x,p,t))
    
    #Decide on Goal State
    x_goal = ifelse(x_t>p[11], p[11], x_t) #If the current capacity is above maximum, allow the capacity to decay to maximum
    
    #Calculate Volumetric Investment Need
    u_m = (x_goal - x_t*(1-p[15][3]))*A_t
    
    #Check and correct for negative investment
    if(u_m < 0)
        return 0
    else
        return u_m
    end
end;

function J_m_η(x,p,t)
    #Note current state and volumetric maintenance investment need
    η_t = ForwardDiff.value(x[15])
    u_m = u_m_η(x,p,t)
    
    #Convert volumetric need to dollars 
    J_m = p[22][3]*((η_t*u_m)^p[23][4])
    
    return J_m
end;

#   Storage Capacity

function u_m_v(x,p,t)
   #Note Current State
    x_t = ForwardDiff.value(x[5])
    c_v_t = ForwardDiff.value(x[13])
    μ_t = ForwardDiff.value(x[11])
    x_max = ForwardDiff.value(p[12][1])
    
    #If it is not an investment priority (like PMA case), no maintenance needed
    if(p[19][2]==0)
        return 0
    else
        #Decide on goal state
        x_goal = ifelse(x_t>x_max,x_max,x_t) #If the current capacity is above maximum, allow the capacity to decay to maximum
    
        #Calculate volumetric investment need
        u_m = (x_goal - x_t*(1-p[15][4]))*μ_t*c_v_t
    
        #check and correct for negative investment
        if(u_m < 0)
            return 0
        else
            return u_m
        end
    
    end
end;

function J_m_v(x,p,t)
    #If it is not an investment priority (like PMA case), no maintenance needed
    if(p[19][2]==0)
        return 0
    else
        #Note volumetric investment need
        u_m = u_m_v(x,p,t)
        
        #Convert volumetric need to dollars
        J_m= p[22][4]*(u_m^p[23][5])
        
        return J_m
    end
end;

#   Surface Processing Capacity

function u_m_w_s(x,p,t)
    #Note Current State
    x_t = ForwardDiff.value(x[7])
    w_max = ForwardDiff.value(w_max_s(x,p,t))
    V_bar_s_t = ForwardDiff.value(Vbar_s(x,p,t))
    μ_s_t = ForwardDiff.value(x[11])
    
    #Decide on goal state
    x_goal = ifelse(x_t>w_max, w_max, x_t) #If the current capacity is above maximum, allow the capacity to decay to maximum
    
    #Calculate volumetric investment need
    u_m = (x_goal - x_t*(1-p[15][5]))*(V_bar_s_t + μ_s_t)
    
    #check and correct for negative investment
    if(u_m < 0)
        return 0
    else
        return u_m
    end
end;

function J_m_w_s(x,p,t)
    #Note volumetric investment need
    u_m = u_m_w_s(x,p,t)
    
    #Convert volumetric need to dollars
    J_m = p[22][5]*(u_m^p[23][6])
    
    return J_m
end;

#   Ground Processing Capacity

function u_m_w_g(x,p,t)
   #Note Current State
    x_t = ForwardDiff.value(x[8])
    w_max = ForwardDiff.value(w_max_g(x,p,t))
    V_bar_g_t = ForwardDiff.value(Vbar_g(x,p,t))
    μ_g_t = ForwardDiff.value(x[12])
    
    #Decide on goal state
    x_goal = ifelse(x_t>w_max, w_max, x_t) #If the current capacity is above maximum, allow the capacity to decay to maximum
    
    #Calculate volumetric investment need
    u_m = (x_goal - x_t*((1-p[15][6])))*(V_bar_g_t + μ_g_t)
    
    #check and correct for negative investment
    if(u_m < 0)
        return 0
    else
        return u_m
    end
end;

function J_m_w_g(x,p,t)
    #Note volumetric investment need
    u_m = u_m_w_g(x,p,t)
    
    #Convert volumetric need to dollars
    J_m = p[22][6]*(u_m^p[23][7])
    
    return J_m
end;

#   1.3.3.2 Expansion Investment Needs (u^{e,need} & J^{e,need})
#   --------------------------------------------------------------

#   Expansionary investment needs are a function of attention, error, and
#   projected demand. We also add a check to ensure that the needs are never
#   negative
# 
# :u
# 
#   ^{e,need}t = max\left(Y{2,t}e{2,t}D^{proj}t, 0\right) $
# 
# :J
# 
#   ^{e,need}{k,t} = \begin{cases} gk \left(\beta{k,t} \etat
#   u^{e,need}{t}\right)^{zk} & \text{if delivery efficiency} \ gk
#   \left(\beta{k,t} u^{e,need}{t}\right)^{zk} & \text{otherwise} \end{cases} $

function u_e_need(x,p,t)
    return max(Y_2(x,p,t)*e_2(x,p,t)*D_proj(x,p,t),0)
end;

function J_k_e_need(x,p,t)
    ##Note Volumetric Expansion Investment Need and current state
    u_e_need_t = ForwardDiff.value(u_e_need(x,p,t))
    η_t = ForwardDiff.value(x[15])
    
    ##Calculate Allocation Percentages for This Year 
    β_k_t = β_k(x,p,t)
    
    ##Initialize J_k vector
    J_k = zeros(length(p[19])+1)
        
    ##Calculate J_k for each k
    for k in 1:length(p[19])
        #error check to ensure that priority or investment need is never negative 
        if(β_k_t[k]<0) 
            print("β_k"*string(k)*"is neg")
        elseif(u_e_need_t<0)
            print("u_e_need is neg")
        end
        
        
        #Account for special functional form of delivery efficiency investment cost function
        if(k==1)
            J_k[k+1] = p[22][2+k]*((β_k_t[k]*η_t*u_e_need_t)^p[23][3+k])
        else
            J_k[k+1] = p[22][2+k]*((β_k_t[k]*u_e_need_t)^p[23][3+k])
        end
    end
    
    ##Calculate remaining investment need for d_bar (long-term demand management)
    J_k[1] = p[22][2]*((((1-sum(β_k_t))*u_e_need_t)/p[17][2])^p[23][3]) ###add consideration for gamma_L
    
    return J_k
end;        

#   1.3.3.3 Re-distribute Excess Infrastructure Investments (\beta_{kt})
#   ----------------------------------------------------------------------

#   This function ensures that all expansionary supply needs are addressed
#   through infrastructure investments according to the PIP's priority
#   (\beta_k). If under the \beta_k priority distribution an infrastructure type
#   k would be raised above its maximum capacity, the excess capacity is
#   distributed among other infrastructure types. Maintenance, by definition,
#   does not move infrastructure capacity beyond its maximum state, so the
#   priority structure only applies for expansion investments.

function β_k(x,p,t)
    ##Note Supply Needs and Current State
    u_e_need_t = ForwardDiff.value(u_e_need(x,p,t))
    A_t = ForwardDiff.value(A(x,p,t))
    η_t = ForwardDiff.value(x[15])
    μ_t = ForwardDiff.value(μ(x,p,t))
    C_v_s_t = ForwardDiff.value(x[13])
    υ_bar_s_t = ForwardDiff.value(x[5])
    Vbar_s_t = ForwardDiff.value(Vbar_s(x,p,t))
    Vbar_g_t = ForwardDiff.value(Vbar_g(x,p,t))
    μ_s_t = ForwardDiff.value(x[11])
    μ_g_t = ForwardDiff.value(x[12])
    w_max_s_t = ForwardDiff.value(w_max_s(x,p,t))
    w_max_g_t = ForwardDiff.value(w_max_g_proj(x,p,t))
    w_s_t = ForwardDiff.value(x[7])
    w_g_t = ForwardDiff.value(x[8])
    η_proj_t = ForwardDiff.value(η_proj(x,p,t))
    υ_bar_proj_s_t = ForwardDiff.value(υbar_proj_s(x,p,t))
    w_proj_s_t = ForwardDiff.value(w_proj_s(x,p,t))
    w_proj_g_t = ForwardDiff.value(w_proj_g(x,p,t))
    
    #Initialize Vectors for default β_k, default β_dbar_0, counting excess, and noting whether an infrastructure has excess 
    β_k_t = copy(p[19])
    β_dbar_0 = 1-sum(β_k_t)
    excess_k = zeros(length(β_k_t))
    β_over = zeros(length(β_k_t))
    β_max = zeros(length(β_k_t))
    
    ##gather excess β
    for k in 1:length(β_k_t)
        if(β_k_t[k] > 0) #if it is currently being targeted for expansion investments
            #Record what would be the supply increase if using the default β
            u_pot = β_k_t[k]*u_e_need_t
            u_max = u_pot + 1 #in case u_max is not initialized, u_pot will not trigger an excess count
            
            #Calculate Maximum Possible Supply Increases for Each Infrastructure
            if(k==1) #delivery efficiency
                u_max = A_t*(p[11] - η_proj_t)
            elseif(k==2) #storage capacity
                u_max = μ_t*C_v_s_t*(p[12][1] - υ_bar_proj_s_t)
            elseif(k==3) #surface processing
                u_max = (Vbar_s_t + μ_s_t)*(w_max_s_t - w_proj_s_t)
            elseif(k==4) #ground processing
                u_max = (Vbar_g_t + μ_g_t)*(w_max_g_t - w_proj_g_t)
            end
            
            #Note the limit to how much beta can be associated with a certain infrastructure type
            β_max[k] = ifelse(u_max<0,0,u_max/u_e_need_t)
            
            #If Potential Excess, Note it 
            if(u_pot > u_max)
                excess_k[k] = (u_pot - u_max)/u_e_need_t
                if(excess_k[k] > β_k_t[k]) #if u_max neg, excess > β, so all is excess 
                    excess_k[k] = β_k_t[k]
                end
                β_over[k]=1
            end
        end
    end
    
    ############Redistribute excess β among available infrastructure types################
    β_k_new = copy(β_k_t) #vector with final new β values
    excess = sum(excess_k) # Total excess to be re-distributed
    incr = zeros(length(β_k_t)) #Increase in β for each infrastructure type
    
    #Calculate magnitude of available β to be proportioned
    avail = β_dbar_0
    for k in 1:length(β_k_t)
        if(β_over[k]==0)
            avail += β_k_t[k]
        end
    end
    
    #Remove Excess from First Identified Infrastructure Types
    for k in 1:length(β_k_t)
        if(β_over[k]==1)
            β_k_new[k] = β_max[k]
        end 
    end
    
    #Continue making attempted re-distributions until no other infrastructure types are pushed over their limit
    new_overs=ifelse(excess>0,true,false)
    β_over2 = copy(β_over)
    
    counter=0
    while(new_overs)
        counter+=1
        new_overs = false
        for k in 1:length(β_k_t)
            if(β_k_t[k] > 0) #if it is currently being targeted for expansion investments
                if(β_over2[k]==0) #if it has not hit its limit
                    #Use Current Proportional Split to Create a Potential New β
                    β_new_pot_k = β_k_t[k]*(1 + (excess/avail))
                    
                    #If it would push the infrastructure type over the limit, note it (β_over and new_overs), set the new β to β_max, and alter the excess and available
                    if(β_new_pot_k > β_max[k])
                        #Note the new exceedence
                        β_over2[k]=1
                        new_overs=true
                        #Set to max
                        β_k_new[k]=β_max[k]
                        #Lower total excess
                        excess -= (β_k_new[k]-β_k_t[k])
                        #Lower available
                        avail -= β_k_t[k]
                    end
                end
            end
        end
        
        if(counter>10)
            print("Error: too many loops")
            new_overs=false
        end
    end
    
    
    #redistribute to available βs and take excess from over βs
    for k in 1:length(β_k_t)
        if(β_k_t[k] > 0) #if it is currently being targeted for expansion investments
            if(β_over2[k]==0) #if it has not hit its limit
                #Distribute remaining available β according to remaining proportional priority 
                β_k_new[k]= β_k_t[k]*(1 + (excess/avail))
            end
        end
    end
    
    return β_k_new
end;

#   1.3.3.4 Total Needed Investments (in dollars) (J^{need})
#   ----------------------------------------------------------

# :J
# 
#   ^{need}{k,t} = J^{e,need}{k,t} + J^{m,need}_{k,t} $

function J_k_need(x,p,t)
    J_k_m_need_t = J_k_m_need(x,p,t)
    J_k_e_need_t = J_k_e_need(x,p,t)
    
    J_k_need_t = J_k_e_need_t .+ J_k_m_need_t
    
    return J_k_need_t
end;

#   1.3.3.5 Maximum Possible Investment (\bar{J})
#   -----------------------------------------------

# :\
# 
#   bar{J}t = \frac{\taub\left(\bar{R}{t+1} -
#   {C}^{o,proj}{t+1}\right)}{\gammar\left(1 + \taubib \right)} + Rt - C^{o}t -
#   \left(1 + ib \right)\taub \tilde{J}^{b}{t} $
# 
#   where, $ \bar{R}{t+1} = P^{proj,1}t \psir \hat{\pi}t = P^{proj,1}t \psir f_t
#   \bar{\pi} $

function J_bar(x,p,t)
    #Determine Maximum Revenue to get in next year if pursued the highest rate increase possible
    if(x[17] > (1-p[18][2])) #check to see if there is room to raise by ψ_r
        max_change = (1-x[17])/x[17]
    else
        max_change = p[18][2]
    end
    
    R_max = P_proj_1(x,p,t)*(1+max_change)*p[2]*x[17]
    
    #Calculate the Terms in the J_bar equation
    incr = (p[31]*(R_max-C_o_proj(x,p,t)))/(p[17][3]*(1+p[31]*p[32]))
    
    NetOp = R(x,p,t)-C_o(x,p,t)

    J_bar = incr + NetOp - (1+p[32])*p[31]*x[18]
    
    return J_bar
end;

#   1.3.3.6 Final Pursued Investment (After Saturation Checking) (J &
#  u_{2,k})
#   --------------------------

function J_k(x,p,t)
    #Note current investment needs and investment constraints
    J_k = copy(J_k_need(x,p,t))
    J_need_t = sum(J_k)
    J_bar_t = ForwardDiff.value(J_bar(x,p,t))
    
    #If needs surpass constraint, proportionally decrease investment need
    if(J_need_t > J_bar_t)
        J_k = J_k.*(J_bar_t/J_need_t)
    end
    
    return J_k
end;

J(x,p,t) = sum(J_k(x,p,t));

function u_2_k(x,p,t)
    u_2_k = zeros(7) 
    J_k_t = J_k(x,p,t)
    η_t = ForwardDiff.value(x[15])
    
    for k in 1:7
        if(J_k_t[k]>0)
            if(k==2)
                u_2_k[k] = (1/η_t)*(((1/p[22][k+1])*(J_k_t[k]))^(1/p[23][k+2]))
            else
                u_2_k[k] = ((1/p[22][k+1])*(J_k_t[k]))^(1/p[23][k+2])
            end
        end
    end
    
    return u_2_k
end;

#   Direct Revenue-Sourced Investments (J^o) & Bond-Sourced Investments
#  (J^b)
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅

function J_o(x,p,t)
    J_t = J(x,p,t)
    J_o_max = R(x,p,t) - C_o(x,p,t) - C_d(x,p,t)
    
    return min(J_o_max, J_t)
end;

function J_b(x,p,t)
    J_t = J(x,p,t)
    J_o_t = J_o(x,p,t)
    
    return max(J_t - J_o_t, 0)
end;

#   Pursued Maintenance (J^m) and Expansion (J^e) Investments
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅

function J_m(x,p,t)
    return min(J_m_need(x,p,t), J(x,p,t)) 
end;

function J_e(x,p,t)
    return max(J(x,p,t) - J_m(x,p,t), 0)
end;

#   1.3.4 Controller Response: Rate-Setting (u_3)
#   –––––––––––––––––––––––––––––––––––––––––––––––

#   
# u_{3,t} = \frac{\hat{C}^d_t}{\bar{\pi} P_t} e_{3t} Y_{3t}
# $
# 
# *Note that f_change cannot rise above $\psi_r
# 
#   , the maximum rate increase constraint

function u_3(x,p,t)
    f_change = (1/(P(x,p,t)*p[2]))*C_d_need(x,p,t)*e_3(x,p,t)*Y_3(x,p,t)
    
    #Check for maximum rate increase limit
    if(f_change+x[17] > 1)
        f_change = 1 - x[17]
    elseif(f_change/x[17] > p[18][2])
        f_change = p[18][2]*x[17]
    end
    
    return f_change
end;

#   1.3.5 Controller Unifying Equation (u)
#   ––––––––––––––––––––––––––––––––––––––––

function u(x,p,t)
    u_2_k_t = u_2_k(x,p,t)
    u_1_t = ForwardDiff.value(u_1(x,p,t))
    u_3_t = ForwardDiff.value(u_3(x,p,t))
    
    u = copy(u_2_k_t)
    push!(u,u_1_t)
    push!(u,u_3_t)
    
    return u
end;

#   1.4 Equations of Motion
#   =========================

#   15 non-dimensional dynamic variables are stored throughout the run of the
#   dynamical system model. Note: xnew = x_{t+1}
# 
#   1. Population (p=\frac{P}{K}): xnew[1] = p_{t+1} = p_t(1+r(1-p_t) +
#   m_i(1-Y^s_t) - m_o Y^s_t)
# 
#   2. Per Capita Demand (\chi=\frac{d}{\bar{\chi}\mu}): xnew[2] = \chi_{t+1} =
#   \chi_t(1 + \delta_d(1-\frac{\chi_t}{\bar{\chi}_t})) + \frac{u^d_t}{\mu_t}
# 
#   3-4. Reservoir Volume (\upsilon=\frac{V}{\bar{V}}): xnew[3] = \upsilon_{t+1}
#   = \upsilon_t + \frac{ q_t - \frac{O_t(\cdot)}{\mu}}{\beta_tC_{vt}}
# 
#   Phoenix Version: \upsilon_{t+1} = \upsilon_t + \frac{\theta_3q^3_t +
#   \frac{Q^b_t}{\mu_3}- \frac{O_t(\cdot)}{\mu_3}}{\beta_tC_v}
# 
#   5-6. Storage Capacity (\beta = \frac{\bar{V}}{C_v\mu}): xnew[5] =
#   \beta_{t+1} = \beta_t (1-\delta_v) + \frac{u^{\bar{v}}_{t}}{C_{vt}\mu_t}
# 
#   **7-8. Processing Capacity (w = \frac{\bar{O}}{\bar{V}+\mu}): xnew[7] =
#   w_{t+1} = w_t (1-\delta_w) + \frac{u^w_t}{\bar{V}_t+\mu_t}
# 
#   9. Surface Inflow (q^s=\frac{Q^s}{\mu^s}): xnew[9] = q^s_{t+1} =
#   q^s_{ac}(q^s_{t}-1) +C^s_{v,t} \sqrt{1-(q^s_{ac})^2}N(0,1)+1
# 
#   10. Ground Inflow (q^g=\frac{Q^g}{\mu^g}): xnew[10] = q^g_{t+1} =
#   q^g_{ac}(q^g_{t}-1) +C^g_{v,t} \sqrt{1-(q^g_{ac})^2}N(0,1)+1
# 
#   11-12. Mean Inflow(\mu): xnew[11] = \mu_{t+1} = \mu_t + u^\mu_t
# 
#   15. Delivery Efficiency (\eta): xnew[15] = \eta_{t+1} =
#   \eta_t(1-\delta_\eta) + \frac{u^{\eta}_t}{A_t}
# 
#   16. Baseline Demand (\bar{\chi}=\frac{\bar{d}}{\mu}): xnew[16] =
#   \bar{\chi}_{t+1} = \bar{\chi}_t - \frac{u^\bar{d}_t}{\mu_t}
# 
#   17. Rates (f=\frac{\pi}{\bar{\pi}}): xnew[17] = f_{t+1} = f_t + \Delta f_t
# 
#   18-end. Planned Investment (u^*): see if-else statements in the function
#   below

#   Water System Equations of Motion
#   ----------------------------------

#   Implement or Store Long-Term Investment
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅

function ImplementOrStoreLTInvest(x,u,p,t)
    tf = convert(Float64, t)
    u_impl_k = zeros(length(u)-2) #vector for total investments that will be implemented
    u_m_k = zeros(length(u)-2) #vector for maintenance investments that will be implemented
    
    ####Note Maintenance Investments for Relevant Infrastructures
    u_m_k[2] = min(u_m_η(x,p,t), copy(u[2]))
    u_impl_k[2] = u_m_k[2]
    u_m_k[3] = min(u_m_v(x,p,t), copy(u[3]))
    u_impl_k[3] = u_m_k[3]
    u_m_k[4] = min(u_m_w_s(x,p,t), copy(u[4]))
    u_impl_k[4] = u_m_k[4]
    u_m_k[5] = min(u_m_w_g(x,p,t), copy(u[5]))
    u_impl_k[5] = u_m_k[5]
    
    ############Store Un-Implemented Investments
    u_plan_k = zeros(planInvestMemorySize(p))
    
    ############Decide Whether a Potential Investment will be Implemented based on τ_i
    counter = 1 #index to keep track of place
    for k in 1:length(u_impl_k)
        if(p[16][k]==1) ####If τ_i is 1 there is no need for planning or separately accounting for maintenance 
            u_impl_k[k] = copy(u[k])
        else
            u_ready = ForwardDiff.value(x[counter+18]) #Note the investments that will be implemented in t 
            
            u_impl_k[k] = u_impl_k[k] + u_ready #implement the next investment in the queue
            
            if(p[16][k]>2)
                for c in 3:p[16][k] #move up past investments in the queue
                    u_next = ForwardDiff.value(x[counter+19])
                    u_plan_k[counter] = u_next
                    counter += 1
                end
            end
            
            u_plan_k[counter] = u[k] - u_m_k[k] #store non-maintenance investments made in year t
            counter += 1
        end
    end
    
    ############Return Implemented Investments & New Stored Investments
    return [u_impl_k, u_plan_k]
end;

#   Update Infrastructure
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅

function UpdateInfrastructure(x,u_impl,C_v_new,p,t) #H is the implemented H for that year
    #####Setup Vectors & Adjust for ForwardDiff
    I_t = [ForwardDiff.value(x[16])]
    push!(I_t, ForwardDiff.value(x[15]))
    push!(I_t, ForwardDiff.value(x[5]))
    push!(I_t, ForwardDiff.value(x[7]))
    push!(I_t, ForwardDiff.value(x[8]))
    push!(I_t, ForwardDiff.value(x[11]))
    push!(I_t, ForwardDiff.value(x[12]))
    I_new = copy(I_t)
    
    ###Note Other Variables that Have to be Adjusted for ForwardDiff
    C_v_old = [ForwardDiff.value(x[13]) ForwardDiff.value(x[14])]
    υ_bar_g_old = ForwardDiff.value(x[6])
    P_t = ForwardDiff.value(P(x,p,t))
    A_t = ForwardDiff.value(A(x,p,t))
    u_k_t=zeros(7)
    for k in 1:7
        u_k_t[k]=ForwardDiff.value(u_impl[k])
    end
    
    #########Inflow
    #Surface
    if(p[25][1] == 2) #sudden change scenario
        if(p[24] == 1) #PMA Case
            I_new[6] = ifelse(t == p[25][3]-1, (1+p[25][2])*(I_t[6]-p[1][3]) + p[1][3], I_t[6]) + u_k_t[6]
        else
            I_new[6] = ifelse(t == p[25][3]-1, (1+p[25][2])*I_t[6], I_t[6]) + u_k_t[6]
        end
    else #no change scenario
        I_new[6] = I_t[6] + u_k_t[6]
    end  
    
    #Ground - currently no groundwater change scenario
    I_new[7] = I_t[7] + u_k_t[7] 
    
    ########Storage Capacity
    ##Surface
    υ_max_s_2 = I_t[3]*((C_v_old[1]*I_t[6])/(C_v_new[1]*I_new[6]))
    H_vbar_t2 = u_k_t[3]/(C_v_new[1]*I_new[6]) #Convert volumetric investment to storage capacity relevant units
    I_new[3] = υ_max_s_2*(1-p[15][4]) + H_vbar_t2;
    
    ##Ground
    υ_max_g_2 = υ_bar_g_old*((C_v_old[2]*I_t[7])/(C_v_new[2]*I_new[7]))
    υ_bar_g_new = υ_max_g_2
    
    ########Processing Capacity
    #Surface
    w_s_2 = I_t[4]*((I_t[6]*(I_t[3]*C_v_old[1]+1))/(I_new[6]*(I_new[3]*C_v_new[1]+1)))
    H_w_s_t2 = u_k_t[4]/(I_new[6]*(I_new[3]*C_v_new[1]+1)) #Convert volumetric investment to processing capacity relevant units
    I_new[4] = w_s_2*(1-p[15][5]) + H_w_s_t2; #SW pumping capacity
    
    #Ground
    w_g_2 = I_t[5]*((I_t[7]*(υ_bar_g_old*C_v_old[2]+1))/(I_new[7]*(υ_bar_g_new*C_v_new[2]+1)))
    H_w_g_t2 = u_k_t[5]/(I_new[7]*(υ_bar_g_new*C_v_new[2]+1)) #Convert volumetric investment to processing capacity relevant units
    I_new[5] = w_g_2*(1-p[15][6]) + H_w_g_t2; #GW pumping capacity
    
    ########Delivery Efficiency
    η_next=I_t[2]*(1-p[15][3]) +  u_k_t[2]/A_t #Convert volumetric investment to delivery efficiency relevant units and add to existing delivery efficiency
    I_new[2]=ifelse(η_next>p[11],p[11],η_next)
    
    ########Long-Term Demand Management
    χbar_2 = I_t[1]*((I_t[6]+I_t[7])/(I_new[6]+I_new[7]))
    H_dbar_t2 = u_k_t[1]/((I_new[6]+I_new[7])*P_t) #Convert volumetric investment to demand management relevant units
    I_new[1] = χbar_2*(1-p[15][2]) - H_dbar_t2;
    
    ########Return Outputs
    return [I_new, υ_bar_g_new]
end;

#   Water Users
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅

function updateUsers(x,I_new,p,t)
    ######Population Dynamics - Logistic Growth
    p_new = x[1]*(1 + p[6]*(1-x[1])) 
    
    ######Per-Capita Demand Dynamics
    d_ST_t = d_ST(x,p,t)
    χ_ST_t = d_ST_t/μ(x,p,t) #Convert d_ST to χ units
    χ_ST_2 = χ_ST_t*((x[11]+x[12])/(I_new[6]+I_new[7])) #Re-calibrate χ_ST to new streamflow context 
    χbar_2 = x[16]*((x[11]+x[12])/(I_new[6]+I_new[7])) #Re-calibrate χbar to new streamflow context
    
    χ_new = χ_ST_2*(1+p[30]*(1-(χ_ST_2/χbar_2))) #Calculate new PC demand in χ units 
    if(χ_new < p[14]/(I_new[6]+I_new[7])) #Check that new PC demand is above minimum per capita demand
        χ_new = p[14]/(I_new[6]+I_new[7])
    end
    
    ##########Return Outupts
    return [p_new χ_new]
end;

#   Water Balance
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅

function runWaterBalance(x,C_v_new,I_new,υ_bar_g_new,p,t)
    #############Reservoir Volume
    #####Surface Storage 
    #Update
    υ_s_2 = x[3]*((x[5]*x[13]*x[11])/(I_new[3]*C_v_new[1]*I_new[6])) #Re-calibrate υ to new storage capacity and inflow context
    υ_new_s = υ_s_2 + (Q_a_s(x,p,t) - (O_s(x,p,t) + O_f(x,p,t) + Q_b(x,p,t)))/(I_new[3]*C_v_new[1]*I_new[6]);
    
    #Enforce Bounds
    if (υ_new_s>1)
        υ_new_s=1
    elseif(υ_new_s < 0)
        υ_new_s=0
    end
    
    #####Ground Storage
    #Update
    υ_g_2 = x[4]*((x[6]*x[14]*x[12])/(υ_bar_g_new*C_v_new[2]*I_new[7])) #Re-calibrate υ to new storage capacity and inflow context
    υ_new_g = υ_g_2 + (Q_a_g(x,p,t) + Q_b(x,p,t) - O_g(x,p,t))/(υ_bar_g_new*C_v_new[2]*I_new[7])
    
    #Enforce Bounds
    if (υ_new_g>1)
        υ_new_g=1
    elseif(υ_new_g<0)
        υ_new_g=0
    end
    
    ##############Inflows for Next Year
    ####Surface
    if(p[24]==1)#PMA Case, inflow is just the mean
        q_new_s=1
    else
        q_new_s = p[4][1]*(x[9]-1) + x[13]*sqrt(1-p[4][1]*p[4][1])*randn()+1;
    end
    
    ####Ground
    if(p[24]==1)#PMA Case, inflow is just the mean
        q_new_g=1
    else
        q_new_g = p[4][2]*(x[10]-1) + x[14]*sqrt(1-p[4][2]*p[4][2])*randn()+1;
    end
    
    
    #########Return Outputs
    return [υ_new_s υ_new_g q_new_s q_new_g]
end;

#   Full Plant EOM
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅

function WaterSystem(x,u_t,p,t)
    #########Initialize Needed Variables
    u_2_k_t = u_t[1:(end-1)]
    u_1_t = u_t[end-1]
    u_3_t = u_t[end]
    xnew = copy(x)
    
    ####################Implement Investments#######################
    ############Long-Term
    u_update = ImplementOrStoreLTInvest(x,u_t,p,t)
    u_impl_k = u_update[1] 
    u_plan_k = u_update[2] 
    
    ######## Store Planned Investments Made During This Time Period & Their Investment Years############
    ##u_plan
    for c in 1:length(u_plan_k)
        xnew[18+c] = u_plan_k[c]
    end
    
    ####################Update C_v#######################
    ## In default case and PMA case, no update to C_v
    #Surface
    xnew[13] = x[13]
    #Ground
    xnew[14] = x[14]
    
    C_v_new = [ForwardDiff.value(xnew[13]) ForwardDiff.value(xnew[14])]
    
    ####################Update Infrastructure#######################
    Infrast_update = UpdateInfrastructure(x,u_impl_k,C_v_new,p,t)
    xnew[16] = Infrast_update[1][1] #LT Demand Management
    xnew[15] = Infrast_update[1][2] #Delivery Efficiency
    xnew[5] = Infrast_update[1][3] #Surface Storage
    xnew[6] = Infrast_update[2] #Ground Storage
    xnew[7] = Infrast_update[1][4] #Surface Processing Capacity
    xnew[8] = Infrast_update[1][5] #Ground Processing Capacity
    xnew[11] = Infrast_update[1][6] #Surface Mean Inflows
    xnew[12] = Infrast_update[1][7] #Ground Mean Inflows
    
    ####################Water Balance#######################
    WB_new = runWaterBalance(x,C_v_new,Infrast_update[1],Infrast_update[2],p,t)
    
    xnew[3] = WB_new[1] #Surface Reservoir Volume
    xnew[4] = WB_new[2] #Ground Aquifer Volume
    xnew[9] = WB_new[3] #Surface Inflows
    xnew[10] = WB_new[4] #Ground Inflows
    
    ####################Water Users#######################
    Users_update = updateUsers(x,Infrast_update[1],p,t)
    
    xnew[1] = Users_update[1]
    xnew[2] = Users_update[2]
    
    ####################Rates & Bond Investments#######################
    #Rates
    xnew[17] = x[17] + u_3_t 
    #Average Bond Investment 
    xnew[18] = (x[18]*(p[31]-1) + J_b(x,p,t))*(1/p[31])
    
    ##########Return New x
    return xnew
end;

#   All Together EOM
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅

function f(xnew,x,p,t)
    ########### Setup
    ##Convert t to a float
    tf = convert(Float64, t)
    
    ###########Run Controller
    u_t = u(x,p,t)
    
    ###########Run Plant
    x_update = WaterSystem(x,u_t,p,t)
    
    for i in 1:length(x_update)
        xnew[i] = x_update[i]
    end
    
    return
end;

#   2. Model Initialization
#   ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

#   2.1 Default Setup (General City)
#   ==================================

#   Parameter & Initial Condition Definitions
#   –––––––––––––––––––––––––––––––––––––––––––

#   The model is initalized using a setup function Default() that is
#   pre-programmed with all default parameter and initial condition values. The
#   user can change any initial condition or parameter through arguments in
#   Default(). Refer to the table below for the parameters and initial
#   conditions used and what name to use when altering the setup.
# 
#   ***Parameters***
# 
#                                Full Parameter Name Model Variable Name                                                                       Definition                   Units Default Value Allowable Range
#   –––––––––––––––––––––––––––––––––––––––––––––––– ––––––––––––––––––– –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– ––––––––––––––––––––––– ––––––––––––– –––––––––––––––
#                        Max Surface Streamflow Mean             \musmax                                   Max mean surface inflow that the city can seek          Bgal/yr OR AFY        200000      [0,\infty)
#                         Max Ground Streamflow Mean             \mugmax                         Max mean ground inflow (recharge) that the city can seek          Bgal/yr OR AFY             0      [0,\infty)
#                             Max Per-Capita Revenue             \pi_max                                    Max revenue that city can extract per citizen       Dollars/person/yr           400     [0, \infty)
#                 Streamflow Variation Scenario Type        \DeltaCvscen                   Variation Change Type (0 = no change, 1 = gradual, 2 = sudden)                Unitless             0       \{0,1,2\}
#                Surface Streamflow Auto-correlation              \rho_s                             1-year lagged auto-correlation in surface streamflow                Unitless           0.2           [0,1]
#                     Ground Inflow Auto-correlation              \rho_g                                  1-year lagged auto-correlation in ground inflow                Unitless             0           [0,1]
#                                  Carrying Capacity              \kappa                                                     Population carrying capacity                Unitless       1500000     (0, \infty)
#                              Intrinsic Growth Rate                   r                            logistically fit, intrinsic growth rate of population                Unitless           0.1      (0,\infty)
#                            Max Delivery Efficiency            \eta_max                                               Max attainable delivery efficiency                Unitless           2.0          [0, 2]
#                       Max Surface Storage Capacity  \bar{\upsilon}maxs    Max feasible surface storage capacity (multiple of inflow standard deviation)                Unitless            12     [0, \infty)
#                        Max Ground Storage Capacity  \bar{\upsilon}maxg     Max feasible ground storage capacity (multiple of inflow standard deviation)                Unitless            40     [0, \infty)
#                    Max Legal Use of Ground Storage                a_gv                           Max legally allowed use of ground storage (proportion)                Unitless             1           [0,1]
#                     Max Legal Use of Ground Inflow                a_gq                            Max legally allowed use of ground inflow (proportion)                Unitless             1           [0,1]
#                   Max Legal Use of Surface Storage                a_sv                          Max legally allowed use of surface storage (proportion)                Unitless             1           [0,1]
#                    Max Legal Use of Surface Inflow                a_sq                           Max legally allowed use of surface inflow (proportion)                Unitless             1           [0,1]
#                                 Min Per Capita Use               d_min                                                   Min possible per-capita demand (Bgal/yr OR AFY)/person          0.04     [0, \infty)
#                       Background Conservation Rate         \delta_dbar                Annual decay rate of long-term conservation measures (background)                Unitless        0.0116           [0,1]
#                         Surface Storage Decay Rate            \delta_v                                    Annual decay rate of surface storage capacity                Unitless         0.001           [0,1]
#                     Delivery Efficiency Decay Rate         \delta_\eta                                         Annual decay rate of delivery efficiency                Unitless         0.001           [0,1]
#                      Surface Processing Decay Rate            \deltaws                                 Annual decay rate of surface processing capacity                Unitless         0.001           [0,1]
#                       Ground Processing Decay Rate            \deltawg                                  Annual decay rate of ground processing capacity                Unitless         0.001           [0,1]
#                 LT Demand Mgmt Implementation Time              \tau_d                                    Time to implement long-term demand management                     yrs             3      [1,\infty)
#                      Delivery Efficiency Impl Time           \tau_\eta                              Time to implement delievery efficiency improvements                     yrs             4      [1,\infty)
#                         Storage Capacity Impl Time              \tau_v                          Time to implement surface storage capacity improvements                     yrs             5      [1,\infty)
#                       Surface Processing Impl Time              \tauws                       Time to implement surface processing capacity improvements                     yrs             3      [1,\infty)
#                        Ground Processing Impl Time              \tauwg                        Time to implement ground processing capacity improvements                     yrs             3      [1,\infty)
#                Surface Flow Augmentation Impl Time            \tau\mus                         Time to implement surface flow augmentation improvements                     yrs             4      [1,\infty)
#                 Ground Flow Augmentation Impl Time            \tau\mug                          Time to implement ground flow augmentation improvements                     yrs             4      [1,\infty)
#                 Short-Term Goal Supply Sufficiency            \gamma_1                             Goal short-term proportion between supply and demand                Unitless             1      [0,\infty)
#                       Long-Term Goal Supply Buffer            \gamma_2                              Goal long-term proportion between supply and demand                Unitless           1.2      [0,\infty)
#                Minimum Debt Service Coverage Ratio            \gamma_3                          Minimum Debt Service Coverage Ratio Allowed in any Year                unitless             2      [0,\infty)
#                                  Max Rate Increase              \psi_r                                      Max possible proportional increase in rates                Unitless          0.06      [0,\infty)
#                       Delivery Efficiency Priority          \beta_\eta              Proportion of long-term investments directed to delivery efficiency                Unitless           0.2           [0,1]
#                          Storage Capacity Priority             \beta_v         Proportion of long-term investments directed to surface storage capacity                Unitless           0.4           [0,1]
#                        Surface Processing Priority             \betaws      Proportion of long-term investments directed to surface processing capacity                Unitless           0.3           [0,1]
#                         Ground Processing Priority             \betawg       Proportion of long-term investments directed to ground processing capacity                Unitless             0           [0,1]
#                      Surface Augmentation Priority           \beta\mus        Proportion of long-term investments directed to surface flow augmentation                Unitless             0           [0,1]
#                       Ground Augmentation Priority           \beta\mug         Proportion of long-term investments directed to ground flow augmentation                Unitless             0           [0,1]
#                 Short-Term Curtailment Sensitivity           \lambda_1                  Sensitivity (inst. friction component) in short-term investment                Unitless            22     [0, \infty)
#                   Long-Term Investment Sensitivity           \lambda_2                   Sensitivity (inst. friction component) in long-term investment                Unitless            22     [0, \infty)
#                           Rate-Setting Sensitivity           \lambda_3                           Sensitivity (inst. friction component) in rate-setting                Unitless            22     [0, \infty)
#                ST Curtailment Activation Threshold          \epsilon_1         Threshold for action (inst. friction component) in short-term investment                Unitless             0     [0, \infty)
#                 LT Investment Activation Threshold          \epsilon_2          Threshold for action (inst. friction component) in long-term investment                Unitless             0     [0, \infty)
#                  Rate-Setting Activation Threshold          \epsilon_3                  Threshold for action (inst. friction component) in rate-setting                Unitless             0     [0, \infty)
#                Operating Cost Function Coefficient                 g_o                           Coefficient in operating costs function (see equation)   Dollars/(persons*Vol)           2.5     [0, \infty)
#                 LT Dem Mgmt Investment Coefficient              g_dbar                                   Coefficient in LT dem mgmt investment function             AFY/Dollars       4.4 E-7      [0,\infty)
#                     Del Eff Investment Coefficient              g_\eta                           Coefficient in delivery efficiency investment funciton             AFY/Dollars        0.0033      [0,\infty)
#            Storage Capacity Investment Coefficient              g_vbar                              Coefficient in storage capacity investment function              AF/Dollars         0.003     [0, \infty)
#            SW Proc Capacity Investment Coefficient                g_ws                   Coefficient in surface processing capacity investment function             AFY/Dollars       0.00015      [0,\infty)
#            GW Proc Capacity Investment Coefficient                g_wg                    Coefficient in ground processing capacity investment function             AFY/Dollars       0.00015      [0,\infty)
#               SW Inflow Aug Investment Coefficient              g_\mus                   Coefficient in surface inflow augmentation investment function             AFY/Dollars        0.0001      [0,\infty)
#               GW Inflow Aug Investment Coefficient              g_\mug                    Coefficient in ground inflow augmentation investment function             AFY/Dollars        0.0001      [0,\infty)
#             Operating Cost Population Scale Factor                z_op                               Population scale factor in operating cost function                unitless         0.563      [0,\infty)
#                 Operating Cost Demand Scale Factor                z_od                                   Demand scale factor in operating cost function                unitless         0.831      [0,\infty)
#                LT Dem Mgmt Investment Scale Factor              z_dbar                          Investment scale factor for long-term demand management                unitless             1      [0,\infty)
#                    Del Eff Investment Scale Factor              z_\eta                                  Investment scale factor for delivery efficiency                unitless          0.82      [0,\infty)
#           Storage Capacity Investment Scale Factor              z_vbar                                     Investment scale factor for storage capacity                unitless             1      [0,\infty)
#           SW Proc Capacity Investment Scale Factor                z_ws                                  Investment scale factor for processing capacity                unitless             1      [0,\infty)
#           GW Proc Capacity Investment Scale Factor                z_wg                                  Investment scale factor for processing capacity                unitless             1      [0,\infty)
#              SW Inflow Aug Investment Scale Factor              z_\mus                               Investment scale factor for sw inflow augmentation                unitless             1      [0,\infty)
#              GW Inflow Aug Investment Scale Factor              z_\mug                               Investment scale factor for gw inflow augmentation                unitless             1      [0,\infty)
#                                 Case Specification                case                        Indicate which case type to follow (0 = default, 1 = PHX)                Unitless             0         \{0,1\}
#                    Mean Surface Inflow Change Type      \Delta\mustype        Type of change to surface mean inflow (0 = none, 1 = gradual, 2 = sudden)                Unitless             0       \{0,1,2\}
#                 Mean Surface Inflow Percent Change        \Delta\muspc                                            Percent change to surface mean inflow                Unitless             0     [0, \infty)
#                    Mean Sufrace Inflow Change Time         \Delta\must      Time that sudden change occurs or time that gradual change will be complete                     yrs             0     [0, \infty)
#                             Max Surface Processing               wmaxs Max attainable surface processing (proportion of storage capacity + mean inflow)                Unitless             1           [0,1]
#                              Max Ground Processing               wmaxg  Max attainable ground processing (proportion of storage capacity + mean inflow)                Unitless           0.5           [0,1]
#                               Base Groundwater Use            \theta_g                       Proportion of Demand that is usually served by groundwater                Unitless           0.0           [0,1]
#                                   Projection Years              \tau_p                              Years of Projection for Long-Term Investment Signal                   Years             5     [0, \infty)
#   ST Dem Mgmt Investment Effectiveness Coefficient              \alpha                     Coefficient for effectiveness in ST conservation investments                unitless           0.5      [0,\infty)
#                                          Bond Life              \tau_b                                                             Life of issued bonds                   years            15      [0,\infty)
#                                 Bond Interest Rate                 i_b                                                    Interest Rate of Issued Bonds                unitless          0.04           [0,1]
#             Proportion of Rates from Fixed Charges             \beta_p           Proportion of the expected revenue to come from fixed per user charges                unitless           0.5           [0,1]
#              Del Eff Investment Dollars Proportion           \phi_\eta                          Proportion of investment dollars to delivery efficiency                unitless           0.6           [0,1]
#           Stor Capac Investment Dollars Proportion              \phi_v                             Proportion of investment dollars to storage capacity                unitless           0.3           [0,1]
#   Surface Proc Capac Investment Dollars Proportion              \phiws                  Proportion of investment dollars to surface processing capacity                unitless           0.3           [0,1]
#    Ground Proc Capac Investment Dollars Proportion              \phiwg                   Proportion of investment dollars to ground processing capacity                unitless             0           [0,1]
#       Surface Inflow Investment Dollars Proportion            \phi\mus                               Proportion of investment dollars to surface inflow                unitless             0           [0,1]
#        Ground Inflow Investment Dollars Proportion            \phi\mug                                Proportion of investment dollars to ground inflow                unitless             0           [0,1]
# 
#   ***Initial Conditions*** | Full Variable Name | Model Variable Name |
#   Definition | Units |Default Value | Allowable Range | | ––––––– | ––– |
#   ––––– | ––- |––––––- | –––––––- | | Population Fill | p0 | Proportion of
#   population to carrying capacity | Unitless | 0.65 | [0, \infty) | | Actual
#   Per-Capita Demand | \chi0 | per-capita demand, accounting for short-term
#   conservation (proportion of mean inflow) | Unitless | 8E-7 | [0, \infty) | |
#   Surface Storage Fill | \upsilons0 | Fill proportion of surface storage |
#   Unitless | 1 | [0, 1] | | Ground Storage Fill | \upsilong0 | Fill proportion
#   of ground storage | Unitless | 1 | [0, 1] | | Surface Storage Capacity |
#   \upsilonbars0 | Surface storage capacity (multiple of inflow standard
#   deviation) | Unitless | 4 | [0, \infty) | | Ground Storage Capacity |
#   \upsilonbarg0 | Ground storage capacity (multiple of inflow standard
#   deviation) | Unitless | 1 | [0, \infty) | | Surface Processing Capacity |
#   ws0 | Surface processing capacity (proportion of storage capacity & mean
#   inflow) | Unitless | 1 | [0, \infty) | | Ground Processing Capacity | wg0 |
#   Ground processing capacity (proportion of storage capacity & mean inflow) |
#   Unitless | 0 | [0, \infty) | | Surface Inflow | qs0 | Surface inflow
#   (proportion of mean) | Unitless | 1 | [0, \infty) | | Ground Inflow | qg0 |
#   Ground inflow (proportion of mean) | Unitless | 0 | [0, \infty) | | Mean
#   Surface Inflow | \mus0 | Mean surface inflow | AFY or Bgal/yr | 200000 | [0,
#   \infty) | | Mean Ground Inflow | \mug0 | Mean ground inflow | AFY or Bgal/yr
#   | 0 | [0, \infty) | | Surface Inflow Variation | Cvs0 | Surface flow
#   coefficient of variation | Unitless | 0.1 | [0, \infty) | | Ground Inflow
#   Variation | Cvg0 | Ground flow coefficient of variation | Unitless | 0.01 |
#   [0, \infty) | | Delivery Efficiency | \eta0 | Delivery efficiency
#   (proportion of withdrawn delivered) | Unitless | 1 | [0, \infty) | | Base
#   Per-Capita Demand | \chibar0 | Base per-capita demand, independent of ST
#   conservation (proportion of mean inflow) | Unitless | 8E-7 | [0, \infty) | |
#   Per-Capita Revenue | f0 | Proportion of per-capita revenue to max (\pimax) |
#   Unitless | 0.5 | [0,1] | | Average Bond Investment | Jbavg_0 | Average
#   annual bond-sourced investment | $/yr | 69000000 | [0,\infty) |

#   Setup Function (Default)
#   ––––––––––––––––––––––––––

function Default(;μ_s_max = 200000, μ_g_max = 0.001, μ_other_1 = 0, μ_other_2 = 0, π_max = 1000, ΔC_v_scen = 0, ρ_s = 0.2, ρ_g = 0.0, κ = 1500000, r = 0.1, η_max = 1.7, υ_bar_max_s = 12.0, 
        υ_bar_max_g = 40.0, a_gv = 1.0, a_gq = 1.0, a_sv = 1.0, a_sq = 1.0, a_q2 = 1.0, a_q3 = 0, a_q4 = 0, a_q5 = 0, d_min = 0.04, δ_dbar = 0.004, δ_v = 0.07, δ_η = 0.07, δ_w_s = 0.07, 
        δ_w_g = 0.07, τ_d = 1.0, τ_v = 1.0, τ_η = 1.0, τ_w_s = 1.0, τ_w_g = 1.0, τ_μ_s = 1.0, τ_μ_g = 1.0, γ_1 = 1, γ_2 = 1.2, γ_3 = 2, ψ_r = 0.15, β_η = 0.3, β_v = 0.3, β_w_s = 0.3, β_w_g = 0.0, 
        β_μ_s = 0.0, β_μ_g = 0.0, λ_1 = 22.0, λ_2 = 22.0, λ_3 = 22.0, ϵ_1 = 0.0, ϵ_2 = 0.0, ϵ_3 = 0.0, 
        g_o = 1, g_dbar = 2000, g_η = 2000, g_vbar = 2000, g_ws = 2000, g_wg = 2000, g_μs = 2000, 
        g_μg = 2000, z_op = 0.5, z_od = 1.03, z_dbar = 1, z_η = 1, z_vbar = 1, z_ws = 1, z_wg = 1, z_μs = 1, z_μg = 1, case = 0, Δμ_s_type = 0, Δμ_s_pc = 0.0, Δμ_s_t = 0, w_max_s = 1.0, 
        w_max_g = 0.0, θ_g = 0, θ_1 = 0, τ_p = 5, α = 0.5, τ_b = 15, i_b = 0.04, ϕ_η = 0.6, ϕ_v = 0.3, ϕ_w_s = 0.3, ϕ_w_g = 0.0, ϕ_μ_s = 0.0, ϕ_μ_g = 0.0, β_p = 0.5, p_0 = 0.65, χ_0 = 8.0E-7, υ_s_0 = 1.0, 
        υ_g_0 = 1.0, υ_bar_s_0 = 4.0, υ_bar_g_0 = 1.0, w_s_0 = 1.0, w_g_0 = 0.0, q_s_0 = 1.0, q_g_0 = 0.0, μ_s_0 = 200000, μ_g_0 = 0.001, C_v_s_0 = 0.1, C_v_g_0 = 0.01, η_0 = 1.0, χbar_0 = 8.0E-7, 
        f_0 = 0.5, J_b_avg_0 = 69000000)
    
    ##########Parameters############
    #1:μ
    μbar = [μ_s_max μ_g_max μ_other_1 μ_other_2]
    
    #2:π_max
    
    #####Other Streamflow Parameters
    #3: ΔC_v
    ΔC_v = ΔC_v_scen 
    
    #4: ρ
    ρ = [ρ_s, ρ_g]
    
    #24: Select Case
    
    #25: Mean Inflow Change Scneario Definition, Δμ
    Δμ_g_type = Δμ_g_pc = Δμ_g_t = Δμ_s_τ = Δμ_g_τ = 0 #not relevant in this model version
    Δμ = [Δμ_s_type Δμ_s_pc Δμ_s_t Δμ_g_type Δμ_g_pc Δμ_g_t Δμ_s_τ Δμ_g_τ]
        
    #####Population Growth 
        
    #5: Carrying Capacity, κ
        
    #6: Intrinsic Growth rate, r
        
    #7: Immigration Rate, m_i - relevant in other model versions
    m_i = 0
        
    #8: Emigration Rate, m_0 - relevant in other model versions
    m_o = 0
        
    #####Hard Infrastructure Operations
        
    #9: Hedging Policy Coefficient, h - relevant in other model versions 
    h = 0
    
    #10: Minimum Allowed Storage Level, υ_min - just set to 0 because υ is a proportion
    υ_min = [0 0]
        
    #11: Maximum Delivery Efficiency, η
        
    #12: Maximum Storage Capacity, β_max
    υ_bar_max = [υ_bar_max_s υ_bar_max_g]
        
    #13: Maximum Legal Allocation, A
    a = [a_gv a_gq a_sv a_sq a_q2 a_q3 a_q4 a_q5]
        
    #27: Maximum Pumping Capacity, w
    w_max = [w_max_s w_max_g]
        
    #####Water Use
    #14: Minimum Per-Capita Use, d_min
    
    #27: Demand Proportions for Base Use, θ
    θ = [θ_g θ_1]
        
    #####Infrastructure Dynamic Characteristics
    #15: Infrastructure Decay, δ
    δ_d = 0 #not relevant in this version
    δ = [δ_d δ_dbar δ_η δ_v δ_w_s δ_w_g]
        
    #16: Implementation Time, τ_i
    τ_i = [τ_d τ_η τ_v τ_w_s τ_w_g τ_μ_s τ_μ_g]
        
    #####Institutional Context
    #17: Scope Rule, Controller Goals/References
    γ = [γ_1 γ_2 γ_3]
    
    #18: Choice Rules/Constraints
    ψ_s = 0.0 #Not relevant in this model version
    ψ = [ψ_s ψ_r]
        
    #19: Investment Allocation Mental Model (β)
    β = [β_η β_v β_w_s β_w_g β_μ_s β_μ_g]
        
    #20,21: Institutional Friction
    λ = [λ_1 λ_2 λ_3]
    ϵ = [ϵ_1 ϵ_2 ϵ_3]
        
    #####Costs
    #22: Cost Function Coefficients, g
    g = [g_o g_dbar g_η g_vbar g_ws g_wg g_μs g_μg]
        
    #23: Cost Function Scale Factors, z
    z = [z_op z_od z_dbar z_η z_vbar z_ws z_wg z_μs z_μg]
    
    #30: Short-term investment effectiveness coefficient, α
    
    #31: Bond LIfe, τ_b
    
    #32: Bond Interest Rate, i_b
    
    #33: Whether γ_3 error is calculated with absolute value
    γ_3_abs = 0 #not relevant for this model version 
    
    #34: Willingness to Cheat on Groundwater Allowance, β_c
    β_c = 0 #not relevant in this model version
    
    #35: Investment ($) Allocation - just used for g_k calculation
    ϕ = [ϕ_η ϕ_v ϕ_w_s ϕ_w_g ϕ_μ_s ϕ_μ_g]
    
    #36: Proportion of the Expected Per Capita Revenue from Fixed Charges
    
    
    #####Special Options
    #28: Memory of Long-Term Integral Controller
    τ_m = 0 #not relevant in this model version 
        
    #29: Projectioin in Long-Term Controller (τ_p) 
    
    #####Create Parameter Vector
    p = [μbar,π_max,ΔC_v,ρ,κ,r,m_i,m_o,h,υ_min,η_max,υ_bar_max,a,d_min,δ,τ_i,γ,ψ,β,λ,ϵ,g,z,case,Δμ,w_max,θ,τ_m,τ_p,α,τ_b,i_b,γ_3_abs,β_c,ϕ,β_p]
    
    ##########Initial Conditions############
    #####Create Initial Conditions Vector
    x_0 = [p_0, χ_0,υ_s_0,υ_g_0,υ_bar_s_0,υ_bar_g_0,w_s_0,w_g_0,q_s_0,q_g_0,μ_s_0,μ_g_0,C_v_s_0,C_v_g_0,η_0,χbar_0,f_0,J_b_avg_0]
    
    ###Add Memory for Planned Investment Capacity
    for k in 1:length(τ_i)
        if(τ_i[k]>1)
            for i in 1:τ_i[k]-1
                append!(x_0,0.0)
            end
        end
    end
    
    return [p, x_0]
end;

#   General Parameter Auxilary Variables
#   ––––––––––––––––––––––––––––––––––––––

#   Size of Planned Investment Capacity Memory (Number of Variables to add
#  to x based on \tau^i_k)
#   --------------------------

function planInvestMemorySize(p)
    num = 0
    for k in 1:length(p[16])
        if(p[16][k]>1)
            num += p[16][k]-1
        end
    end
    
    return floor(Int,num)
end;

#   2.2 Phoenix Municipal Area (PMA) setups
#   =========================================

#   PMA Parameters & Initial Conditions
#   –––––––––––––––––––––––––––––––––––––

#   ***Parameters***
# 
#                                Full Parameter Name Model Variable Name                                                                  Definition                   Units Default Value Allowable Range
#   –––––––––––––––––––––––––––––––––––––––––––––––– ––––––––––––––––––– ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– ––––––––––––––––––––––– ––––––––––––– –––––––––––––––
#                                    Mean SRP Inflow             \mu_SRP                                     Mean inflow into PMA through SRP canals                     AFY        900000      [0,\infty)
#                                    Mean CAP Inflow             \mu_CAP                                     Mean inflow into PMA through CAP canals                     AFY        650491      [0,\infty)
#                            Mean Groundwater Inflow               \mug0                                            Mean Groundwater Inflow into PMA                     AFY        690602      [0,\infty)
#                             Max Per Capita Revenue             \pi_max                                                      Max Per Capita Revenue                    $/yr          1000      [0,\infty)
#                                  Carrying Capacity              \kappa                                                Population carrying capacity                Unitless       1500000     (0, \infty)
#                              Intrinsic Growth Rate                   r                       logistically fit, intrinsic growth rate of population                Unitless           0.1      (0,\infty)
#                            Max Delivery Efficiency            \eta_max                                          Max attainable delivery efficiency                Unitless           2.0          [0, 2]
#                    Max Legal Use of Ground Storage                a_gv                      Max legally allowed use of ground storage (proportion)                Unitless             1           [0,1]
#                     Max Legal Use of Ground Inflow                a_gq                       Max legally allowed use of ground inflow (proportion)                Unitless             1           [0,1]
#                   Max Legal Use of Surface Storage                a_sv                     Max legally allowed use of surface storage (proportion)                Unitless             0           [0,1]
#                                     SRP Allocation               a_SRP                                  Proportion of SRP inflow allocated to city                Unitless       0.30883           [0,1]
#                                     CAP Allocation               a_CAP                                  Proportion of CAP inflow allocated to city                Unitless       0.41168           [0,1]
#                        CAP Low Priority Allocation             aCAPlow                     Proportion of Low Priority CAP inflow allocated to city                Unitless        0.5324           [0,1]
#                       CAP High Priority Allocation            aCAPhigh                    Proportion of High Priority CAP inflow allocated to city                Unitless        0.3894           [0,1]
#                     SRP NCS & Gatewater Allocation              aSRPNG       Proportion of SRP inflow allocation attributable to NCS and gatewater                Unitless       0.06367           [0,1]
#                                 Min Per Capita Use               d_min                                              Min possible per-capita demand (Bgal/yr OR AFY)/person          0.04     [0, \infty)
#                          Hard Infrastructure Decay              \delta                                             Hard infrastructure decay value                Unitless          0.05           [0,1]
#                       Background Conservation Rate         \delta_dbar           Annual decay rate of long-term conservation measures (background)                Unitless        0.0003           [0,1]
#                 LT Demand Mgmt Implementation Time              \tau_d                               Time to implement long-term demand management                     yrs             1      [1,\infty)
#            Hard Infrastructure Implementation Time              \tau_i                                    Hard infrastructure implementation times                     yrs             3      [1,\infty)
#                       Long-Term Goal Supply Buffer            \gamma_2                         Goal long-term proportion between supply and demand                Unitless           1.2      [0,\infty)
#                Minimum Debt Service Coverage Ratio            \gamma_3                     Minimum Debt Service Coverage Ratio Allowed in any Year                unitless             2      [0,\infty)
#                                  Max Rate Increase              \psi_r                                 Max possible proportional increase in rates                Unitless          0.06      [0,\infty)
#                       Delivery Efficiency Priority          \beta_\eta         Proportion of long-term investments directed to delivery efficiency                Unitless           0.2           [0,1]
#                        Surface Processing Priority             \betaws Proportion of long-term investments directed to surface processing capacity                Unitless             0           [0,1]
#                         Ground Processing Priority             \betawg  Proportion of long-term investments directed to ground processing capacity                Unitless           0.7           [0,1]
#                 Short-Term Curtailment Sensitivity           \lambda_1             Sensitivity (inst. friction component) in short-term investment                Unitless            22     [0, \infty)
#                   Long-Term Investment Sensitivity           \lambda_2              Sensitivity (inst. friction component) in long-term investment                Unitless            22     [0, \infty)
#                           Rate-Setting Sensitivity           \lambda_3                      Sensitivity (inst. friction component) in rate-setting                Unitless            22     [0, \infty)
#                ST Curtailment Activation Threshold          \epsilon_1    Threshold for action (inst. friction component) in short-term investment                Unitless             0     [0, \infty)
#                 LT Investment Activation Threshold          \epsilon_2     Threshold for action (inst. friction component) in long-term investment                Unitless             0     [0, \infty)
#                  Rate-Setting Activation Threshold          \epsilon_3             Threshold for action (inst. friction component) in rate-setting                Unitless             0     [0, \infty)
#                Operating Cost Function Coefficient                 g_o                      Coefficient in operating costs function (see equation)   Dollars/(persons*Vol)        0.1435     [0, \infty)
#                 LT Dem Mgmt Investment Coefficient              g_dbar                              Coefficient in LT dem mgmt investment function             AFY/Dollars          5948      [0,\infty)
#             Operating Cost Population Scale Factor                z_op                          Population scale factor in operating cost function                unitless        0.5581      [0,\infty)
#                 Operating Cost Demand Scale Factor                z_od                              Demand scale factor in operating cost function                unitless        1.0303      [0,\infty)
#                LT Dem Mgmt Investment Scale Factor              z_dbar                     Investment scale factor for long-term demand management                unitless             1      [0,\infty)
#                    Del Eff Investment Scale Factor              z_\eta                             Investment scale factor for delivery efficiency                unitless       1.01266      [0,\infty)
#              Proc Capacity Investment Scale Factor                 z_w                             Investment scale factor for processing capacity                unitless       1.04019      [0,\infty)
#                    Mean Surface Inflow Change Type      \Delta\mustype   Type of change to surface mean inflow (0 = none, 1 = gradual, 2 = sudden)                Unitless             2       \{0,1,2\}
#                 Mean Surface Inflow Percent Change        \Delta\muspc                                       Percent change to surface mean inflow                Unitless        -0.284     [0, \infty)
#                    Mean Sufrace Inflow Change Time         \Delta\must Time that sudden change occurs or time that gradual change will be complete                     yrs            14     [0, \infty)
#                               Base Groundwater Use            \theta_g                  Proportion of Demand that is usually served by groundwater                Unitless         0.024           [0,1]
#                              SRP Demand Proportion            \theta_1                                   Proportion of Demand that is SRP eligible                Unitless           0.5           [0,1]
#                                   Projection Years              \tau_p                         Years of Projection for Long-Term Investment Signal                   Years             5     [0, \infty)
#   ST Dem Mgmt Investment Effectiveness Coefficient              \alpha                Coefficient for effectiveness in ST conservation investments                unitless           0.5      [0,\infty)
#                                          Bond Life              \tau_b                                                        Life of issued bonds                   years            15      [0,\infty)
#                                 Bond Interest Rate                 i_b                                               Interest Rate of Issued Bonds                unitless          0.04           [0,1]
#             Proportion of Rates from Fixed Charges             \beta_p      Proportion of the expected revenue to come from fixed per user charges                unitless           0.5           [0,1]
#              Del Eff Investment Dollars Proportion           \phi_\eta                     Proportion of investment dollars to delivery efficiency                unitless        0.6605           [0,1]
#   Surface Proc Capac Investment Dollars Proportion              \phiws             Proportion of investment dollars to surface processing capacity                unitless        0.2964           [0,1]
#    Ground Proc Capac Investment Dollars Proportion              \phiwg              Proportion of investment dollars to ground processing capacity                unitless        0.0331           [0,1]
#                           Initial SRP Availability               ASRP0                                       Initial Volume of Available SRP Water                     AFY     200275.18      [0,\infty)
# 
#   ***Initial Conditions***
# 
#           Full Variable Name Model Variable Name                                                                            Definition    Units Default Value Allowable Range
#   –––––––––––––––––––––––––– ––––––––––––––––––– ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– –––––––– ––––––––––––– –––––––––––––––
#              Population Fill                 p_0                                         Proportion of population to carrying capacity Unitless       0.86466     [0, \infty)
#     Actual Per-Capita Demand              \chi_0 per-capita demand, accounting for short-term conservation (proportion of mean inflow) Unitless    8.96274E-8     [0, \infty)
#          Ground Storage Fill          \upsilong0                                                     Fill proportion of ground storage Unitless       0.53257          [0, 1]
#   Ground Processing Capacity                 wg0             Ground processing capacity (proportion of storage capacity & mean inflow) Unitless      0.005315     [0, \infty)
#     Surface Inflow Variation               Cvs_0                                                 Surface flow coefficient of variation Unitless         0.001     [0, \infty)
#      Ground Inflow Variation               Cvg_0                                                  Ground flow coefficient of variation Unitless         0.001     [0, \infty)
#          Delivery Efficiency              \eta_0                               Delivery efficiency (proportion of withdrawn delivered) Unitless        0.9742     [0, \infty)
#       Base Per-Capita Demand           \chibar_0    Base per-capita demand, independent of ST conservation (proportion of mean inflow) Unitless    8.96274E-8     [0, \infty)
#           Per-Capita Revenue                 f_0                                     Proportion of per-capita revenue to max (\pi_max) Unitless       0.23836           [0,1]
#      Average Bond Investment             Jbavg_0                                                Average annual bond-sourced investment     $/yr      69694375      [0,\infty)

#   Determine Investment Cost Function Coefficients from Other Parameters
#  (PMA-Specific Function)
#   ––––––––––––––––––––––––––

function set_g(μ_s_0,μ_g_0,μ_tot,μ_SRP,μ_CAP,a_SRP,a_CAP,a_SRP_NG,a_gq,a_gv,θ_1,V_g_0, A_l_g_0,V_bar_g_0,
        κ_endog,κ,p_0,χbar_0,χ_0,f_0,π_bar,g_o,z_op,z_od,J_b_avg_0,τ_b,i_b,η_0,w_g_0,w_s_0,δ_η,δ_w_s,
        δ_w_g,z_η, z_w,ϕ_η,ϕ_w,A_SRP_0)
    
    #Determine Maintenance Investment Need ($)
    a_q_0 = (a_SRP*μ_SRP + a_CAP*(μ_s_0-μ_SRP) + a_gq*μ_g_0)/μ_tot
    
    #Set Initial Conditions
    P_0 = p_0*κ
    D_0 = χ_0*P_0*μ_tot
    R_0 = f_0*π_bar*P_0
    C_o_0 = g_o*(P_0^z_op)*(D_0^z_od)
    C_d_0 = (1+τ_b*i_b)*J_b_avg_0
    J_o_0 = C_d_0
    
    #Determine Capacity Need (H) Based on δ
    A_l_s_SRP = A_SRP_0
    A_l_s_CAP = μ_CAP*a_CAP
    A_l_s_0 = A_l_s_SRP + A_l_s_CAP
    A_w_s_0 = w_s_0*μ_s_0
    A_w_g_0 = w_g_0*(μ_g_0+V_bar_g_0)
    
    A_s=min(A_l_s_0,A_w_s_0)
    A_g=min(A_l_g_0,A_w_g_0)
    A_0=A_s+A_g
    
    u_m_η = (η_0)*A_0*δ_η
    u_m_ws = w_s_0*μ_s_0*δ_w_s 
    u_m_wg = w_g_0*(μ_g_0+V_bar_g_0)*δ_w_g
    
    #Determine Available Maintenance Investment (J_m)
    J_m = J_b_avg_0 + J_o_0
    
    J_m_η = ϕ_η*J_m
    J_m_w = ϕ_w*J_m
    
    #Calculate g as a function of H_m, U_m, and z
    g_η = J_m_η/((η_0*u_m_η)^z_η)
    g_w = J_m_w/((u_m_ws^z_w)+(u_m_wg^z_w)) #treating the g's for processing capacity as the same across surface/ground
    
    return[g_η g_w g_w]
end;

#   City of Phoenix Setup
#   –––––––––––––––––––––––

function Phoenix(;μ_SRP = 900000,  μ_CAP = 448663, μ_g_0 = 690602, π_max = 1000, κ = 1686528, r = 0.0878, η_max = 1.5, a_gv = 1, a_gq = 0.05357, a_sv = 0.0, a_SRP = 0.30883, 
        a_CAP = 0.41168, a_CAP_low=0.5324,a_CAP_high=0.3894,a_SRP_NG=0.06367, d_min = 0.072, δ=0.05, δ_dbar = 0.0003, τ=3.0, γ_2 = 1.2, γ_3 = 2, ψ_r = 0.15, β_η = 0.2, 
        β_w_s = 0, β_w_g = 0.7, λ_1 = 110, λ_2 = 22, λ_3 = 22, ϵ_1 = 0.0, ϵ_2 = 0.0, ϵ_3 = 0.0, g_o = 0.1435, g_dbar = 5948,  z_op = 0.5581, z_od = 1.0303, z_dbar = 1, z_η = 1.01266,  
        z_w=1.04019, Δμ_s_type = 2, Δμ_s_pc = -0.284, Δμ_s_t = 14, θ_g = 0.024, θ_1 = 0.5, τ_p = 5, α = 0.5, τ_b = 15, i_b=0.04, ϕ_η = 0.6605, ϕ_w_s = 0.2964, ϕ_w_g = 0.0331, 
        β_p = 0.6322, QC=0, A_SRP_0=200275.18, p_0 = 0.86466, χ_0 = 8.96274E-8, υ_g_0 = 0.53257, υ_bar_s_0 = 7.415E-9, w_s_0 = nothing, w_g_0 = 0.005315, C_v_s_0 = 0.001,C_v_g_0=0.001, 
        η_0 = 0.9742, χbar_0 = 8.96274E-8, f_0 = 0.23836, J_b_avg_0 = 69694375)
    
    ####Adjust J_b_avg_0 to τ_b and ι if necessary
    C_d_0 = J_b_avg_0*(1+15*0.04) #with default settings
    J_b_avg_0 = C_d_0/(1+τ_b*i_b) 
    
    ####μ        
    μ_s_0 = μ_SRP + μ_CAP
    μ_tot = μ_s_0 + μ_g_0
    
    ####Set υ_g_0, β_g_0 and β_g_max
    V_bar_g_0 = (200*a_gq*μ_g_0) #double the 100 years of safe yield alloance
    υ_bar_g_0 = V_bar_g_0/(C_v_g_0*μ_g_0)
    υ_bar_g_max = υ_bar_g_0
    V_g_0 = υ_g_0*V_bar_g_0
    
    ####W_max
    P_0 = p_0*κ
    D_0 = χ_0*P_0*μ_tot
    #CAP Avail
    A_l_s_CAP = μ_CAP*a_CAP
    #SRP Avail
    A_l_s_SRP = A_SRP_0
    #Total Surface
    A_l_s_0 = A_l_s_SRP + A_l_s_CAP
    #Max Surface Processing
    w_max_s = A_l_s_0/(υ_bar_s_0*C_v_s_0*μ_s_0 + μ_s_0)
    #Legal GW Avail
    A_l_g_0 = a_gq*μ_g_0 + a_gv*max(0,V_g_0-(100*a_gq*μ_g_0))
    #Max Ground Processing
    w_max_g = (A_l_g_0*1.6)/(V_bar_g_0+μ_g_0) #1.6 accounts for the desire to have surplus capacity for peak intrannual demands
    
    ####Initial Surface Processing
    if(w_s_0 == nothing)
        w_s_0 = w_max_s
    end
    
    ####Implementation Times, τ
    τ_η =  τ_w_s = τ_w_g = τ
    
    ####Hard Infrastructure Decay Rates, δ
    δ_η = δ_w_s = δ_w_g = δ
    
    ####Investment Cost Function Normalizing Coefficients, g
    z_ws = z_w
    z_wg = z_w
    if(QC==1) #Queen Creek is calibrated to 2016 data to account for H20 aquisition
        p_0_QC=74842/κ #2016
        χbar_0_QC = 16344.01/(μ_tot*p_0_QC*κ) #2016
        χ_0_QC=χbar_0_QC
        f_0_QC=23734654/(p_0_QC*κ*π_max) #2016
        J_b_avg_0_QC=3492917.5 #2016
        
        g=set_g(μ_s_0,μ_g_0,μ_tot,μ_SRP,μ_CAP,a_SRP,a_CAP,a_SRP_NG,a_gq,a_gv,θ_1,V_g_0, A_l_g_0, V_bar_g_0, 0,κ,p_0_QC,χbar_0_QC,
            χ_0_QC,f_0_QC,π_max,g_o,z_op,z_od,J_b_avg_0_QC,τ_b,i_b,η_0,w_g_0,w_s_0,δ_η,δ_w_s, δ_w_g,z_η, z_w, ϕ_η, (ϕ_w_s+ϕ_w_g),A_SRP_0)
    else
        g=set_g(μ_s_0,μ_g_0,μ_tot,μ_SRP,μ_CAP,a_SRP,a_CAP,a_SRP_NG,a_gq,a_gv,θ_1,V_g_0, A_l_g_0, V_bar_g_0, 0,κ,p_0,χbar_0,χ_0,
            f_0,π_max,g_o,z_op,z_od,J_b_avg_0,τ_b,i_b,η_0,w_g_0,w_s_0,δ_η,δ_w_s, δ_w_g,z_η, z_w, ϕ_η, (ϕ_w_s+ϕ_w_g),A_SRP_0)
    end
    g_η=g[1]
    g_ws=g[2]
    g_wg=g[3]
    
    
    return Default(μ_s_max = μ_s_0, μ_g_max = μ_g_0, μ_other_1=μ_SRP, μ_other_2=μ_CAP, π_max = π_max, κ=κ, r = r, η_max = η_max, υ_bar_max_s = 6.45E-9, υ_bar_max_g = υ_bar_g_max, a_gv = a_gv, 
        a_gq = a_gq, a_sv = a_sv, a_sq = a_SRP, a_q2 = a_CAP, a_q3=a_CAP_low, a_q4=a_CAP_high, a_q5 = a_SRP_NG,  d_min = d_min, δ_dbar = δ_dbar, δ_v = 0.0, δ_η = δ_η, δ_w_s = δ_w_s, 
        δ_w_g = δ_w_g, τ_d = 1.0, τ_η = τ_η, τ_w_s = τ_w_s, τ_w_g = τ_w_g, γ_2 = γ_2, γ_3 = γ_3, ψ_r = ψ_r, β_η = β_η, β_v = 0, β_w_s = β_w_s, β_w_g = β_w_g, β_μ_s = 0.0, β_μ_g = 0.0, 
        λ_1 = λ_1, λ_2 = λ_2, λ_3 = λ_3, ϵ_1 = ϵ_1, ϵ_2 = ϵ_2, ϵ_3 = ϵ_3, g_o = g_o, g_dbar = g_dbar, g_η = g_η, g_vbar = 0, g_μs = 0, g_μg = 0, g_ws = g_ws, g_wg = g_wg, z_op = z_op, 
        z_od = z_od, z_dbar = z_dbar, z_η = z_η, z_vbar = 1, z_ws = z_ws, z_wg = z_wg, z_μs = 1, z_μg = 1, case = 1, Δμ_s_type = Δμ_s_type, Δμ_s_pc = Δμ_s_pc, Δμ_s_t = Δμ_s_t, w_max_s = w_max_s, 
        w_max_g = w_max_g, θ_g = θ_g, θ_1 = θ_1, τ_p = τ_p, α=α, τ_b = τ_b, i_b=i_b, ϕ_η = ϕ_η, ϕ_v = 0, ϕ_w_s = ϕ_w_s, ϕ_w_g = ϕ_w_g, ϕ_μ_s = 0, ϕ_μ_g = 0, β_p=β_p, p_0 = p_0, χ_0 = χ_0, 
        υ_s_0 = 0.0, υ_g_0 = υ_g_0, υ_bar_s_0 = υ_bar_s_0, υ_bar_g_0 = υ_bar_g_0, w_s_0 = w_s_0, w_g_0 = w_g_0, q_s_0 = 1.0, q_g_0 = 1.0, μ_s_0 = μ_s_0, μ_g_0 = μ_g_0, C_v_s_0 = C_v_s_0, 
        C_v_g_0 = C_v_g_0, η_0 = η_0, χbar_0 = χbar_0, f_0 = f_0, J_b_avg_0=J_b_avg_0) 
end;

#   City of Scottsdale Setup
#   ––––––––––––––––––––––––––

function Scottsdale(;π_max = 1000, κ = 242300, r = 0.143, η_max = 1.5, a_gv = 1, a_gq = 0.03114, a_sv = 0.0, a_SRP = 0.03067, a_CAP = 0.18076, a_CAP_low=0.0472, a_CAP_high = 0.2055,
        a_SRP_NG=0, d_min = 0.074, δ=0.05, δ_dbar = 0.0003, τ=3.0, γ_2 = 1.2, γ_3 = 2, ψ_r = 0.15, β_η = 0.2, β_w_s = 0, β_w_g = 0.7, λ_1 = 110, λ_2 = 22, λ_3 = 22, ϵ_1 = 0.0, ϵ_2 = 0.0, 
        ϵ_3 = 0.0, g_o = 0.4316, z_op = 0.5581, z_od = 1.0303, z_η = 1.01266, z_w = 1.04019, Δμ_s_type = 2, Δμ_s_pc = -0.284, Δμ_s_t = 14, θ_g = 0.075, θ_1 = 0.17, τ_p = 5, α = 0.5, 
        τ_b = 15, i_b=0.04, ϕ_η = 0.6605, ϕ_w_s = 0.2409, ϕ_w_g = 0.0886, β_p = 0.3586, A_SRP_0=13632.80, p_0 = 0.89948, χ_0 = 1.86072E-7, υ_g_0 = 0.52182, w_s_0 = nothing, 
        w_g_0 = 0.006977, η_0 = 1.0562, χbar_0 = 1.86072E-7, f_0 = 0.39811, J_b_avg_0 = 12500000)

    
    
    return Phoenix(π_max = π_max, κ=κ, r = r, η_max = η_max, a_gv = a_gv, a_gq = a_gq, a_sv = a_sv, a_SRP = a_SRP, a_CAP = a_CAP, a_CAP_low=a_CAP_low,a_CAP_high=a_CAP_high, a_SRP_NG=a_SRP_NG, 
        d_min = d_min, δ=δ, δ_dbar = δ_dbar, τ=τ, γ_2 = γ_2, γ_3 = γ_3, ψ_r = ψ_r, β_η = β_η, β_w_s = β_w_s, β_w_g = β_w_g, λ_1 = λ_1, λ_2 = λ_2, λ_3 = λ_3, ϵ_1 = ϵ_1, ϵ_2 = ϵ_2, ϵ_3 = ϵ_3,
        g_o = g_o, z_op = z_op, z_od = z_od, z_η = z_η, z_w = z_w, Δμ_s_type = Δμ_s_type, Δμ_s_pc = Δμ_s_pc, Δμ_s_t = Δμ_s_t, θ_g = θ_g, θ_1 = θ_1, τ_p = τ_p, α=α, τ_b = τ_b, i_b=i_b, 
        ϕ_η = ϕ_η, ϕ_w_s = ϕ_w_s, ϕ_w_g = ϕ_w_g, β_p=β_p, QC=0, A_SRP_0=A_SRP_0, p_0 = p_0, χ_0 = χ_0, υ_g_0 = υ_g_0, w_s_0 = w_s_0, w_g_0 = w_g_0, η_0 = η_0, χbar_0 = χbar_0, f_0 = f_0,
        J_b_avg_0=J_b_avg_0)
end;

#   Town of Queen Creek Setup
#   –––––––––––––––––––––––––––

function QueenCreek(;π_max = 1000, κ = 101553, r = 0.24, η_max = 1.5, a_gv = 1, a_gq = 0.02135, a_sv = 0.0, a_SRP = 0.0, a_CAP = 0.01038, a_CAP_low=0.0594, a_CAP_high = 0.0013,a_SRP_NG=0,
        d_min = 0.063, δ=0.05, δ_dbar = 0.0003, τ=3.0, γ_2 = 1.2, γ_3 = 2, ψ_r = 0.15, β_η=0.04, β_w_s=0.04, β_w_g=0.9, λ_1 = 110, λ_2 = 22, λ_3 = 22, ϵ_1 = 0.0, ϵ_2 = 0.0, ϵ_3 = 0.0, 
        g_o = 0.8346, z_op = 0.5581, z_od = 1.0303, z_η = 1.01266, z_w = 1.04019, Δμ_s_type = 2, Δμ_s_pc = -0.284, Δμ_s_t = 14, θ_g = 0.075, θ_1 = 0.0, τ_p = 5, α = 0.5, τ_b = 15, 
        i_b=0.04, ϕ_η = 0.6605, ϕ_w_s = 0.0088, ϕ_w_g = 0.3207, β_p = 0.5801, A_SRP_0=0, p_0 = 0.31705, χ_0 = 1.55365E-7, υ_g_0 = 0.59817, w_s_0 = 0.0005, w_g_0 =0.006737, η_0 = 0.9339, 
        χbar_0 = 1.55365E-7, f_0 = 0.24106, J_b_avg_0 = 1322425.63)
    
    return Phoenix(π_max = π_max, κ=κ, r = r, η_max = η_max, a_gv = a_gv, a_gq = a_gq, a_sv = a_sv, a_SRP = a_SRP, a_CAP = a_CAP, a_CAP_low=a_CAP_low,a_CAP_high=a_CAP_high, a_SRP_NG=a_SRP_NG, 
        d_min = d_min, δ=δ, δ_dbar = δ_dbar, τ=τ, γ_2 = γ_2, γ_3 = γ_3, ψ_r = ψ_r, β_η = β_η, β_w_s = β_w_s, β_w_g = β_w_g, λ_1 = λ_1, λ_2 = λ_2, λ_3 = λ_3, ϵ_1 = ϵ_1, ϵ_2 = ϵ_2, ϵ_3 = ϵ_3,
        g_o = g_o, z_op = z_op, z_od = z_od, z_η = z_η, z_w = z_w, Δμ_s_type = Δμ_s_type, Δμ_s_pc = Δμ_s_pc, Δμ_s_t = Δμ_s_t, θ_g = θ_g, θ_1 = θ_1, τ_p = τ_p, α=α, τ_b = τ_b, i_b=i_b, 
        ϕ_η = ϕ_η, ϕ_w_s = ϕ_w_s, ϕ_w_g = ϕ_w_g, β_p=β_p, QC=1, A_SRP_0=A_SRP_0, p_0 = p_0, χ_0 = χ_0, υ_g_0 = υ_g_0, w_s_0 = w_s_0, w_g_0 = w_g_0, η_0 = η_0, χbar_0 = χbar_0, f_0 = f_0,
        J_b_avg_0=J_b_avg_0)
end;

#   3.Running the Model
#   ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

#   3.1 Specify Input
#   ===================

 function create_UWIIM(setup::Any=Default(); f=f)
    p = copy(setup[1])
    x_0 = Float64.(setup[2])
    
    model = DiscreteDynamicalSystem(f, x_0, p)
    
    return model
end;

#   3.2 Create Output
#   ===================

#   Every time the model is run (run_UWIIM()) an array of four components is
#   returned to visualize the trajectory
# 
#     1. Dataframe of state and auxilary variables over time
# 
#     2. Non-dimensional Time Series Plots (vector of 13 plots) (see
#        sub-section b)
# 
#     3. Dimensional Time Series Plots (vector of 13 plots) (see
#        sub-section b)

#   a.) Generate Dataframe of Trajectory Variables
#   ––––––––––––––––––––––––––––––––––––––––––––––––

#   Reported Output Variables
#   ---------------------------

#                                         Full Variable Name Model Variable Name                                                                                        Definition      Units
#   –––––––––––––––––––––––––––––––––––––––––––––––––––––––– ––––––––––––––––––– ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– ––––––––––
#                                                  Time Step                   t                                                                           Time Step (initial = 0)      years
#                                                       Year                year                                                                                 Year of time step      years
#                               Population (Non-Dimensional)                   p                                           Service population as a proportion of carrying capacity   Unitless
#                                                 Population                   P                                                                                Service population    persons
#                        Per Capita Demand (Non-Dimensional)                   χ                                            Per capita demand as a proportion of total mean inflow   Unitless
#               Baseline Per Capita Demand (Non-Dimensional)                χbar                                   Baseline per capita demand as a proportion of total mean inflow   Unitless
#                                              Annual Demand                   D                                               Total annual system demand at the start of the year        AFY
#                                           Projected Demand              D_proj                                                                    Projected annual system demand        AFY
#                                            Baseline Demand               D_bar                                                                       Baseline demand in the year        AFY
#                                    Demand Post-Curtailment                D_ST                                                       System demand after curtailment in the year        AFY
#                         Per Capita Demand Post-Curtailment                d_ST                                                   Per capita demand after curtailment in the year AFY/person
#             Surface Water Reservoir Fill (Non-Dimensional)                 υ_s                                                      Proportional fill of surface water reservoir   Unitless
#                 Groundwater Aquifer Fill (Non-Dimensional)                 υ_g                                                          Proportional fill of groundwater aquifer   Unitless
#                 Surface Storage Capacity (Non-Dimensional)               υbars                       Surface water storage capacity as a proportion of inflow standard deviation   Unitless
#                  Ground Storage Capacity (Non-Dimensional)               υbarg                         Groundwater storage capacity as a proportion of inflow standard deviation   Unitless
#                                       Stored Surface Water                 V_s                               Volume of water stored in surface water reservoirs at start of year         AF
#                                         Stored Groundwater                 V_g                                   Volume of water stored in groundwater aquifers at start of year         AF
#                                   Surface Storage Capacity              Vbar_s                                                                    Surface water storage capacity         AF
#                                    Ground Storage Capacity              Vbar_g                                                                     Ground water storage capacity         AF
#       Projected Groundwater Aquifer Fill (Non-Dimensional)              υprojg                                                Projected proportional fill of groundwater aquifer   Unitless
#              Surface Processing Capacity (Non-Dimensional)                 w_s                  Surface water processing capacity as a proportion of storage capacity and inflow   Unitless
#               Ground Processing Capacity (Non-Dimensional)                 w_g                    Groundwater processing capacity as a proportion of storage capacity and inflow   Unitless
#                        Maximum Surface Processing Capacity               wmaxs          Maximum surface water processing capacity as a proportion of storage capacity and inflow   Unitless
#                         Maximum Ground Processing Capacity               wmaxg            Maximum groundwater processing capacity as a proportion of storage capacity and inflow   Unitless
#                      Projected Surface Processing Capacity              wsproj        Projected surface water processing capacity as a proportion of storage capacity and inflow   Unitless
#                           Surface Inflow (Non-Dimensional)                 q_s                                               Surface water inflow as a proportion of mean inflow   Unitless
#                            Ground Inflow (Non-Dimensional)                 q_g                                                 Groundwater inflow as a proportion of mean inflow   Unitless
#                                   Inflow (Non-Dimensional)                   q                                                       Total inflow as a proportion of mean inflow   Unitless
#                                             Surface Inflow                 Q_s                                                                  Surface water inflow in the year        AFY
#                                              Ground Inflow                 Q_g                                                                    Groundwater inflow in the year        AFY
#                                               Total Inflow                   Q                                                                          Total inflow in the year        AFY
#                                           Available Inflow                 Q_a                                                         Total inlow available for the city to use        AFY
#                         Available Inflow (Non-Dimensional)                 q_a                    Total inlow available for the city to use as a proportion of total mean inflow   Unitless
#                                   Available Surface Inflow                 Qas                                          Total surface water inflow available for the city to use        AFY
#                                    Available Ground Inflow                 Qag                                            Total groundwater inflow available for the city to use        AFY
#                                               Banked Water                 Q_b                                    Surface inflow water stored in groundwater storage in the year        AFY
#                                        Mean Surface Inflow                 μ_s                                                                  Mean annual surface water inflow        AFY
#                                         Mean Ground Inflow                 μ_g                                                                    Mean annual groundwater inflow        AFY
#                                          Mean Total Inflow                   μ                                                                          Mean annual total inflow        AFY
#                    Surface Inflow Coefficient of Variation                 Cvs                                          Coefficient of variation for surface water annual inflow   Unitless
#                     Ground Inflow Coefficient of Variation                 Cvg                                            Coefficient of variation for groundwater annual inflow   Unitless
#                                        Delivery Efficiency                   η                                            Delivery efficiency as a proportion of available water   Unitless
#                              Projected Delivery Efficiency              η_proj                                  Projected delivery efficiency as a proportion of available water   Unitless
#                                             Total Outflows                   O                                      All water releases and uses from storage sources in the year        AFY
#                           Total Outflows (Non-Dimensional)                   o All water releases and uses from storage sources in the year as a proportion of total mean inflow   Unitless
#                                    Demand-Related Outflows                 O_d                                                Water uses from storage to meet demand in the year        AFY
#                  Demand-Related Outflows (Non-Dimensional)                 O_d           Water uses from storage to meet demand in the year as a proportion of total mean inflow   Unitless
#                                    Demand-Related Outflows                 O_d                                                Water uses from storage to meet demand in the year        AFY
#                  Demand-Related Outflows (Non-Dimensional)                 O_d           Water uses from storage to meet demand in the year as a proportion of total mean inflow   Unitless
#                                           Surface Outflows                 O_s                                                   All surface water releases and uses in the year        AFY
#                                           Surface Outflows                 o_s                    All surface water releases and uses in the year as a proportion of mean inflow   Unitless
#                                            Ground Outflows                 O_g                                                     All groundwater releases and uses in the year        AFY
#                                            Ground Outflows                 o_g                      All groundwater releases and uses in the year as a proportion of mean inflow   Unitless
#                                             Flood Releases                 O_f                                                   Water released from storage to prevent overflow        AFY
#                                             Flood Releases                 o_f              Water released from storage to prevent overflow as a proportion of total mean inflow   Unitless
#                  Projected Legally Available Surface Water            Als_proj                                                         Projected legally available surface water        AFY
#                    Projected Legally Available Groundwater            Alg_proj                                                           Projected legally available groundwater        AFY
#                          Projected Available Surface water              Asproj                                                                 Projected available surface water        AFY
#                            Projected Available Groundwater              Agproj                                                                   Projected available groundwater        AFY
#                                    Available Surface Water                 A_s                                                  Total available surface water to use in the year        AFY
#                                      Available Groundwater                 A_g                                                    Total available groundwater to use in the year        AFY
#                                            Available Water                   A                                                          Total available water to use in the year        AFY
#                            Legally Available Surface Water                 Als                                                 Total legally available surface water in the year        AFY
#                              Legally Available Groundwater                 Alg                                                   Total legally available groundwater in the year        AFY
#                        Technically Available Surface Water                 Ats                                             Total technically available surface water in the year        AFY
#                          Technically Available GroundWater                 Atg                                               Total technically available groundwater in the year        AFY
#                                                     Supply                   S                                                                          Total Supply in the year        AFY
#                                           Projected Supply              S_proj                                                                            Projected total supply        AFY
#                                Shortage Before Curtailment               ω_pre                                        Ratio of supply deficit to total demand before curtailment   Unitless
#                                 Shortage After Curtailment              ω_post                                         Ratio of supply deficit to total demand after curtailment   Unitless
#                                              Safety Factor                  SF                                                             Ratio of Supply to Demand in the year   Unitless
#                                    Projected Safety Factor             SF_proj                                                               Projected ratio of supply to demand   Unitless
#                                Debt Service Coverage Ratio                DSCR                                      Ratio of net revenue to debt service requirement in the year   Unitless
#                                           Short-Term Error                 e_1                                   Error in the short-term curtailment action situation/controller   Unitless
#                                            Long-Term Error                 e_2                                     Error in the long-term investment action situation/controller   Unitless
#                                         Rate-Setting Error                 e_3                                             Error in the rate-setting action situation/controller   Unitless
#                                     Short-Term Curtailment                 u_1      Short-term curtailment pursued in the year as a proportion of total mean inflow (terms of χ)   Unitless
#                                       Short-Term Attention                 Y_1                               Attention in the short-term curtailment action situation/controller   Unitless
#                                        Long-Term Attention                 Y_2                                 Attention in the long-term investment action situation/controller   Unitless
#                                     Rate-Setting Attention                 Y_3                                         Attention in the rate-setting action situation/controller   Unitless
#                                         Per Capita Revenue                   f                    Annual per capita revenue as a proportion of maximum annual per capita revenue   Unitless
#                                                    Revenue                   R                                                               Total revenue generated in the year          $
#                                            Operating Costs                 C_o                                                              Operating costs required in the year       $/yr
#                                               Debt Service                 C_d                                                                 Debt service required in the year       $/yr
#                                             Investment ($)                   J                                                    Total investment in infrastructure in the year       $/yr
#                          Needed Maintenance Investment ($)              Jmneed                                          Needed infrastructure maintenance investment in the year       $/yr
#                                 Maintenance Investment ($)                 J_m                                                    Maintenance investment implemented in the year       $/yr
#                                Expansionary Investment ($)                 J_e                               Expansionary (increase capacity) investment implemented in the year       $/yr
#                                Average Bond Investment ($)               Jbavg                                                    Average bond-sourced investment over all years       $/yr
#                                        Bond Investment ($)                 J_b                                                               Bond-sourced investment in the year       $/yr
#                              Direct Revenue Investment ($)                 J_o                                                      Direct investment of net revenue in the year       $/yr
#                             Maximum Investment Allowed ($)               J_bar                                                Maximum investment that can be pursued in the year       $/yr
#                              Expansionary Investment (AFY)              ueneed                               Expansionary (increase capacity) investment implemented in the year        AFY
#                    Delivery Efficiency Investment Priority                 β_η                    Proportion of expansionary investment in the year going to delivery efficiency   Unitless
#            Surface Processing Capacity Investment Priority                 βws            Proportion of expansionary investment in the year going to surface processing capacity   Unitless
#             Ground Processing Capacity Investment Priority                 βwg             Proportion of expansionary investment in the year going to ground processing capacity   Unitless
#                        Implemented Demand Investment (AFY)           uimpldbar                                              Investment in baseline demand management in the year        AFY
#           Implemented Delivery Efficiency Investment (AFY)              uimplη                                                     Investment in delivery efficiency in the year        AFY
#      Implemented Surface Storage Capacity Investment (AFY)           uimplvbar                                          Investment in surface water storage capacity in the year        AFY
#   Implemented Surface Processing Capacity Investment (AFY)            uimplw_s                                       Investment in surface water processing capacity in the year        AFY
#    Implemented Ground Processing Capacity Investment (AFY)            uimplw_g                                         Investment in groundwater processing capacity in the year        AFY
#                Implemented Surface Inflow Investment (AFY)            uimplμ_s                                                    Investment in surface water inflow in the year        AFY
#                 Implemented Ground Inflow Investment (AFY)            uimplμ_g                                                      Investment in groundwater inflow in the year        AFY
# 
#   ***PMA City-Unique Variables*** | Full Variable Name | Model Variable Name |
#   Definition | Units | | ––––––– | ––– | ––––– | ––- | | SRP Use | O1 | Water
#   used from SRP | AFY | | SRP Use (Non-Dimensional) | o1 | Water used from SRP
#   as a proportion of mean inflow | Unitless | | CAP Use | O2 | Water used from
#   CAP | AFY | | CAP Use (Non-Dimensional) | o2 | Water used from CAP as a
#   proportion of mean inflow | AFY | | SRP Inflow | Q1 | Inflow from SRP into
#   the PMA in the year | AFY | | CAP Inflow | Q2 | Inflow from CAP into the PMA
#   in the year | AFY |

#   Output Variable Dataframe Function
#   ------------------------------------

function createVarsDF(UWIIM,p,num_t,year_0)
   ##Generate Trajectory
    tr = trajectory(UWIIM, num_t)
    
    ##Collect Into DataFrame
    t = collect(0:num_t)
    col_copy = copy(collect(columns(tr)))

    p_ex, χ, υ_s, υ_g, υ_bar_s, υ_bar_g, w_s, w_g, q_s_ex, q_g_ex, μ_s, μ_g, C_v_s, C_v_g, η, χbar, f, J_b_avg = columns(tr)
    
    μ_ex = zeros(length(tr.data))
    V_s_ex = zeros(length(tr.data))
    V_g_ex = zeros(length(tr.data))
    Vbar_s_ex = zeros(length(tr.data))
    Vbar_g_ex = zeros(length(tr.data))
    υ_proj_g_ex = zeros(length(tr.data))
    year_ex = zeros(length(tr.data))
    D_ex = zeros(length(tr.data))
    d_ST_ex = zeros(length(tr.data))
    D_ST_ex = zeros(length(tr.data))
    D_proj_ex = zeros(length(tr.data))
    D_bar_ex = zeros(length(tr.data))
    O_ex = zeros(length(tr.data))
    o_ex = zeros(length(tr.data))
    O_d_ex = zeros(length(tr.data))
    o_d_ex = zeros(length(tr.data))
    O_s_ex = zeros(length(tr.data))
    o_s_ex = zeros(length(tr.data))
    O_g_ex = zeros(length(tr.data))
    o_g_ex = zeros(length(tr.data))
    O_f_ex = zeros(length(tr.data))
    o_f_ex = zeros(length(tr.data))
    q_ex = zeros(length(tr.data))
    Q_ex = zeros(length(tr.data))
    Q_a_ex = zeros(length(tr.data))
    q_a_ex = zeros(length(tr.data))
    Q_a_s_ex = zeros(length(tr.data))
    Q_a_g_ex = zeros(length(tr.data))
    Q_s_ex = zeros(length(tr.data))
    Q_g_ex = zeros(length(tr.data))
    P_ex = zeros(length(tr.data))
    κbar_ex = zeros(length(tr.data))
    S_ex = zeros(length(tr.data))
    S_proj_ex = zeros(length(tr.data))
    SF_ex = zeros(length(tr.data))
    SF_proj_ex = zeros(length(tr.data))
    DSCR_ex=zeros(length(tr.data)) 
    e_1_ex = zeros(length(tr.data))
    e_2_ex = zeros(length(tr.data))
    e_3_ex = zeros(length(tr.data))
    u_1_ex = zeros(length(tr.data))
    R_ex = zeros(length(tr.data))
    Y_1_ex = zeros(length(tr.data))
    Y_2_ex = zeros(length(tr.data))
    Y_3_ex = zeros(length(tr.data))
    C_o_ex = zeros(length(tr.data))
    C_d_ex = zeros(length(tr.data))
    J_ex = zeros(length(tr.data)) 
    J_bar_ex = zeros(length(tr.data))
    A_s_ex = zeros(length(tr.data)) 
    A_g_ex = zeros(length(tr.data)) 
    A_l_s_ex = zeros(length(tr.data)) 
    A_l_g_ex = zeros(length(tr.data)) 
    A_w_s_ex = zeros(length(tr.data)) 
    A_w_g_ex = zeros(length(tr.data)) 
    A_ex = zeros(length(tr.data))
    w_max_s_ex = zeros(length(tr.data))
    w_max_g_ex = zeros(length(tr.data))
    J_m_need_ex = zeros(length(tr.data))
    u_e_need_ex = zeros(length(tr.data))
    u_impl_dbar_ex = zeros(length(tr.data))
    u_impl_η_ex = zeros(length(tr.data))
    u_impl_vbar_ex = zeros(length(tr.data))
    u_impl_w_s_ex = zeros(length(tr.data))
    u_impl_w_g_ex = zeros(length(tr.data))
    u_impl_μ_s_ex = zeros(length(tr.data))
    u_impl_μ_g_ex = zeros(length(tr.data))
    J_m_ex = zeros(length(tr.data))
    J_e_ex = zeros(length(tr.data))
    J_o_ex = zeros(length(tr.data))
    J_b_ex = zeros(length(tr.data))
    A_l_g_proj_ex=zeros(length(tr.data))
    A_g_proj_ex=zeros(length(tr.data))
    A_s_proj_ex=zeros(length(tr.data))
    A_l_s_proj_ex=zeros(length(tr.data))
    η_proj_ex = zeros(length(tr.data))
    w_s_proj_ex = zeros(length(tr.data))
    β_η_ex = zeros(length(tr.data))
    β_w_s_ex = zeros(length(tr.data))
    β_w_g_ex = zeros(length(tr.data))
    Q_b_ex = zeros(length(tr.data))
    
    
    #Note Phoenix Specific Variables
    if(p[24]==1)
        O_ex_1 = zeros(length(tr.data))
        O_ex_2 = zeros(length(tr.data))
        o_ex_1 = zeros(length(tr.data))
        o_ex_2 = zeros(length(tr.data))
        Q_1_ex = zeros(length(tr.data))
        Q_2_ex = zeros(length(tr.data))
    end
    
    for i in 1:length(tr.data)
        μ_ex[i] = μ(tr.data[i],p,t[i])
        year_ex[i] = year_0 + (i-1)
        V_s_ex[i] = V_s(tr.data[i],p,t[i])
        V_g_ex[i] = V_g(tr.data[i],p,t[i])
        Vbar_s_ex[i] = Vbar_s(tr.data[i],p,t[i])
        Vbar_g_ex[i] = Vbar_g(tr.data[i],p,t[i])
        υ_proj_g_ex[i] = υ_and_η_proj_g(tr.data[i],p,t[i])[1]
        D_ex[i] = D(tr.data[i],p,t[i])
        d_ST_ex[i] = d_ST(tr.data[i],p,t[i])
        D_ST_ex[i] = D_ST(tr.data[i],p,t[i])
        D_proj_ex[i] = D_proj(tr.data[i],p,t[i])
        D_bar_ex[i] = D_bar(tr.data[i],p,t[i])
        O_ex[i] = O(tr.data[i],p,t[i]) 
        o_ex[i] = O_ex[i]/μ_ex[i]
        O_d_ex[i] = O_d(tr.data[i],p,t[i])
        o_d_ex[i] = O_d_ex[i]/μ_ex[i]
        O_s_ex[i] = O_s(tr.data[i],p,t[i]) 
        o_s_ex[i] = O_s_ex[i]/μ_ex[i]
        O_g_ex[i] = O_g(tr.data[i],p,t[i]) 
        o_g_ex[i] = O_g_ex[i]/μ_ex[i]
        O_f_ex[i] = O_f(tr.data[i],p,t[i]) 
        o_f_ex[i] = O_f_ex[i]/μ_ex[i]
        q_ex[i] = q(tr.data[i],p,t[i]) 
        Q_ex[i] = Q(tr.data[i],p,t[i]) 
        Q_a_ex[i] = Q_a(tr.data[i],p,t[i]) 
        q_a_ex[i] = q_a(tr.data[i],p,t[i]) 
        Q_a_s_ex[i] = Q_a_s(tr.data[i],p,t[i]) 
        Q_a_g_ex[i] = Q_a_g(tr.data[i],p,t[i])  
        Q_s_ex[i] = q_s_ex[i]*μ_s[i]
        Q_g_ex[i] = q_g_ex[i]*μ_g[i]
        Q_b_ex[i]=Q_b(tr.data[i],p,t[i])
        P_ex[i] = P(tr.data[i],p,t[i])
        S_ex[i] = S(tr.data[i],p,t[i])
        S_proj_ex[i] = S_proj(tr.data[i],p,t[i])
        SF_ex[i] = M_1(tr.data[i],p,t[i])
        SF_proj_ex[i] = M_2(tr.data[i],p,t[i])
        DSCR_ex[i] = M_3(tr.data[i],p,t[i])
        e_1_ex[i] = e_1(tr.data[i],p,t[i])
        e_2_ex[i] = e_2(tr.data[i],p,t[i])
        e_3_ex[i] = e_3(tr.data[i],p,t[i])
        u_1_ex[i] = u_1(tr.data[i],p,t[i])
        R_ex[i] = R(tr.data[i],p,t[i])
        Y_1_ex[i] = Y_1(tr.data[i],p,t[i])
        Y_2_ex[i] = Y_2(tr.data[i],p,t[i])
        Y_3_ex[i] = Y_3(tr.data[i],p,t[i])
        C_o_ex[i] = C_o(tr.data[i],p,t[i])
        C_d_ex[i] = C_d(tr.data[i],p,t[i])
        J_ex[i] = J(tr.data[i],p,t[i])
        J_bar_ex[i] = J_bar(tr.data[i],p,t[i])    
        A_s_ex[i] = A_s(tr.data[i],p,t[i])
        A_g_ex[i] = A_g(tr.data[i],p,t[i])
        A_s_proj_ex[i]=A_proj_s(tr.data[i],p,t[i])
        A_g_proj_ex[i]=A_proj_g(tr.data[i],p,t[i])
        A_l_s_ex[i] = A_l_s(tr.data[i],p,t[i])
        A_l_s_proj_ex[i] = A_proj_l_s(tr.data[i],p,t[i])
        A_l_g_proj_ex[i] = A_proj_l_g(tr.data[i],p,t[i])
        A_l_g_ex[i] = A_l_g(tr.data[i],p,t[i])
        A_w_s_ex[i] = A_w_s(tr.data[i],p,t[i])
        A_w_g_ex[i] = A_w_g(tr.data[i],p,t[i])
        A_ex[i] = A(tr.data[i],p,t[i])
        w_max_s_ex[i] = w_max_s(tr.data[i],p,t[i])
        w_max_g_ex[i] = w_max_g(tr.data[i],p,t[i])
        J_m_need_ex[i] = J_m_need(tr.data[i],p,t[i]) 
        J_m_ex[i] = J_m(tr.data[i],p,t[i])
        J_e_ex[i] = J_e(tr.data[i],p,t[i])
        J_b_ex[i] = J_b(tr.data[i],p,t[i])
        J_o_ex[i] = J_o(tr.data[i],p,t[i])
        u_e_need_ex[i] = u_e_need(tr.data[i],p,t[i])
        η_proj_ex[i] = η_proj(tr.data[i],p,t[i])
        w_s_proj_ex[i] = w_proj_s(tr.data[i],p,t[i])
        β_k_t = β_k(tr.data[i],p,t[i])
        β_η_ex[i] = β_k_t[1]
        β_w_s_ex[i] = β_k_t[3]
        β_w_g_ex[i] = β_k_t[4]
        
        #Record Implemented Expansionary Investments
        u_t = u(tr.data[i],p,t[i])
        u_impl_k_ex = ImplementOrStoreLTInvest(tr.data[i],u_t,p,t[i])[1]
        u_impl_dbar_ex[i] = u_impl_k_ex[1]
        u_impl_η_ex[i] = u_impl_k_ex[2]
        u_impl_vbar_ex[i] = u_impl_k_ex[3]
        u_impl_w_s_ex[i] = u_impl_k_ex[4]
        u_impl_w_g_ex[i] = u_impl_k_ex[5]
        u_impl_μ_s_ex[i] = u_impl_k_ex[6]
        u_impl_μ_g_ex[i] = u_impl_k_ex[7]
        
        #PMA Specific Auxiliary Variables of Interest
        if(p[24]==1) 
            O_ex_1[i] = O_1(tr.data[i],p,t[i])
            O_ex_2[i] = O_2(tr.data[i],p,t[i])
            o_ex_1[i] = O_ex_1[i]/(μ_s[i]+μ_g[i])
            o_ex_2[i] = O_ex_2[i]/(μ_s[i]+μ_g[i])
            Q_1_ex[i] = p[1][3]*p[13][4]
            
            Q_CAP = (Q_s(tr.data[i],p,t[i]) - p[1][3])
            CAP_short = p[1][4] - Q_CAP
            NIA_avail = max(0,70022-CAP_short)
            high_avail = Q_CAP - NIA_avail
            Q_2_ex[i] = p[13][6]*NIA_avail + p[13][7]*high_avail
        end
    end
    
    ##Shortage Calculations
    ω_pre = ifelse.(D_bar_ex.> S_ex,(D_bar_ex-S_ex)./D_bar_ex,0)
    ω_post = ifelse.(D_ST_ex.> S_ex,(D_ST_ex-S_ex)./D_ST_ex,0)
    
    ##Aggregate Variables into Single Dataframe
    if(p[24]==1)
        vars = DataFrame(t=t,year=year_ex,p=p_ex, χ=χ, υ_s=υ_s, υ_g=υ_g, υ_bar_s=υ_bar_s, υ_bar_g=υ_bar_g, w_s=w_s, w_g=w_g, q_s=q_s_ex, q_g=q_g_ex, q=q_ex, Q_s=Q_s_ex, Q_g=Q_g_ex, Q=Q_ex, μ_s=μ_s, μ_g=μ_g, μ=μ_ex,
            C_v_s=C_v_s, C_v_g=C_v_g, η=η, χbar=χbar, f=f,  V_s = V_s_ex, V_g = V_g_ex, Vbar_s = Vbar_s_ex, Vbar_g = Vbar_g_ex, D=D_ex, D_proj=D_proj_ex,υ_proj_g=υ_proj_g_ex, O=O_ex, o=o_ex, O_d = O_d_ex, o_d = o_d_ex, 
            O_s=O_s_ex, o_s=o_s_ex, O_g=O_g_ex, o_g=o_g_ex, O_f=O_f_ex, o_f=o_f_ex, O_1=O_ex_1, o_1=o_ex_1, O_2=O_ex_2, o_2=o_ex_2, P=P_ex, S=S_ex, S_proj=S_proj_ex,Q_a = Q_a_ex, q_a = q_a_ex, 
            Q_a_s = Q_a_s_ex, Q_a_g = Q_a_g_ex, SF=SF_ex, SF_proj=SF_proj_ex,DSCR=DSCR_ex, e_1=e_1_ex, e_2 = e_2_ex, e_3 = e_3_ex, u_1=u_1_ex, R=R_ex, Y_1=Y_1_ex, Y_2=Y_2_ex, Y_3 = Y_3_ex, C_o = C_o_ex, 
            C_d = C_d_ex, J=J_ex, ω_pre=ω_pre,ω_post=ω_post, Q_1=Q_1_ex, Q_2=Q_2_ex, A_s=A_s_ex,A_g=A_g_ex,A=A_ex, A_l_s=A_l_s_ex, A_l_g=A_l_g_ex, A_w_s=A_w_s_ex,
            A_w_g=A_w_g_ex, w_max_s = w_max_s_ex, w_max_g = w_max_g_ex, J_m_need=J_m_need_ex, u_impl_dbar =u_impl_dbar_ex,u_impl_η =u_impl_η_ex,u_impl_vbar =u_impl_vbar_ex, 
            u_impl_w_s =u_impl_w_s_ex,u_impl_w_g =u_impl_w_g_ex, u_impl_μ_s =u_impl_μ_s_ex,u_impl_μ_g =u_impl_μ_g_ex, J_m=J_m_ex, J_e=J_e_ex, J_b_avg=J_b_avg,J_bar=J_bar_ex,J_b=J_b_ex,J_o=J_o_ex, 
            D_bar=D_bar_ex,D_ST = D_ST_ex, d_ST = d_ST_ex,u_e_need=u_e_need_ex,A_l_g_proj=A_l_g_proj_ex,A_g_proj=A_g_proj_ex, A_s_proj=A_s_proj_ex,η_proj=η_proj_ex,
            A_l_s_proj=A_l_s_proj_ex, β_η=β_η_ex,β_w_s=β_w_s_ex,β_w_g=β_w_g_ex,Q_b=Q_b_ex,w_s_proj=w_s_proj_ex)
    else 
        vars = DataFrame(t=t,year=year_ex,p=p_ex, χ=χ, υ_s=υ_s, υ_g=υ_g, υ_bar_s=υ_bar_s, υ_bar_g=υ_bar_g, w_s=w_s, w_g=w_g, q_s=q_s_ex, q_g=q_g_ex, q=q_ex, Q_s=Q_s_ex, Q_g=Q_g_ex, Q=Q_ex, μ_s=μ_s, μ_g=μ_g, μ=μ_ex,
            C_v_s=C_v_s, C_v_g=C_v_g, η=η, χbar=χbar, f=f, V_s = V_s_ex, V_g = V_g_ex, Vbar_s = Vbar_s_ex,υ_proj_g=υ_proj_g_ex,Vbar_g = Vbar_g_ex, D=D_ex, D_proj=D_proj_ex, O=O_ex, o=o_ex, O_d = O_d_ex, o_d = o_d_ex, 
            O_s=O_s_ex, o_s=o_s_ex, O_g=O_g_ex, o_g=o_g_ex, O_f=O_f_ex, o_f=o_f_ex, P=P_ex, S=S_ex, S_proj=S_proj_ex, Q_a = Q_a_ex, q_a = q_a_ex, Q_a_s = Q_a_s_ex, Q_a_g = Q_a_g_ex, SF=SF_ex, SF_proj=SF_proj_ex, 
            DSCR=DSCR_ex, e_1=e_1_ex, e_2 = e_2_ex, e_3 = e_3_ex, u_1=u_1_ex, R=R_ex, Y_1=Y_1_ex, Y_2=Y_2_ex, Y_3 = Y_3_ex, C_o = C_o_ex, C_d = C_d_ex, J=J_ex, ω_pre=ω_pre,ω_post=ω_post,
            A_s=A_s_ex,A_g=A_g_ex,A=A_ex,A_l_s=A_l_s_ex, A_l_g=A_l_g_ex, A_w_s=A_w_s_ex,A_w_g=A_w_g_ex,w_max_s = w_max_s_ex, w_max_g = w_max_g_ex,J_m_need=J_m_need_ex,
            u_impl_dbar =u_impl_dbar_ex,u_impl_η =u_impl_η_ex,u_impl_vbar =u_impl_vbar_ex, u_impl_w_s =u_impl_w_s_ex,u_impl_w_g =u_impl_w_g_ex,
            u_impl_μ_s =u_impl_μ_s_ex,u_impl_μ_g =u_impl_μ_g_ex, J_m=J_m_ex, J_e=J_e_ex,J_b_avg=J_b_avg,J_bar=J_bar_ex,J_b=J_b_ex,J_o=J_o_ex, D_bar=D_bar_ex,D_ST = D_ST_ex, d_ST = d_ST_ex,
            u_e_need=u_e_need_ex,A_l_g_proj=A_l_g_proj_ex,A_g_proj=A_g_proj_ex,A_s_proj=A_s_proj_ex,η_proj=η_proj_ex,A_l_s_proj=A_l_s_proj_ex, β_η=β_η_ex,β_w_s=β_w_s_ex,β_w_g=β_w_g_ex)
    end

    return vars
end;

#   b.) Generate Time Series Plots
#   ––––––––––––––––––––––––––––––––

#   Description of Plots
#   ----------------------

#   All Time Series Plotting Functions return a vector of 13 plots:
# 
#   i. Shortage Plot (Supply, Demand, Shortage)
# 
#   ii. Flows Plot (Inflows & Use)
# 
#   iii. Storage Plot (Fill Volume & Storage Capacity)
# 
#   iv. Error
# 
#   v. Attention
# 
#   vi. Financial Flows (Revenue, Costs, Investments)
# 
#   vii. Per-Capita Revenue
# 
#   viii. Delivery Efficiency
# 
#   ix. Per-Capita Demand (Base & Actual)
# 
#   x. Processing Capacity
# 
#   xi. Population
# 
#   xii. Mean Inflow (including augmentation)
# 
#   xiii. All plots

#   Non-dimensional Time Series
#   -----------------------------

function timeSeriesPlot_nondim(vars,p,x_0)
   
    ##Plot Test Plots
    #Shortage, Demand, Supply, Population
    plt_short=plot(vars.t, [vars.S./vars.μ vars.χ.*vars.P vars.χbar.*vars.P vars.ω_pre vars.ω_post], labels = ["Supply" "Demand" "Demand_base" "Shortage_PreCons" "Shortage_w/Cons"], xlabel = "Year", ylabel = "Flow/μ", 
        linecolor = [:green :blue :indigo :pink :red], title = "Shortage",legend=:outerright)
    
    #Flows
    if(p[24]==1)
        plt_flows = plot(vars.t, [vars.q_a vars.o_d], labels=["In_avail" "Use"], ylims = (0,vars.q_a[1]*1.5), xlabel = "Year", ylabel="Flow/μ", title = "Inflows & Use", 
        legend = :outerright, linecolor = [:darkgreen :darkgreen], linestyle=[:solid :dot])
    else
        plt_flows = plot(vars.t, [vars.q_a vars.Q_a_s./vars.μ p[13][4].*(1 .+ vars.C_v_s) p[13][4].*(1 .- vars.C_v_s) vars.Q_a_g./vars.μ vars.o_d vars.o_s vars.o_g vars.o_f], 
            labels= ["In_avail" "In_s" "+σ_s" "-σ_s" "In_g" "Use_all" "Use_s" "Use_g" "Use_f"], ylims = (0,vars.q_a[1]*1.5), xlabel = "Year", ylabel="Flow/μ", title = "Inflows & Use", legend = :outerright, 
            linecolor = [:black :green :green :green :turquoise :black :green :turquoise :brown], linestyle=[:solid :solid :dash :dash :solid :dot :dot :dot :dot])
    end
    
    #Storage Volume & Capacity
    plt_stor = plot(vars.t, [vars.υ_s.*vars.υ_bar_s vars.υ_bar_s], labels = ["Stor_Vol" "Stor_Capac"], ylims = (0,p[12][1].*1.5), xlabel = "Year", ylabel="Vol/(μ*C_v)", title = "Reservoir Storage", 
        legend = :outerright, linecolor = [:grey :black])
    
    #Error
    plt_e = plot(vars.t,[vars.e_1 vars.e_2 vars.e_3], labels = ["Short-Term" "Invest" "Rates"], xlabel = "Year", ylabel = "Error", title = "Error", 
        legend = :outerright, ylims = (-3,3), linecolor = [:blue :green :brown])
    
    #Attention
    plt_ρ = plot(vars.t,[vars.Y_1 vars.Y_2 vars.Y_3],labels = ["Short-Term" "Invest" "Rates"], xlabel="Year", ylabel = "Attention", 
        title = "Attention", legend=:outerright, linecolor = [:blue :green :brown])

    #Financial Flows
    plt_fin = plot(vars.t,[vars.R.*0.000001 vars.C_o.*0.000001 vars.C_d.*0.000001 vars.J.*0.000001./p[16][2]], labels = ["Rev" "Op Costs" "Debt Serv" "Inv Long"], 
        xlabel="Year", ylabel = "Dollars (M)", title = "Financial Flows", legend=:outerright, ylims=(0,(p[2]*p[5].*0.000001)), 
        linecolor = [:brown :pink :blue :green])
    
    #Rate-Setting
    plt_f = plot(vars.t,vars.f, labels = "f", xlabel="Time(years)", ylabel = "Per-Capita Rates/Max Per-Capita Rates", title = "Rates - Max Rates Ratio", legend=:outerright, 
        ylims=(0,1.2), linecolor = :brown)
    
    #Delivery Efficiency
    plt_η = plot(vars.t, [vars.η vars.u_impl_η./vars.A], xlabel = "Year",labels = ["η" "Invest_η"],ylabel="Del. Eff (%Outflow)", title = "Delivery Efficiency", 
        linecolor = :purple, linestyle = [:solid :dot], legend=:outerright, ylims = (0,2.2))
    hline!([p[11]], labels="η_max", linestyle=:dash, linecolor=:purple)
    
    #Per-Capita Demand
    plt_d = plot(vars.t, [vars.χ vars.χbar (vars.u_1./vars.μ).*(-1) (vars.u_impl_dbar./vars.μ).*(-1) p[14]./vars.μ], xlabel = "Year", ylabel="Flows/cap/μ",
        labels=["Dem_PC" "BaseDem_PC" "Conserv_ST" "Conserv_LT" "Dem_PC_min"], title="Per-Capita Demand", ylims = (0,1.1*vars.χbar[1]), 
        legend=:outerright, linecolor = [:blue :indigo :blue :indigo :indigo], linestyle = [:solid :solid :dot :dot :dash])
    
    #Pumping Capacity
    plt_w = plot(vars.t, [vars.w_s vars.w_g vars.u_impl_w_s./(vars.Vbar_s.+vars.Q_a_s) vars.u_impl_w_g./(vars.Vbar_s.+vars.Q_a_s) vars.w_max_s vars.w_max_g], xlabel = "Year",
        labels = ["w_s" "Invest_w_s" "w_g" "Invest_w_g" "w_s_max" "w_g_max"], ylabel="Surface Pump Capac (%Storage)", title = "Processing Capacity", 
        linecolor = [:green :turquoise :green :turquoise :green :turquoise], linestyle = [:solid :solid :dot :dot :dash :dash], legend=:outerright, ylims = (0,1.2))
    
    #Mean Inflow
    plt_μ = plot(vars.t, [vars.μ_s./vars.μ[1] vars.μ_g./vars.μ[1]], xlabel = "Year", labels = ["μ_s" "μ_g"], ylabel = "Mean Inflow/Total Initial Mean Inflow (μ_0)", title = "Mean River Inflow", 
        linecolor = [:green :turquoise], linestyle = :solid, legend = :outerright, ylims = (0,1.5))
    hline!([p[1][1]/vars.μ[1]], labels="μ_s_max", linestyle=:dash, linecolor=:green)
    hline!([p[1][2]/vars.μ[1]], labels="μ_g_max", linestyle=:dash, linecolor=:turquoise)
    
    #Population
    plt_P = plot(vars.t, vars.P.*0.001, ylabel="1000 Persons", xlabel = "Year", title = "Population", labels="Pop", ylims=(0,p[5]*1.2*0.001), 
        linecolor = :orange, linestyle = [:solid :dash :dot], legend = :outerright)
    
    #Aggregate Plots
    plt = plot(plt_short, plt_flows, plt_stor, plt_e, plt_ρ, plt_fin, plt_f, plt_η, plt_d, plt_w, plt_P, plt_μ, size = (1800,1500), layout = (4,3))
    
    #Save in Plot List
    plt_list = [plt_short plt_flows plt_stor plt_e plt_ρ plt_fin plt_f plt_η plt_d plt_w plt_P plt_μ plt]
    
    return plt_list
end;

#   Dimensional Time Series
#   -------------------------

function timeSeriesPlot_dim(vars,p,x_0,units)
    ##Note Key Times; If Phoenix scenario, note CAP shock time. If general, not time that demand reaches μ
    if(p[24]==1)
        t_CAP = p[25][3]
        year_CAP = vars.year[1] + t_CAP
    end
    
    ##Plot Test Plots
    #Shortage, Demand, Supply, Population
    if(units == "AF")
        plt_short=plot(vars.year, [vars.S.*0.001 vars.D_ST.*0.001 vars.χbar.*vars.P.*vars.μ.*0.001 vars.ω_pre.*vars.D_bar.*0.001 vars.ω_post.*vars.D_ST.*0.001], 
            labels = ["Supply" "Demand" "Demand_base" "Shortage_preCons" "Shortage_w/Cons"], xlabel = "Year", ylabel = "KAFY", 
            linecolor = [:green :blue :indigo :pink :red], title = "Supply & Demand",legend=:outerright, legendtitle="Supply & Demand")
    else
        plt_short=plot(vars.year, [vars.S.*0.000000001 vars.D_ST.*0.000000001 vars.χbar.*vars.P.*vars.μ.*0.000000001 vars.ω_pre.*vars.D_bar.*0.000000001 vars.ω_post.*vars.D_ST.*0.000000001], 
            labels = ["Supply" "Demand" "Demand_base" "Shortage_preCons" "Shortage_w/Cons"], xlabel = "Year", ylabel = "Bgal/yr", linecolor = [:green :blue :indigo :pink :red], 
            title = "Supply & Demand",legend=:outerright, legendtitle="Supply & Demand")
    end
    if(p[24]==1)
        vline!([year_CAP],labels="t_CAP",linestyle=:dash, linecolor = :magenta)
    end
    
    #Flows
    if(p[24] == 1)
        plt_flows = plot(vars.year, [vars.Q_a.*0.001 vars.Q_1.*0.001 vars.Q_2.*0.001 vars.Q_a_g.*0.001 vars.O_d.*0.001 vars.O_1.*0.001 vars.O_2.*0.001 vars.O_g.*0.001], 
            labels=["In_all" "In_SRP" "In_CAP" "In_GW" "Use_all" "Use_SRP" "Use_CAP" "Use_GW"], ylims = (0,1.2*0.001*vars.Q_a[1]), xlabel = "Year", ylabel = "KAFY", title = "Inflows & Use", legend = :outerright, 
            linecolor = [:black :darkgreen :purple :turquoise :black :darkgreen :purple :turquoise], linestyle=[:solid :solid :solid :solid :dot :dot :dot :dot], legendtitle="Inflows & Use")
        vline!([year_CAP],labels="t_CAP",linestyle=:dash, linecolor = :magenta)
    else
        if(units == "AF")
            plt_flows = plot(vars.year, [vars.Q_a.*0.001 vars.Q_a_s.*0.001 (1 .+ vars.C_v_s).*vars.μ_s.*p[13][4].*0.001 (1 .- vars.C_v_s).*vars.μ_s.*p[13][4].*0.001 vars.Q_a_g.*0.001 (1 .+ vars.C_v_g).*vars.μ_g.*p[13][2].*0.001 (1 .- vars.C_v_g).*vars.μ_g.*p[13][2].*0.001 vars.O_d.*0.001 vars.O_s.*0.001 vars.O_g.*0.001 vars.O_f.*0.001], 
                labels=["In_all" "In_s" "+σ_s" "-σ_s" "In_g" "+σ_g" "-σ_g" "Use_all" "Use_s" "Use_g" "Flood"], ylims = (0,1.5*0.001*vars.μ[1]), xlabel = "Year", title = "Inflows & Use", ylabel = "KAFY", legend = :outerright, 
                linecolor = [:black :green :green :green :turquoise :turquoise :turquoise :black :green :turquoise :brown], linestyle=[:solid :solid :dash :dash :solid :dash :dash :dot :dot :dot :dot], legendtitle="Inflows & Use")
        else
            plt_flows = plot(vars.year, [vars.Q_a.*0.000000001 vars.Q_a_s.*0.000000001 (1 .+ vars.C_v_s).*vars.μ_s.*p[13][4].*0.000000001 (1 .- vars.C_v_s).*vars.μ_s.*p[13][4].*0.000000001 vars.Q_a_g.*0.000000001 (1 .+ vars.C_v_g).*vars.μ_g.*p[13][2].*0.000000001 (1 .- vars.C_v_g).*vars.μ_g.*p[13][2].*0.000000001 vars.O_d.*0.000000001 vars.O_s.*0.000000001 vars.O_g.*0.000000001 vars.O_f.*0.000000001], 
                labels=["In_all" "In_s" "+σ_s" "-σ_s" "In_g" "+σ_g" "-σ_g" "Use_all" "Use_s" "Use_g" "Flood"], ylims = (0,1.5*0.000000001*vars.μ[1]), xlabel = "Year", title = "Inflows & Use", ylabel = "Bgal.yr", 
                legend = :outerright, linecolor = [:black :green :green :green :turquoise :turquoise :turquoise :black :green :turquoise :brown], 
                linestyle=[:solid :solid :dash :dash :solid :dash :dash :dot :dot :dot :dot], legendtitle="Inflows & Use")
        end
    end
    
    #Storage Volume & Capacity
    if(p[24] ==1)
        SY = p[13][2]*p[1][2]*100 #calculate the safe-yield storage volume
        
        plt_stor = plot(vars.year, [vars.V_g.*0.000001 vars.Vbar_g.*0.000001 vars.υ_proj_g.*vars.Vbar_g.*0.000001], labels = ["Stor_Vol" "Stor_Capac" "Stor_Vol_proj"], ylims = (0,vars.Vbar_g[1]*0.000001*1.2), 
                xlabel = "Year", ylabel = "MAF", title = "Local Aquifer Storage", legend = :outerright, linecolor = [:grey :black :brown], legendtitle="Storage Capacity")
        vline!([year_CAP],labels="t_CAP",linestyle=:dash, linecolor = :magenta)
        hline!([SY.*0.000001],labels="SY",linestyle=:dash,linecolor=:red)
    else
        if(units == "AF")
            plt_stor = plot(vars.year, [vars.V_s.*0.000001 vars.Vbar_s.*0.000001 vars.V_g.*0.000001 vars.Vbar_g.*0.000001], labels = ["Stor_Vol" "Stor_Capac" "Stor_Vol" "Stor_Capac"], 
                ylims = (0,max(p[12][1]*vars.C_v_s[1]*vars.μ_s[1], p[12][2]*vars.C_v_g[1]*vars.μ_g[1])*0.000001*1.5), xlabel = "Year", ylabel = "MAF", title = "Storage", 
                legend = :outerright, linecolor = [:grey :black :grey :black], legendtitle="Storage Capacity")
        else
            plt_stor = plot(vars.year, [vars.V_s.*0.000000001 vars.Vbar_s.*0.000000001 vars.V_g.*0.000000001 vars.Vbar_g.*0.000000001], labels = ["Stor_Vol" "Stor_Capac" "Stor_Vol" "Stor_Capac"], 
                ylims = (0,max(p[12][1]*vars.C_v_s[1]*vars.μ_s[1], p[12][2]*vars.C_v_g[1]*vars.μ_g[1])*0.000000001*1.5), xlabel = "Year", ylabel = "Bgal", title = "Storage", 
                legend = :outerright, linecolor = [:grey :black :grey :black], legendtitle="Storage Capacity")
        end
    end
    
    #Error
    plt_e = plot(vars.year,[vars.e_1 vars.e_2 vars.e_3], labels = ["Short-Term" "Invest" "Rates"], xlabel = "Year", ylabel = "Error", title = "Error", 
        legend = :outerright, ylims = (-3,3), linecolor = [:blue :grey :green :brown], legendtitle="Error")
    if(p[24]==1)
        vline!([year_CAP],labels="t_CAP",linestyle=:dash, linecolor = :magenta)
    end
    
    
    #Attention
    plt_Y = plot(vars.year,[vars.Y_1 vars.Y_2 vars.Y_3],labels = ["Short-Term" "Long-Term" "Rates"], xlabel="Year", ylabel = "Attention", 
        title = "Attention", legend=:outerright, linecolor = [:blue :green :brown], legendtitle="Attention")
    if(p[24]==1)
        vline!([year_CAP],labels="t_CAP",linestyle=:dash, linecolor = :magenta)
    end

    #Financial Flows
    plt_fin = plot(vars.year,[vars.R.*0.000001 vars.C_o.*0.000001 vars.C_d.*0.000001 vars.J_bar.*0.000001 vars.J.*0.000001 vars.J_m_need.*0.000001], 
        labels = ["Rev" "Op Costs" "Debt Serv" "MaxInvest" "Inv Long" "Inv_Maint"], 
        xlabel="Year", ylabel = "Dollars (M)", title = "Financial Flows", legend=:outerright, ylims=(0,(maximum(vars.R)*0.000001*2)), legendtitle="Financial Flows", 
        linecolor = [:purple :red :orange :black :green :brown])
    if(p[24]==1)
        vline!([year_CAP],labels="t_CAP",linestyle=:dash, linecolor = :magenta)
    end
    
    #Rate-Setting
    plt_f = plot(vars.year,vars.f.*p[2], labels = "f", xlabel="Time(years)", ylabel = "Per-Capita Rates (Dollars/yr)", title = "Rates - Max Rates Ratio", legend=:outerright, 
        ylims=(0,1.2*maximum(vars.f.*p[2])), linecolor = :brown, legendtitle="Rates")
    hline!([p[2]], labels = "Rate_PC_max", linestyle=:dash, linecolor=:brown) 
    if(p[24]==1)
        vline!([year_CAP],labels="t_CAP",linestyle=:dash, linecolor = :magenta)
    end
    
    #Delivery Efficiency
    plt_η = plot(vars.year, [vars.η vars.u_impl_η./vars.A], xlabel = "Year",labels = ["η" "Invest_η"],ylabel="Del. Eff (%Outflow)", title = "Delivery Efficiency", 
        linecolor = :purple, linestyle = [:solid :dot], legend=:outerright, ylims = (0,p[11]*1.2), legendtitle="Delivery Efficiency")
    hline!([p[11]], labels="η_max", linestyle=:dash, linecolor=:purple)
    if(p[24]==1)
        vline!([year_CAP],labels="t_CAP",linestyle=:dash, linecolor = :magenta)
    end
    
    #Per-Capita Demand
    if(units == "AF")
        plt_d = plot(vars.year, [vars.d_ST.*892.15 vars.χbar.*vars.μ.*892.15 vars.u_1.*(-1).*vars.μ.*892.15 vars.u_impl_dbar.*vars.μ.*(-1).*892.15], xlabel = "Year", ylabel="PC Dem (GPCD)",
            labels=["Dem_PC" "BaseDem_PC" "Conserv_ST" "Conserv_LT"], title="Per-Capita Demand", ylims = (0,1.1*vars.χbar[1]*892.15*vars.μ[1]), legend=:outerright, linecolor = [:blue :indigo :blue :indigo], 
            linestyle = [:solid :solid :dot :dot], legendtitle="Per-Capita Demand")
        hline!([p[14]*892.15], labels="Dem_PC_min", linestyle=:dash, linecolor=:indigo)
    else
       plt_d = plot(vars.year, [vars.d_ST./365 vars.χbar.*vars.μ./365 vars.u_1.*(-1)./365 vars.u_impl_plan_dbar.*(-1)./365], xlabel = "Year", ylabel="PC Dem (GPCD)", 
            labels=["Dem_PC" "BaseDem_PC" "Conserv_ST" "Conserv_LT"], title="Per-Capita Demand", ylims = (0,1.1*vars.χbar[1]/365*vars.μ[1]), legend=:outerright, linecolor = [:blue :indigo :blue :indigo], 
            linestyle = [:solid :solid :dot :dot], legendtitle="Per-Capita Demand")
        hline!([p[14]/365], labels="Dem_PC_min", linestyle=:dash, linecolor=:indigo)
    end
    if(p[24]==1)
        vline!([year_CAP],labels="t_CAP",linestyle=:dash, linecolor = :magenta)
    end
    
    #Pumping Capacity
    if(units == "AF")
        plt_w = plot(vars.year, [vars.w_s.*(vars.υ_bar_s.*vars.C_v_s .+ 1).*vars.μ_s.*0.001 vars.w_g.*(vars.υ_bar_g.*vars.C_v_g .+ 1).*vars.μ_g.*0.001 vars.u_impl_w_s.*0.001 vars.u_impl_w_g*0.001 vars.w_max_s.*(vars.υ_bar_s.*vars.C_v_s .+ 1).*vars.μ_s.*0.001 vars.w_max_g.*(vars.υ_bar_g.*vars.C_v_g .+ 1).*vars.μ_g.*0.001], 
            xlabel = "Year", ylabel = "KAF/yr", labels = ["w_s" "w_g" "Invest_w_s" "Invest_w_g" "w_s_max" "w_g_max"], title = "Processing Capacity", linecolor = [:green :turquoise :green :turquoise :green :turquoise], 
            linestyle = [:solid :solid :dot :dot :dash :dash], legend=:outerright, ylims = (0, max(maximum(vars.A_l_s*1.2*0.001),maximum(vars.A_l_g.*1.2*0.001))), legendtitle="Processing Capacity")
    else
        plt_w = plot(vars.year, [vars.w_s.*(vars.υ_bar_s.*vars.C_v_s .+ 1).*vars.μ_s.*0.000000001 vars.w_g.*(vars.υ_bar_g.*vars.C_v_g .+ 1).*vars.μ_g.*0.000000001 vars.u_impl_w_s.*0.000000001 vars.u_impl_w_g*0.000000001 vars.w_max_s.*(vars.υ_bar_s.*vars.C_v_s .+ 1).*vars.μ_s.*0.000000001 vars.w_max_g.*(vars.υ_bar_g.*vars.C_v_g .+ 1).*vars.μ_g.*0.000000001], 
            xlabel = "Year", ylabel = "Bgal/yr", labels = ["w_s" "w_g" "Invest_w_s" "Invest_w_g" "w_s_max" "w_g_max"], title = "Processing Capacity", linecolor = [:green :turquoise :green :turquoise :green :turquoise], 
            linestyle = [:solid :solid :dot :dot :dash :dash], legend=:outerright, ylims = (0, max(maximum(vars.w_max_s.*(vars.V_s .+ vars.μ_s)),maximum(vars.w_max_g.*(vars.Vbar_g .+ vars.μ_g)))*0.000000001*1.2), legendtitle="Processing Capacity")
    end
    if(p[24]==1)
        vline!([year_CAP],labels="t_CAP",linestyle=:dash, linecolor = :magenta)
    end
    
    #Mean Inflow
    if(units=="AF")
        plt_μ = plot(vars.year, [vars.μ_s.*0.001 vars.μ_g.*0.001], xlabel = "Year", labels = ["μ_s" "μ_g"], ylabel = "Mean Inflow (KAFY)", title = "Mean Inflow", linecolor = [:darkgreen :turquoise], linestyle = :solid, 
        legend = :outerright, ylims = (0,max(p[1][1],p[1][2])*1.2*0.001), legendtitle="Mean Inflow")
    else
        plt_μ = plot(vars.year, [vars.μ_s.*0.000000001 vars.μ_g.*0.000000001], xlabel = "Year", labels = ["μ_s" "μ_g"], ylabel = "Mean Inflow (Bgal.yr)", title = "Mean Inflow", linecolor = [:darkgreen :turquoise], 
            linestyle = :solid, legend = :outerright, ylims = (0,max(p[1][1],p[1][2])*1.2*0.000000001), legendtitle="Mean Inflow")
    end
    if(p[24]==1)
        vline!([year_CAP],labels="t_CAP",linestyle=:dash, linecolor = :magenta)
    end
    
    #Population
    plt_P = plot(vars.year, vars.P.*0.001, ylabel="Persons (K)", xlabel = "Year", title = "Population", labels="Pop", 
        ylims=(0,p[5]*1.2*0.001), linecolor = :orange, linestyle = [:solid :dash :dot], legend = :outerright, legendtitle="Population")
    if(p[24]==1)
        vline!([year_CAP],labels="t_CAP",linestyle=:dash, linecolor = :magenta)
    end
    
    #Aggregate Plots
    plt = plot(plt_short, plt_flows, plt_stor, plt_e, plt_Y, plt_fin, plt_f, plt_η, plt_d, plt_w, plt_P, plt_μ, size = (1800,1500), layout = (4,3))
    
    #Save in Plot List
    plt_list = [plt_short plt_flows plt_stor plt_e plt_Y plt_fin plt_f plt_η plt_d plt_w plt_P plt_μ plt]
    
    return plt_list
end;

#   c.) Aggregate Output
#   ––––––––––––––––––––––

function UWIIM_output(model; setup=Default(), t_run=100, year_0=2010, units="AF")
    p=setup[1]
    x_0=Float64.(setup[2])
    
    varsDF = createVarsDF(model,p,t_run,year_0)
    
    plt_ts_nd = timeSeriesPlot_nondim(varsDF,p,x_0)
    
    plt_ts_dim = timeSeriesPlot_dim(varsDF,p,x_0,units)
    
    output = [varsDF, plt_ts_nd, plt_ts_dim]
    
    return output
end;

#   3.3 One Line to Rule Them All
#   ===============================

#   In addition to defining the setup for the model, the user has additional
#   settings they can set for a given model run
# 
#     1. t_run (default = 100): the number of years the trajectory should
#        run
# 
#     2. year_0 (default = 2010): for plotting, define what the initial
#        year is
# 
#     3. units (default = "AF"): specify "AF" or "gal" for water units

function run_UWIIM(setup::Any=Default(); t_run=50, year_0=2010, units="AF", f=f)
    model = create_UWIIM(setup;f=f)
    
    output = UWIIM_output(model; setup=setup, t_run=t_run, year_0=year_0, units=units)
    
    return output
end;

#   4. Example Ouptputs
#   ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

#   To run any of the setups below, simply uncomment the line of code and add
#   whatever parameter or initial conditions changes you desire (see tables in
#   ReadMe for guidance on changing parameter and initial conditions)

#   With Phoenix Setup
#   ====================

#Generate Output (based on 1 trajectory with Phoenix setup) 
#output_test_PHX = run_UWIIM(Phoenix(); t_run=50,year_0=2010,units="AF");
#output_test_PHX[3][end] #this line shows the dimensional (water units of AF) time series for the model run

#   With Scottsdale Setup
#   =======================

#Generate Output (based on 1 trajectory with Phoenix setup)
#output_test_Sc = run_UWIIM(Scottsdale(); t_run=50,year_0=2010,units="AF")
#output_test_Sc[3][end] #this line shows the dimensional (water units of AF) time series for the model run

#   With Queen Creek Setup
#   ========================

#Generate Output
#output_test_QC = run_UWIIM(QueenCreek(Δμ_s_pc=0); t_run=50,year_0=2010,units="AF")
#output_test_QC[3][end] #this line shows the dimensional (water units of AF) time series for the model run