#   Read Me
#   ≡≡≡≡≡≡≡≡≡
# 
#   This model version is set up to allow for you to use without fear of messing
#   with the main model file. Section 1 contains the dynamic model definition.
#   For a detailed description of the equations and function definitions used,
#   refer to Section 1. For additional explanation, see the Tables and
#   description in the "UWIIM General Model Paper" (Overleaf). All state
#   variables and their dynamic equations are expressed in the "Equations of
#   Motion" in section 1.4.
# 
#   Section 2 contains information for how the model is initialized, including
#   the setup for the general default city, and city-specific setups. Each setup
#   is defined in terms of a function that takes as its input, alterations to
#   the originally intended parameters and initial conditions.
# 
#   Section 3 contains information for how to run the model, including the
#   SINGLE function that allows one to run a model with a certain defined setup
#   (see Section 2) and receive a report (output) on a single generated
#   trajectory. This output contains a (1) dataframe of all state and auxiliary
#   variables of interest over time, (2) multiple time series plots of the
#   variables in their reduced dimensional form, (3) multiple time series plots
#   of the variables in their common dimensional form (i.e. AFY), and (4) phase
#   plots of the major state variables in their reduced dimensional form.
# 
#   Sections 4 and 5 contain example outputs for the Phoenix and Indy scenarios.

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

#   1.1 Definitions of Auxilary Water Variables (Plant)
#   =====================================================

#   Mean Inflow (\mu)
#   –––––––––––––––––––

# :\
# 
#   mut = \mu^st + \mu^g_t + \mu^p $

μ(x,p,t) = x[11]+x[12]+p[1][3];

#   Streamflow, Inflow (Q)
#   ––––––––––––––––––––––––

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
#   = q^st\mu^st + q^gt\mu^gt + q^p_t\mu^p $

function Q(x,p,t)
    return Q_s(x,p,t) + Q_g(x,p,t) + Q_p(x,p,t)
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

#   Purchased Water Inflows (Q^p)
#   -------------------------------

# :Q
# 
#   ^pt = O^pt $
# 
# :q
# 
#   ^pt = \frac{Q^pt}{\mu^p} $

function Q_p(x,p,t)
    return O_p(x,p,t)
end;

function q_p(x,p,t)
    return ifelse(p[1][3]==0, 0, Q_p(x,p,t)/p[1][3])
end;

#   Available Inflows (Q^a)
#   -------------------------

#   Surface Available Inflows
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅

# :Q
# 
#   ^{a,s} = \begin{cases} a^{SRP}\mu^1t + a^{CAP}(Q^st - \mu^1) & \text{if}
#   \quad cases = 1 \ a^{s,q}Q^s_t & \text{otherwise} \end{cases} $

function Q_a_s(x,p,t) 
    if(p[24]==1)
        return Q_a_SRP(x,p,t) + Q_a_CAP(x,p,t)
    else
        return p[13][4]*Q_s(x,p,t)
    end
end;

#   CAP and SRP Available Inflows
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅

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
    Q_CAP = Q_s(x,p,t) - p[1][4]
    short = p[1][5] - Q_CAP
    NIA_avail = max(0,70022-short)
    high_avail = Q_CAP - NIA_avail
    
    return p[13][6]*NIA_avail + p[13][7]*high_avail
end;

# :\
# 
#   tilde{Q}^{SRP}t = min(\mu1a^{q,1}, \frac{Dt}{\etat} \xi1 + \mu1 a^{q,NCS}) $

function Q_a_SRP(x,p,t)
    max_prop_use = (D_bar(x,p,t)/x[15])*p[27][2] #ifelse(p[27][2]==0,0,(((1/p[27][2])-1)^(-1))*(A_l_CAP + A_g_t))
    NCS = p[13][8]*p[1][4]
    max_allocation = p[1][4]*p[13][4]
    
    return min(max_allocation, max_prop_use + NCS)  
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

Q_a(x,p,t) = Q_a_s(x,p,t) + Q_a_g(x,p,t) + Q_p(x,p,t);

# :q
# 
#   ^at = \frac{\tilde{Q}t}{\mu_t} $

q_a(x,p,t) = Q_a(x,p,t)/μ(x,p,t);

#   City Portion of Inflow (a^q)
#   ------------------------------

#   Default Version (All inflow going to the city)
# 
# :a
# 
#   ^q = \frac{\tilde{Q}}{\mu_t} $

function a_q(x,p,t)
    return Q_a(x,p,t)/μ(x,p,t)
end;

#   System Population (P)
#   –––––––––––––––––––––––

#   Total Population (P)
#   ----------------------
# 
#   Default $ Pt = pt \kappa_t $

P(x,p,t) = x[1]*κ(x,p,t);

#   Water-based Carrying Capacity (Kbar)
#   --------------------------------------
# 
# :\
# 
#   bar{\kappa}t = \frac{1}{\bar{\chi}t} a^q $

κbar(x,p,t) = (1/x[16])*a_q(x,p,t);

#   Actual Carrying Capacity (\kappa)
#   -----------------------------------
# 
#   Default (exogenous \kappa)
# 
# :\
# 
#   kappat = \begin{cases} \kappa & \text{if} \quad \kappa{endog} = 0 \
#   \bar{\kappa} & \text{otherwise} \end{cases} $

function κ(x,p,t)
    if(p[5][1]==0)
        return p[5][2]
    else
        return κbar(x,p,t)
    end
end;

#   Projected Population (P_{proj})
#   ---------------------------------
# 
#   *For Long-Term Investment Controller if Projection Setting is Turned On
# 
# :P
# 
#   ^{proj}t = \begin{cases} \frac{Kt}{\frac{Kt-Pt}{Pt}exp(-r\taup) + 1} &
#   \text{if} \quad P{assume} = 0 \ \frac{Kt(r+mi)}{r - \frac{rPt -
#   Kt(r+mi)}{Pt} exp(-(r+mi)\tau_p)} & \text{otherwise} \end{cases} $

function P_proj(x,p,t)
    P_t = P(x,p,t)
    κ_t = κ(x,p,t)
    
    if(p[29][4]==0)
        return κ_t/(((κ_t-P_t)/P_t)*exp(-p[6]*p[29][2]) + 1)
    else
        num = κ_t*(p[6]+p[7][2])
        denom = p[6] - ((p[6]*P_t - num)/P_t)*exp(-(p[6]+p[7][2])*p[29][2])
        
        return num/denom
    end
end;

function P_proj_1(x,p,t)
    P_t = P(x,p,t)
    κ_t = κ(x,p,t)
    
    if(p[29][4]==0)
        return κ_t/(((κ_t-P_t)/P_t)*exp(-p[6]) + 1)
    else
        num = κ_t*(p[6]+p[7][2])
        denom = p[6] - ((p[6]*P_t - num)/P_t)*exp(-(p[6]+p[7][2]))
        
        return num/denom
    end
end;

#   System Demand (D)
#   –––––––––––––––––––

# :D_t
# 
#   = dt Pt = \chit \mut P_t $
# 
# :\
# 
#   bar{D}t = \bar{d}t Pt = \bar{\chi}t \mut Pt $

D(x,p,t) = x[2]*P(x,p,t)*μ(x,p,t);

D_bar(x,p,t) = x[16]*P(x,p,t)*μ(x,p,t);

#   Post-ST Conservation Demand
#   -----------------------------

d_ST(x,p,t) = (x[2]-H_d(x,p,t))*μ(x,p,t);

D_ST(x,p,t) = d_ST(x,p,t)*P(x,p,t);

#   Projected Demand
#   ------------------
# 
#   \hat{D}_t = \hat{P}\bar{d}_t

D_proj(x,p,t) = P_proj(x,p,t)*x[16]*μ(x,p,t)*(1-p[15][2])^p[29][2];

D_proj_1(x,p,t) = P_proj_1(x,p,t)*x[16]*μ(x,p,t)*(1-p[15][2]);

#   Physical Volumes (V & \bar{V})
#   ––––––––––––––––––––––––––––––––

#   Storage Capacity
#   ------------------

# :\
# 
#   bar{V}t = \bar{\upsilon}t c^v{t} \mut $

Vbar_s(x,p,t) = x[5]*x[13]*x[11];

Vbar_g(x,p,t) = x[6]*x[14]*x[12];

#   Volume
#   --------

# :V_t
# 
#   = \upsilont \bar{V}t $

V_s(x,p,t) = x[3]*Vbar_s(x,p,t);

V_g(x,p,t) = x[4]*Vbar_g(x,p,t);

#   Projected Volume (\upsilon_{proj})
#   ------------------------------------

function υ_proj_s(x,p,t)
    return x[3]
end;

function υ_and_η_proj_g(x,p,t)
    A_l_s_now = A_l_s(x,p,t)
    Inflow = Q_a_g(x,p,t) + Q_b(x,p,t)
    
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
    A_t_s_now = A_t_s(x,p,t)
    A_t_s_prev = copy(A_t_s_now)
    A_t_s_next = copy(A_t_s_now)
    A_prev=0
    
    if(p[29][2]>0)
        #Number of investments that can be implemented in the projection period
        possibleInvests_η = floor(Int,min(p[29][2],p[16][2]-1)) 
        possibleInvests_w_s = floor(Int,min(p[29][2],p[16][4]-1)) 
        
        #Loop for Each Projection Year and Calculate Projected Groundwater Use
        for y in 1:p[29][2]
            #0: Add Groundwater Use to Sum
            O_proj += O_g_next
            
            #1: New Surface Water Processing Availability
            A_t_s_prev = copy(A_t_s_next)
            if(possibleInvests_w_s>0)
                A_t_s_next += x[plan_index_k(x,p,t)[4]+(y-1)]
                possibleInvests_w_s -= 1
            end
            
            #2: New Surface Water Availability
            new_SW_prev = copy(new_SW)
            O_s_add_prev = copy(O_s_add_next)
            if(A_l_s_now>=A_t_s_next)
                new_SW = A_t_s_next - A_t_s_prev
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
                H_impl_η=x[plan_index_k(x,p,t)[2]+(y-1)]
                η_next = η_prev + (H_impl_η/A_prev) 
                possibleInvests_η -= 1
                
                η_next = ifelse(η_next>p[11],p[11],η_next)
            else
                η_next = copy(η_prev)
            end
            
            #5: New Groundwater Use
            O_g_prev = copy(O_g_next)
            O_g_next = (η_prev/η_next)*(O_g_prev + O_s_add_prev) - O_s_add_next
        end
        
        υ_proj_t=max(100*p[13][2]*p[1][2],(V_g(x,p,t) + Inflow*p[29][2] - O_proj))/Vbar_g(x,p,t)
        η_proj_t=η_next
        return [υ_proj_t η_proj_t]
    else
        return [x[4] x[15]]
    end
end;

#   Available Volume (A)
#   ––––––––––––––––––––––

#   Total Available
#   -----------------

# :A_t
# 
#   = A^st + A^gt + A^p_t $

A(x,p,t) = A_s(x,p,t) + A_g(x,p,t) + A_p(x,p,t);

#   Available by Types of Sources
#   -------------------------------

# :A
# 
#   ^st = min(A^{l,s}t + A^{T,s}_t) $

A_s(x,p,t) = min(A_l_s(x,p,t), A_t_s(x,p,t));

# :A
# 
#   ^gt = min( A^{l,g}t + A^{T,g}_t) $

A_g(x,p,t) = min(A_l_g(x,p,t), A_t_g(x,p,t));

# :A
# 
#   ^p_t = \bar{\mu}^p $

A_p(x,p,t) = p[1][3];

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
    if(p[34]==0)
        if(p[24]==1)
            A_l_g_v = p[13][1]*max(0,V_g(x,p,t)-p[13][2]*p[1][2]*100) ##only count available V if above 100-yr SY amount 
        else
            A_l_g_v = p[13][1]*V_g(x,p,t)
        end
    else
        A_l_g_v = V_g(x,p,t)
    end
    
    
    A_l_g_q = p[13][2]*Q_g(x,p,t);
    
    return A_l_g_v + A_l_g_q
end;

#   Technically Available Volume for Withdrawal (A^T)
#   ---------------------------------------------------

#   The technological limit is steered by processing capacity (i.e., treatment,
#   pumping, etc.). This specifies the limit to what the system can process from
#   its physical inflow.
# 
# :A
# 
#   ^Tt = wt (\bar{V} + \mu_t) $

A_t_s(x,p,t) = x[7]*(Vbar_s(x,p,t) + x[11]);

A_t_g(x,p,t) = x[8]*(Vbar_g(x,p,t) + x[12]);

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

function η_proj(x,p,t)
    return υ_and_η_proj_g(x,p,t)[2]
end

function Vbar_proj_s(x,p,t)
    Vbar_t = Vbar_s(x,p,t)
    
    ##Determine the number of planned investments that will be implemented by the end of the projection period
    if(p[29][2]==0) #no projection 
        return Vbar_t
    elseif(p[16][3]==1) #immediate implementation
        return Vbar_t
    else
        possibleInvests = floor(Int,min(p[29][2],p[16][3]-1)) 
        if(possibleInvests==1)#only implement 1 planned investment
            return Vbar_t + x[plan_index_k(x,p,t)[3]]
        else
            return Vbar_t + sum(x[plan_index_k(x,p,t)[3]:(plan_index_k(x,p,t)[3]+possibleInvests-1)])
        end
    end
end;

function υbar_proj_s(x,p,t)
    return Vbar_proj_s(x,p,t)/(x[13]*x[11])
end

function w_proj_s(x,p,t)
    Vbar_t = Vbar_s(x,p,t)
    Vbar_proj_t = Vbar_proj_s(x,p,t)
    μ_t = copy(x[11])
    μ_proj_t = μ_t
    A_w_t = (x[7]*(Vbar_t+μ_t))
    
    ##Determine the number of planned investments that will be implemented by the end of the projection period
    if(p[29][2]==0) #no projection 
        return x[7]
    elseif(p[16][4]==1) #immediate implementation
        return A_w_t/(Vbar_proj_t + μ_proj_t)
    else
        possibleInvests = floor(Int,min(p[29][2],p[16][4]-1)) 
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
    if(p[29][2]==0) #no projection 
        return x[8]
    elseif(p[16][5]==1) #immediate implementation
        return A_w_t/(Vbar_proj_t + μ_proj_t)
    else
        possibleInvests = floor(Int,min(p[29][2],p[16][5]-1))  
        if(possibleInvests==1)#only implement 1 planned investment
            return (A_w_t + x[plan_index_k(x,p,t)[5]])/(Vbar_proj_t + μ_proj_t)
        else
            return (A_w_t + sum(x[plan_index_k(x,p,t)[5]:(plan_index_k(x,p,t)[5]+possibleInvests-1)]))/(Vbar_proj_t + μ_proj_t)
        end
    end
end;

#   Total Projected Available Volume
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅

#   Same format as A_t but need to project \upsilon. Assume available purchased
#   water remains the same.
# 
# :A
# 
#   ^{proj}t = min(A^{l,proj}t, A^{T,proj}_t) $

A_proj(x,p,t) = A_proj_s(x,p,t) + A_proj_g(x,p,t) + A_p(x,p,t);

A_proj_s(x,p,t) = min(A_proj_l_s(x,p,t), A_proj_T_s(x,p,t));

A_proj_g(x,p,t) = min(A_proj_l_g(x,p,t), A_proj_T_g(x,p,t));

#   Projected Legally Available Volume (A^{proj,l})
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅

#   Follows the same algorithm as the actual legal availability assuming that
#   planned investments are implemented, the projected groundwater volume, and
#   the same mean inflow

function A_proj_l_s(x,p,t) 
    A_proj_l_sv = p[13][3]*υ_proj_s(x,p,t)*Vbar_proj_s(x,p,t)
    
    if(p[24]==1) #Phoenix Scenario
        #CAP
        Q_CAP = x[11] - p[1][4]
        short = p[1][5] - Q_CAP
        NIA_avail = max(0,70022-short)
        high_avail = Q_CAP - NIA_avail
        CAP_avail = p[13][6]*NIA_avail + p[13][7]*high_avail
        
        #SRP
        max_prop_use = (D_proj(x,p,t)/η_proj(x,p,t))*p[27][2]
        NCS = p[13][8]*p[1][4]
        max_allocation = p[1][4]*p[13][4]
        
        SRP_avail = min(max_allocation, max_prop_use + NCS) 
        
        return A_proj_l_sv + SRP_avail + CAP_avail
    elseif(p[29][3]==0)
        return  A_proj_l_sv + p[13][4]*x[11]
    else
        return  A_proj_l_sv + p[13][4]*x[11]*(1-x[13])
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
    if(p[29][3]==0)
        if(p[24]==1)
            A_l_g_q = p[13][2]*x[12]
        else
            A_l_g_q = p[13][2]*x[12];
        end
    else
        A_l_g_q = p[13][2]*x[12]*(1-x[14]);
    end
    
    return A_l_g_v + A_l_g_q
end; 

#   Projected Maximum Processing Capacity

function w_max_g_proj(x,p,t)
    A_l_g_t = A_proj_l_g(x,p,t)
    w_max_g = (A_l_g_t)/(Vbar_g(x,p,t)+x[12]) 
    
    return w_max_g
end;

#   Projected Technically Available Volume (A^{proj,T})
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅

# :A
# 
#   ^{T,proj}t = \hat{w}t(\hat{\bar{V}}t + \mut) $

function A_proj_T_s(x,p,t)
    return w_proj_s(x,p,t)*(Vbar_proj_s(x,p,t)+x[11])
end;

A_proj_T_g(x,p,t) = w_proj_g(x,p,t)*(Vbar_g(x,p,t)+x[12]);

#   Outflows, Releases (O)
#   ––––––––––––––––––––––––

#   Total Use/Outflows (O)
#   ------------------------

#   Outflows either are used to satisfy demand using surface water (O^s_t),
#   groundwater (O^g_t), purchased water (O^p_t), or release flood water (O^f_t)
# 
# :O_t
# 
#   = O^dt + O^ft $

O(x,p,t) = O_d(x,p,t) + O_f(x,p,t);

#   Demand-Related Outflows (O^d)
#   -------------------------------

# :O
# 
#   ^dt = O^st + O^gt + O^pt $

O_d(x,p,t) = O_s(x,p,t) + O_g(x,p,t) + O_p(x,p,t);

#   Surface Water Demand-Related Outflows (O^s)
#   ---------------------------------------------

# :O
# 
#   ^st = \begin{cases} O^1t + O^2t & \text{if} \quad cases = 1 \ min(A^st,
#   \frac{\tilde{D}t}{\etat}(1-\theta^g) & \text{otherwise} \end{cases} $

function O_s(x,p,t)
    if(p[24]==1)
        return O_1(x,p,t) + O_2(x,p,t)
    else
        need = (D_ST(x,p,t)/x[15])*(1-p[27][1])
        return min(A_s(x,p,t), need)
    end
end;

#   SRP Use (PHX Version)
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅

# :O
# 
#   ^1t = min(\mu1a^{SRP}, A^{T,s}t, \frac{\tilde{D}t}{\etat}\xi^{SRP} +
#   \mu1a^{NCS}, \frac{\tilde{D}t}{\etat}(1-\xi^{g})) $

function O_1(x,p,t) 
    #calculate need after considering annual groundwater use
    need = (D_ST(x,p,t)/x[15])*(1-p[27][1])
    
    #calculate available considering legal and technical
    tech_avail = A_t_s(x,p,t)
    allocation = p[1][4]*p[13][4]
    dem_op = (D_ST(x,p,t)/x[15])*p[27][2] + p[13][8]*p[1][4]
    
    #return min of need and available
    return min(allocation, tech_avail, dem_op, need)
end;

#   CAP Use (PHX Version)
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅

# :O
# 
#   ^2t = min(\tilde{Q}^{CAP}t, \frac{\tilde{D}t}{\etat} - O^1t, A^{T,s}t) $

function O_2(x,p,t)
    #calculate need
    need =(D_ST(x,p,t)/x[15])*(1-p[27][1])
    need_left = need - O_1(x,p,t)
    
    #calculate available considering legal and technical
    legal_avail = Q_a_CAP(x,p,t)
    tech_avail = A_t_s(x,p,t)
    
    return min(legal_avail, tech_avail, need_left)
end;

#   Banked Water (Q^b) (PHX Version)
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅

# :Q
# 
#   ^b = \begin{cases} (Q^st - \mu^1)a^{CAP} - O^2t & \text{if} \quad cases = 1
#   \ 0 & \text{otherwise} \end{cases} $

function Q_b(x,p,t) 
    if(p[24]==1)
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
    need_left = need - O_s(x,p,t)
    
    return min(need_left, A_g(x,p,t))
end;

#   Purchased Water Demand-Related Outflows (O^p_t)
#   -------------------------------------------------

# :O
# 
#   ^pt = min(A^pt, \frac{Dt}{\etat} - (O^st + O^gt)) $

function O_p(x,p,t) 
    need = D_ST(x,p,t)/x[15]
    need_left = need - (O_s(x,p,t) + O_g(x,p,t))
    
    return min(A_p(x,p,t),  need_left);
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

#   System Supply (S)
#   –––––––––––––––––––

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

#   1.2 Signal & Error Calculations
#   =================================

#   Signal - Safety Factor - Short-Term Conservation
#   ––––––––––––––––––––––––––––––––––––––––––––––––––
# 
# :SF_t
# 
#   = \frac{St}{Dt} $

M_sf(x,p,t) = S(x,p,t)/D(x,p,t);

#   Signal - Safety Fator - Baseline Demand, Long-Term Projection
#   –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

M_sf_l(x,p,t) = S_proj(x,p,t)/D_proj(x,p,t);

#   Signal - Debt Service Capacity Ratio (DSCR)
#   –––––––––––––––––––––––––––––––––––––––––––––
# 
# :DSCR
# 
#   = \frac{Rt-C^ot}{\hat{C}^{d}_t} $
# 
# :\
# 
#   hat{C}^dt = \tilde{J}^bt(1+\taubib - \frac{1}{\taub} - ib) +
#   \hat{J}^bt(\frac{1}{\taub} + i_b) $

function C_d_need(x,p,t)
    return x[18]*(1+p[31]*p[32]-(1/p[31])-p[32]) + J_b(x,p,t)*((1/p[31])+p[32])
end;

function M_dscr(x,p,t)
    C_d_need_t = C_d_need(x,p,t)
    
    return (R_proj(x,p,t) - C_o_proj(x,p,t))/C_d_need_t
end;

#   Error - Short-Term (Proportional)
#   –––––––––––––––––––––––––––––––––––
# 
#   Proportional Error Just to Meet Demand (goal = 1)
# 
# :e_t
# 
#   ^s = 1-SF_t $

e_s(x,p,t) = p[17][1] - M_sf(x,p,t);

#   Error - Long-Term (Integral)
#   ––––––––––––––––––––––––––––––

#   Long-Term Error
#   -----------------

#   Proportional Error with reference to long-term buffer (\gamma). If using the
#   projection form of long-term error, the safety factor is projected out
#   according to \tau_p instead of using the current safety factor
# 
# :e
# 
#   ^lt = \gamma - \hat{SF}t $

e_l(x,p,t) = p[17][2] - M_sf_l(x,p,t);

#   Integrating Long-Term Error
#   -----------------------------

#   Now, add the proportional error to the memory of past error to create the
#   integral of long-term error
# 
# :e
# 
#   ^it = e^i{t-1} + e_t^l $

function e_i(x,p,t)
    e_l_t = e_l(x,p,t)
    
    return e_l_t + ifelse(p[28][1] == 0, x[planInvestMemorySize(p) + 19], sum(x[(planInvestMemorySize(p) + 20):(length(x)-p[28][2])]))
end;

#   Error - Rate-Setting (Proportional)
#   –––––––––––––––––––––––––––––––––––––
# 
#   *Implied that goal rate-need ratio is 1.
# 
# :e
# 
#   ^rt = 1 - RNt $

function e_r(x,p,t)
    #if(p[33]==1)
    #    return abs(p[17][3] - M_dscr(x,p,t))
    #else
        return p[17][3] - M_dscr(x,p,t)
    #end
end;

#   1.3 Definition of Controller Steps
#   ====================================

#   1.3.1 Attention
#   –––––––––––––––––

# :Y_t
# 
#   = \frac{1}{1 + \mathrm{exp}(-\lambda(e_t - \epsilon))} $

function Y_s(x,p,t)
    e_st = e_s(x,p,t)
    
    return 1/(1+exp(-p[20][1]*(e_st - p[21][1])))
end;

function Y_l(x,p,t)
    e_it = e_l(x,p,t)#e_i(x,p,t)
    
    return 1/(1+exp(-p[20][2]*(e_it - p[21][2])))
end;

function Y_r(x,p,t)
    e_rt = e_r(x,p,t)
    if(p[33]==1)
        return 1/(1+exp(-p[20][3]*(abs(e_rt) - p[21][3])))
    else
        return 1/(1+exp(-p[20][3]*(e_rt - p[21][3])))
    end
end;

#   1.3.2 Revenue & Cost Functions
#   ––––––––––––––––––––––––––––––––

#   Total Revenue Available to City (R)
#   -------------------------------------
# 
# :R_t
# 
#   = Pt\hat{\pi}t(\beta^{(\pi)}p +
#   (1-\beta^{(\pi)}p)\frac{\tilde{d}t}{\bar{d}t}) $

#R(x,p,t) = p[2]*x[17]*P(x,p,t);
R(x,p,t) = P(x,p,t)*p[2]*x[17]*(p[36]+((1-p[36])*(d_ST(x,p,t)/(x[16]*μ(x,p,t)))));

#   Projected Revenue (Next Year)
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅

R_proj(x,p,t) = P_proj_1(x,p,t)*x[17]*p[2];

#   Maxmimum Potential Revenue (R_max)
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅
# 
# :R_
# 
#   {max} = (1 + \psi^r) \hat{P}{t+1} \hat{\pi}t $

function R_max(x,p,t)
    if(x[17] > (1-p[18][2]))
        max_change = (1 - x[17])/x[17]
        return  P_proj_1(x,p,t)*(1+max_change)*p[2]*x[17]
    else
        return P_proj_1(x,p,t)*(1+p[18][2])*p[2]*x[17]
    end
end;

#   Operating Costs (C_0)
#   -----------------------
# 
#   The operating costs are a function of the system's infrastructure states,
#   service population size, and demand.
# 
# :C
# 
#   ^ot = go Pt^{z{op}}\tilde{D}t^{z{od}} $

C_o(x,p,t) = p[22][1]*(P(x,p,t)^p[23][1])*(D_bar(x,p,t)^p[23][2]);

#   Projected Operating Costs
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅

C_o_proj(x,p,t) = p[22][1]*(P_proj_1(x,p,t)^p[23][1])*(D_proj_1(x,p,t)^p[23][2]);

#   Debt Service (C_D)
#   --------------------
# 
#   C^d_t = \tilde{J}^b_{t}(1+\tau_bi_b)

C_d(x,p,t) = x[18]*(1+p[31]*p[32]);

#   1.3.3 Investment Action Situations/Controllers Response
#   –––––––––––––––––––––––––––––––––––––––––––––––––––––––––

#   1.3.3.1 Maintenance Investment Needs
#   --------------------------------------

#   Total Maintenance Investment Need
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅

function J_m_need(x,p,t)
    J_m_need = J_m_η(x,p,t) + J_m_v(x,p,t) + J_m_w_s(x,p,t) + J_m_w_g(x,p,t)
    return J_m_need
end;

function J_k_m_need(x,p,t)
    J_k_m = zeros(7)
    
    J_k_m[2] = J_m_η(x,p,t)
    J_k_m[3] = J_m_v(x,p,t)
    J_k_m[4] = J_m_w_s(x,p,t)
    J_k_m[5] = J_m_w_g(x,p,t)
    
    return J_k_m
end;

#   Delivery Efficiency Maintenance Investment Need
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅

function H_m_η(x,p,t)
    ###Adjust for ForwardDiff Objects
    x_t = ForwardDiff.value(x[15])
    A_t = ForwardDiff.value(A(x,p,t))
    
    x_goal = ifelse(x_t>p[11], p[11], x_t) #If the current capacity is above maximum, allow the capacity to decay to maximum
    
    H_m = (x_goal - x_t*(1-p[15][3]))*A_t
    
    if(H_m < 0)
        return 0
    else
        return H_m
    end
end;

function J_m_η(x,p,t)
    η_t = ForwardDiff.value(x[15])
    
    H_m = H_m_η(x,p,t)
    
    J_m = p[22][3]*((η_t*H_m)^p[23][4])
    
    return J_m
end;

#   Storage Capacity Maintenance Investment Need
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅

function H_m_v(x,p,t)
    ###Adjust for ForwardDiff Objects
    x_t = ForwardDiff.value(x[5])
    c_v_t = ForwardDiff.value(x[13])
    μ_t = ForwardDiff.value(μ(x,p,t))
    
    if(p[19][2]==0)
        return 0
    else
        x_goal = ifelse(x_t>p[12][1],p[12][1],x_t) #If the current capacity is above maximum, allow the capacity to decay to maximum
    
        H_m = (x_goal - x[5]*(1-p[15][4]))*μ_t*c_v_t
    
        if(H_m < 0)
            return 0
        else
            return H_m
        end
    
    end
end;

function J_m_v(x,p,t)
    if(p[19][2]==0)
        return 0
    else
        H_m = H_m_v(x,p,t)
        
        J_m= p[22][4]*(H_m^p[23][5])
        
        return J_m
    end
end;

#   Surface Processing Capacity Maintenance Investment Need
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅

function H_m_w_s(x,p,t)
    ###Adjust for ForwardDiff Objects
    x_t = ForwardDiff.value(x[7])
    w_max = ForwardDiff.value(w_max_s(x,p,t))
    V_bar_s_t = ForwardDiff.value(Vbar_s(x,p,t))
    μ_s_t = ForwardDiff.value(x[11])
    
    x_goal = ifelse(x_t>w_max, w_max, x_t) #If the current capacity is above maximum, allow the capacity to decay to maximum
    
    #No more τ_i because instantly implemented maintenance
    H_m = (x_goal - x_t*(1-p[15][5]))*(V_bar_s_t + μ_s_t)
    
    if(H_m < 0)
        return 0
    else
        return H_m
    end
end;

function J_m_w_s(x,p,t)
    H_m = H_m_w_s(x,p,t)
    
    J_m = p[22][5]*(H_m^p[23][6])
    
    return J_m
end;

#   Ground Processing Capacity Maintenance Investment Need
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅

function H_m_w_g(x,p,t)
    ###Adjust for ForwardDiff Objects
    x_t = ForwardDiff.value(x[8])
    w_max = ForwardDiff.value(w_max_g(x,p,t))
    V_bar_g_t = ForwardDiff.value(Vbar_g(x,p,t))
    μ_g_t = ForwardDiff.value(x[12])
    
    x_goal = ifelse(x_t>w_max, w_max, x_t) #If the current capacity is above maximum, allow the capacity to decay to maximum
    
    #No more τ_i because instantly implemented maintenance
    H_m = (x_goal - x_t*((1-p[15][6])))*(V_bar_g_t + μ_g_t)
    
    if(H_m < 0)
        return 0
    else
        return H_m
    end
end;

function J_m_w_g(x,p,t)
    H_m = H_m_w_g(x,p,t)
    
    J_m = p[22][6]*(H_m^p[23][7])
    
    return J_m
end;

#   1.3.3.2 Long-Term Investment
#   ------------------------------
# 
#   The amount of revenue allocated to long term investments (J_t) is a function
#   of the net available revenue after accounting for costs (R-C^o_t) and the
#   generated attention (Y).
# 
# :J_t
# 
#   = J(e{jt},x{jt}) \quad \text{where } j \in {l} $

#   Maximum Possible Investment (\bar{J})
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅

function J_bar(x,p,t)
    incr = (p[31]*(R_max(x,p,t)-C_o_proj(x,p,t)))/(p[17][3]*(1+p[31]*p[32]))
    
    NetOp = R(x,p,t)-C_o(x,p,t)

    J_bar = incr + NetOp - (1+p[32])*p[31]*x[18]
    
    return J_bar
end;

function J_bar_o(x,p,t)
    J_o_max = R(x,p,t) - C_o(x,p,t) - C_d(x,p,t)
    return J_o_max
end;

#   Supply Need (\hat{H})
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅

function H_m_need(x,p,t)
    return H_m_η(x,p,t) + H_m_v(x,p,t) + H_m_w_s(x,p,t) + H_m_w_g(x,p,t)
end;

function H_e_need(x,p,t)
    return max(Y_l(x,p,t)*e_l(x,p,t)*D_proj(x,p,t),0)
end;

function H_need(x,p,t)
    return H_m_need(x,p,t) + H_e_need(x,p,t)
end;

#   Re-distribute Excess Infrastructure Investments (\beta_{kt})
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅

#   function βk(x,p,t) ##Note Supply Needs Hneedt =
#   ForwardDiff.value(Heneed(x,p,t)) ####Now β just for expansion At =
#   ForwardDiff.value(A(x,p,t)) ηt = ForwardDiff.value(x[15]) μt =
#   ForwardDiff.value(μ(x,p,t)) Cvst = ForwardDiff.value(x[13]) υbarst =
#   ForwardDiff.value(x[5]) Vbarst = ForwardDiff.value(Vbars(x,p,t)) Vbargt =
#   ForwardDiff.value(Vbarg(x,p,t)) μst = ForwardDiff.value(x[11]) μgt =
#   ForwardDiff.value(x[12]) wmaxst = ForwardDiff.value(wmaxs(x,p,t)) wmaxgt =
#   ForwardDiff.value(wmaxgproj(x,p,t)) wst = ForwardDiff.value(x[7]) wgt =
#   ForwardDiff.value(x[8]) ηprojt = ForwardDiff.value(ηproj(x,p,t)) υbarprojst
#   = ForwardDiff.value(υbarprojs(x,p,t)) wprojst =
#   ForwardDiff.value(wprojs(x,p,t)) wprojgt = ForwardDiff.value(wprojg(x,p,t))
# 
#   #Initialize Vectors for default β_k, default β_dbar_0, counting excess, and noting whether an infrastructure has excess 
#   β_k_t = copy(p[19])
#   β_dbar_0 = 1-sum(β_k_t)
#   excess_k = zeros(length(β_k_t))
#   β_over = zeros(length(β_k_t))
#   β_max = zeros(length(β_k_t))
#   
#   ##gather excess β
#   for k in 1:length(β_k_t)
#       if(β_k_t[k] > 0) #if it is currently being targeted for expansion investments
#           #Record what would be the supply increase if using the default β
#           H_pot = β_k_t[k]*H_need_t
#           H_max = H_pot + 1 #in case H_max is not initialized, H_pot will not trigger an excess count
#           
#           #Calculate Maximum Possible Supply Increases for Each Infrastructure
#           if(k==1) #delivery efficiency
#               H_max = A_t*(p[11] - η_proj_t)
#           elseif(k==2) #storage capacity
#               H_max = μ_t*C_v_s_t*(p[12][1] - υ_bar_proj_s_t)
#           elseif(k==3) #surface processing
#               H_max = (Vbar_s_t + μ_s_t)*(w_max_s_t - w_proj_s_t)
#           elseif(k==4) #ground processing
#               H_max = (Vbar_g_t + μ_g_t)*(w_max_g_t - w_proj_g_t)
#           end
#           
#           #Note the limit to how much beta can be associated with a certain infrastructure type
#           β_max[k] = H_max/H_need
#           
#           #If Potential Excess, Note it 
#           if(H_pot > H_max)
#               excess_k[k] = (H_pot - H_max)/H_need_t
#               if(excess_k[k] > β_k_t[k]) #if H_max neg, excess > β, so all is excess 
#                   excess_k[k] = β_k_t[k]
#               end
#               β_over[k]=1
#           end
#       end
#   end
#   
#   ##redistribute excess β among available infrastructure types
#   β_k_new = zeros(length(β_k_t))
#   excess = sum(excess_k)
#   avail = β_dbar_0
#   
#   
#   #Calculate magnitude of available β
#   for k in 1:length(β_k_t)
#       if(β_over[k]==0)
#           avail = avail + β_k_t[k]
#       end
#   end
#   
#   #redistribute to available βs and take excess from over βs
#   for k in 1:length(β_k_t)
#       if(β_k_t[k] > 0) #if it is currently being targeted for expansion investments
#           if(β_over[k]==1) #if it has excesses
#               β_k_new[k] = β_k_t[k] - excess_k[k]
#           else #if it can receive excess
#               #Max Increase
#               max_incr = β_max[k] - β_k_t[k] 
#               
#               if(max_incr<0)
#                   print("WRONG")
#               end
#               
#               #Potential Increase
#               pot_incr = (β_k_t[k]/avail)*excess
#               
#               #Check if it will cause β_k to exceed β_max
#               
#               
#               
#               β_k_new[k] = β_k_t[k] + prop_excess*excess
#           end
#       end
#   end
#   
#   return β_k_new
# 
#   end;

function β_k(x,p,t)
    ##Note Supply Needs
    H_need_t = ForwardDiff.value(H_e_need(x,p,t)) ####Now β just for expansion
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
            H_pot = β_k_t[k]*H_need_t
            H_max = H_pot + 1 #in case H_max is not initialized, H_pot will not trigger an excess count
            
            #Calculate Maximum Possible Supply Increases for Each Infrastructure
            if(k==1) #delivery efficiency
                H_max = A_t*(p[11] - η_proj_t)
            elseif(k==2) #storage capacity
                H_max = μ_t*C_v_s_t*(p[12][1] - υ_bar_proj_s_t)
            elseif(k==3) #surface processing
                H_max = (Vbar_s_t + μ_s_t)*(w_max_s_t - w_proj_s_t)
            elseif(k==4) #ground processing
                H_max = (Vbar_g_t + μ_g_t)*(w_max_g_t - w_proj_g_t)
            end
            
            #Note the limit to how much beta can be associated with a certain infrastructure type
            β_max[k] = ifelse(H_max<0,0,H_max/H_need_t)
            
            #If Potential Excess, Note it 
            if(H_pot > H_max)
                excess_k[k] = (H_pot - H_max)/H_need_t
                if(excess_k[k] > β_k_t[k]) #if H_max neg, excess > β, so all is excess 
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
            print("LOTS OF LOOPS")
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

#   Needed New Investments, Beyond Maintenance (\hat{J})
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅

function J_k_new_need(x,p,t)
    ##Attention Manipulated Supply Need & Adjust for ForwardDiff Objects
    H_e_need_t = ForwardDiff.value(H_e_need(x,p,t))
    
    ##Calculate Allocation Percentages for This Year (Based on Re-distribution algorithm)
    β_k_t = β_k(x,p,t)
    
    ##Initialize J_k vector
    J_k = zeros(length(p[19])+1)
        
    ##Calculate J_k for each k
    for k in 1:length(p[19])
        #error check in model runs 
        if(β_k_t[k]<0) 
            print("β_k"*string(k)*"is neg")
        elseif(H_e_need_t<0)
            print("H_need is neg")
        end
        
        J_k[k+1] = p[22][2+k]*((β_k_t[k]*H_e_need_t)^p[23][3+k])
    end
    
    ##For d_bar (long-term demand management)
    J_k[1] = p[22][2]*((((1-sum(β_k_t))*H_e_need_t)/p[17][2])^p[23][3]) ###add consideration for gamma_L
    
    return J_k
end;        

function J_k_need(x,p,t)
    J_k_m_need_t = J_k_m_need(x,p,t)
    J_k_new_need_t = J_k_new_need(x,p,t)
    
    J_k_need_t = J_k_new_need_t .+ J_k_m_need_t
    
    return J_k_need_t
end;

#Long-Term Investment
function J_need(x,p,t)
    ##Total Investment Needs
    J_need_t = sum(J_k_need(x,p,t))
     
    return J_need_t
end;

#   Saturation (Based on Budget Potential) (sat(\hat{J}))
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅

function J_k_pot(x,p,t)
    J_k_pot = copy(J_k_need(x,p,t))
    J_bar_t = ForwardDiff.value(J_bar(x,p,t))
    J_need_t = sum(J_k_pot)
    
    if(J_need_t > J_bar_t)
        J_k_pot = J_k_pot.*(J_bar_t/J_need_t)
    end
    
    return J_k_pot
end;

function J_pot(x,p,t)
    return sum(J_k_pot(x,p,t))
end;

#   Actual New Investments (After Saturation from Actual Rate Changes)
#  (sat_2(sat_1(\hat{J}))
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅

function J(x,p,t)
    return J_pot(x,p,t)
end;

function J_o(x,p,t)
    J_t = J(x,p,t)
    J_o_max = J_bar_o(x,p,t)
    
    return min(J_o_max, J_t)
end;

function J_b(x,p,t)
    J_t = J(x,p,t)
    J_o_t = J_o(x,p,t)
    
    return max(J_t - J_o_t, 0)
end;

function J_k(x,p,t)
    return J_k_pot(x,p,t)
end;

function J_m(x,p,t)
    return min(J_m_need(x,p,t), J(x,p,t)) 
end;

function J_e(x,p,t)
    return max(J(x,p,t) - J_m(x,p,t), 0)
end;

#   1.3.3.3 Short-Term Investment Allocation
#   ------------------------------------------
# 
#   There is only 1 infrastructure option (demand management) in the short-term
# 
# :H
# 
#   ^dt = \chit Y^st \alpha \psis (1-\frac{\chi{min}}{\chit}) $

function H_d(x,p,t)
    χ_min = p[14]/μ(x,p,t)
    
    if(p[18][1]==0)
        H_d = x[2]*Y_s(x,p,t)*p[30]*(1-(χ_min/x[2]))
    else
        H_d = x[2]*Y_s(x,p,t)*p[18][1]*p[30]*(1-(χ_min/x[2]))
    end
    
    if (H_d > (x[2] - χ_min)) #if investment would push the per-capita demand below the minimum demand, lower the investment to the most allowed
        H_d = x[2] - χ_min
    end
    
    return H_d
end;

#   1.3.3.4 Infrastructure Performance Metric Functions
#   -----------------------------------------------------

#   Note Long Term Investments (J) Made in a Year ($) and Count Year as
#  Investment Year (n)
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅

J_init(x,p,t) = J(x,p,t);

function n_i(x,p,t)
    J_e_init_t = J_init(x,p,t) - J_m(x,p,t)  
    
    return ifelse((J_e_init_t/R(x,p,t))>0.005, 1, 0)
end;

#   Note Long Term Investments (H) Implemented in a Year
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅

function H_impl_k(H,x,p,t)
    return H
end;

#   1.3.3.5 Long-Term Investment Allocations
#   ------------------------------------------

function H_k(x,p,t)
    H_k = zeros(7) 
    J_k_t = J_k(x,p,t)
    η_t = ForwardDiff.value(x[15])
    
    for k in 1:7
        if(J_k_t[k]>0)
            if(k==2)
                H_k[k] = (1/η_t)*(((1/p[22][k+1])*(J_k_t[k]))^(1/p[23][k+2]))
            else
                H_k[k] = ((1/p[22][k+1])*(J_k_t[k]))^(1/p[23][k+2])
            end
        end
    end
    
    return H_k
end;

#   1.3.4 Rate-Setting Action Situation Response
#   ––––––––––––––––––––––––––––––––––––––––––––––

#   Rate Increase/Decrease (f_change)
#   -----------------------------------
# 
#   
# \Delta f_t = \frac{\hat{C}^d_t}{\bar{\pi} P_t} e_{3t} Y_{3t}
# $
# 
# *Note that f_change cannot rise above $\psi^r
# 
#   , the maximum rate increase constraint

function f_change(x,p,t)
    f_change = (1/(P(x,p,t)*p[2]))*C_d_need(x,p,t)*e_r(x,p,t)*Y_r(x,p,t)
    
    #Check for maximum rate increase limit
    if(f_change+x[17] > 1)
        f_change = 1 - x[17]
    elseif(f_change/x[17] > p[18][2])
        f_change = p[18][2]*x[17]
    end
    
    return f_change
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
#   18-24. Planned Investment (u^*): see if-else statements in the function
#   below
# 
#   25-31. Implementation Year (t^*): see if-else statements in the function
#   below
# 
#   32. Error Integral (e_i): xnew[32] = e^i_{t+1} = e^i_{t} + e^l_t OR xnew[32]
#   = e^{l,1}_{t+1} = e^{l,0}_t (for limited memory version)
# 
#   33-41. Long-Term Attention Rollover (\rho^{l,m}_{t-1}): xnew[31] =
#   \rho^{l,m}_{t+1} = \rho^l_t

#   Controller Unifying Equation (U)
#   ----------------------------------

function u(x,p,t)
    U_k = H_k(x,p,t)
    U_d = ForwardDiff.value(H_d(x,p,t))
    U_r = ForwardDiff.value(f_change(x,p,t))
    
    U = copy(U_k)
    push!(U,U_d)
    push!(U,U_r)
    
    return U
end;

#   Water System (Plant) Equations of Motion
#   ------------------------------------------

#   Implement or Store Long-Term Investment
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅

function ImplementOrStoreLTInvest(x,u,p,t)
    tf = convert(Float64, t)
    H_impl_k = zeros(length(u)-2)
    H_m_k = zeros(length(u)-2)
    
    ####Implement Maintenance
    H_m_k[2] = min(H_m_η(x,p,t), copy(u[2]))
    H_impl_k[2] = H_m_k[2]
    H_m_k[4] = min(H_m_w_s(x,p,t), copy(u[4]))
    H_impl_k[4] = H_m_k[4]
    H_m_k[5] = min(H_m_w_g(x,p,t), copy(u[5]))
    H_impl_k[5] = H_m_k[5]
    
    ############Store Un-Implemented Investments
    H_plan_k = zeros(planInvestMemorySize(p))
    
    ############Decide Whether a Potential Investment will be Implemented based on τ_i
    counter = 1 #index to keep track of place
    for k in 1:length(H_impl_k)
        if(p[16][k]==1) ####If τ_i is 1 there is no need for planning or separately accounting for maintenance 
            H_impl_k[k] = copy(u[k])
        else
            H_ready = ForwardDiff.value(x[counter+18]) #Note the investment that will be implemented in t 
            
            H_impl_k[k] = H_impl_k[k] + H_ready #implement the next investment in the queue
            
            if(p[16][k]>2)
                for c in 3:p[16][k] #move up past investments in the queue
                    H_next = ForwardDiff.value(x[counter+19])
                    H_plan_k[counter] = H_next
                    counter += 1
                end
            end
            
            H_plan_k[counter] = u[k] - H_m_k[k] #store non-maintenance investments made in year t
            counter += 1
        end
    end
    
    ############Return Implemented Investments & New Stored Investments
    return [H_impl_k, H_plan_k]
end;

#   Update Infrastructure
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅

function UpdateInfrastructure(x,H_impl,C_v_new,p,t) #H is the implemented H for that year
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
    H=zeros(7)
    for k in 1:7
        H[k]=ForwardDiff.value(H_impl[k])
    end
    
    #########Inflow
    #Surface
    if(p[25][1] == 2) #sudden change scenario
        if(p[24] == 1)
            I_new[6] = ifelse(t == p[25][3]-1, (1+p[25][2])*(I_t[6]-p[1][4]) + p[1][4], I_t[6]) + H[6]
        else
            I_new[6] = ifelse(t == p[25][3]-1, (1+p[25][2])*I_t[6], I_t[6]) + H[6]
        end
    elseif(p[25][1] == 1) #gradual change scenario
        if(t < p[25][3]-1 || t >= p[25][3]+p[25][7]-1)#if specified time for the gradual change to start has not occurred or has passed, do not change inflow
            I_new[6] = I_t[6] + H[6]
        else
            if(p[24]==1)
                b_g = (p[25][2]*p[1][5])/p[25][7]
                
                I_new[6] = p[1][4] + ((I_t[6]-p[1][4]) + b_g) + H[6]
            else
                I_new[6] = I_t[6]*(-p[25][2]) + H_impl_k[6]
            end
        end
    else #no change scenario
        I_new[6] = I_t[6] + H_impl_k[6]
    end  
    
    #Ground
    if(p[25][4] == 2) #sudden change scenario
        if(t==p[25][6]-1) #if at year of sudden change
            if(p[24]==1) #if PHX case
                if(p[25][5] == 0) #endogenous GW change based on CAP offset
                    CAP_need = (1-p[13][5])*(I_t[6]-p[1][4])*(-p[25][2])
                    I_new[7] = I_t[7]-CAP_need + H[7] #- ifelse(p[25][2]<0,p[25][10],0) #ag CAP use
                else
                    I_new[7]=I_t[7]*(-p[25][5]) + H[7]
                end
            else
                I_new[7]=I_t[7]*(-p[25][5]) + H[7]
            end
        else
            I_new[7]=I_t[7] + H[7]
        end
    elseif(p[25][4]==1) #gradual change scenario
        if(p[24]==1) #PHX case
            if(p[25][5]==0) #endogenous GW change based on CAP effect
                if(t < p[25][3]-1 || t >= p[25][3]+p[25][7]-1)
                    I_new[7] = I_t[7] + H[7]
                else
                    CAP_need = (1-p[13][5])*(-b_g)
                    I_new[7] = I_t[7]-CAP_need + H[7] - (p[25][10]/p[25][7]) #gradual ag transition 
                end
                    
            else
                I_new[7]=I_t[7]*(-p[25][5]) + H[7]
            end
        else 
            I_new[7]=I_t[7]*(-p[25][5]) + H[7]
        end
    else #No Change Scenario
       I_new[7] = I_t[7] + H[7]     
    end  
    
    ########Storage Capacity
    ##Surface
    υ_max_s_2 = I_t[3]*((C_v_old[1]*I_t[6])/(C_v_new[1]*I_new[6]))
    H_vbar_t2 = H[3]/(C_v_new[1]*I_new[6])
    I_new[3] = υ_max_s_2*(1-p[15][4]) + H_vbar_t2;
    
    ##Ground
    υ_max_g_2 = υ_bar_g_old*((C_v_old[2]*I_t[7])/(C_v_new[2]*I_new[7]))
    υ_bar_g_new = υ_max_g_2
    
    ########Processing Capacity
    #Surface
    w_s_2 = I_t[4]*((I_t[6]*(I_t[3]*C_v_old[1]+1))/(I_new[6]*(I_new[3]*C_v_new[1]+1)))
    H_w_s_t2 = H[4]/(I_new[6]*(I_new[3]*C_v_new[1]+1))
    I_new[4] = w_s_2*(1-p[15][5]) + H_w_s_t2; #SW pumping capacity
    
    #Ground
    w_g_2 = I_t[5]*((I_t[7]*(υ_bar_g_old*C_v_old[2]+1))/(I_new[7]*(υ_bar_g_new*C_v_new[2]+1)))
    H_w_g_t2 = H[5]/(I_new[7]*(υ_bar_g_new*C_v_new[2]+1))
    I_new[5] = w_g_2*(1-p[15][6]) + H_w_g_t2; #GW pumping capacity
    
    ########Delivery Efficiency
    η_next=I_t[2]*(1-p[15][3]) +  H[2]/A_t
    I_new[2]=ifelse(η_next>p[11],p[11],η_next)
    
    ########Long-Term Demand Management
    χbar_2 = I_t[1]*((I_t[6]+I_t[7]+p[1][3])/(I_new[6]+I_new[7]+p[1][3]))
    H_dbar_t2 = H[1]/((I_new[6]+I_new[7]+p[1][3])*P_t)
    I_new[1] = χbar_2*(1-p[15][2]) - H_dbar_t2;
    
    ########Return Outputs
    return [I_new, υ_bar_g_new]
end;

#   Water Users
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅

function updateUsers(x,I_new,p,t)
    ######Population Dynamics
    if(p[5][1]==1)
        p_2 = x[1]*(I_new[1]/x[16])
    else
        p_2=x[1]
    end
    
    m_i = ifelse(p[7][1]==0,0,p[7][2])
    m_o = ifelse(p[8][1]==0,0,p[8][2])
    p_new = p_2*(1 + p[6]*(1-p_2) + m_i*(1-Y_s(x,p,t)) - m_o*Y_s(x,p,t));
    
    ######Per-Capita Demand Dynamics
    d_ST_t = d_ST(x,p,t)
    χ_ST_t = d_ST_t/μ(x,p,t)
    χ_ST_2 = χ_ST_t*((x[11]+x[12]+p[1][3])/(I_new[6]+I_new[7]+p[1][3]))
    χbar_2 = x[16]*((x[11]+x[12]+p[1][3])/(I_new[6]+I_new[7]+p[1][3]))
    
    χ_new = χ_ST_2*(1+p[30]*(1-(χ_ST_2/χbar_2)))
    if(χ_new < p[14]/(I_new[6]+I_new[7]+p[1][3]))
        χ_new = p[14]/(I_new[6]+I_new[7]+p[1][3])
    end
    
    ##########Return Outupts
    return [p_new χ_new]
end;

#   Water Balance
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅

function runWaterBalance(x,C_v_new,I_new,υ_bar_g_new,p,t)
    #############Reservoir Volume
    #####Surface
    #Update
    υ_s_2 = x[3]*((x[5]*x[13]*x[11])/(I_new[3]*C_v_new[1]*I_new[6]))
    υ_new_s = υ_s_2 + (Q_a_s(x,p,t) - (O_s(x,p,t) + O_f(x,p,t) + Q_b(x,p,t)))/(I_new[3]*C_v_new[1]*I_new[6]);
    
    #Enforce Bounds
    if (υ_new_s>1)
        υ_new_s=1
    elseif(υ_new_s < 0)
        υ_new_s=0
    end
    
    #####Ground
    #Update
    υ_g_2 = x[4]*((x[6]*x[14]*x[12])/(υ_bar_g_new*C_v_new[2]*I_new[7]))
    υ_new_g = υ_g_2 + (Q_a_g(x,p,t) + Q_b(x,p,t) - O_g(x,p,t))/(υ_bar_g_new*C_v_new[2]*I_new[7])
    
    #Enforce Bounds
    if (υ_new_g>1)
        υ_new_g=1
    elseif(υ_new_g<0)
        υ_new_g=0
    end
    
    ##############Inflows for Next Year
    ####Surface
    if(p[24]==1)
        q_new_s=1
    else
        q_new_s = p[4][1]*(x[9]-1) + x[13]*sqrt(1-p[4][1]*p[4][1])*randn()+1;
    end
    
    ####Ground
    if(p[24]==1)
        q_new_g=1
    else
        q_new_g = p[4][2]*(x[10]-1) + x[14]*sqrt(1-p[4][2]*p[4][2])*randn()+1;
    end
    
    
    #########Return Outputs
    return [υ_new_s υ_new_g q_new_s q_new_g]
end;

#   Keep Track of Error
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅

function storeError(x,p,t)
    ###Setup Storage Vector
    e_new = copy(x[(planInvestMemorySize(p) + 19):end])
    
    ###Long-Term Integral 
    e_new[1] = e_i(x,p,t)
    
    #########Return Outputs
    return e_new
end;

#   Full Plant EOM
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅

function WaterSystem(x,u,p,t)
    #########Initialize Needed Variables
    u_k = u[1:(end-1)]
    u_d = u[end-1]
    u_r = u[end]
    xnew = copy(x)
    
    ####################Implement Investments#######################
    ############Long-Term
    H_update = ImplementOrStoreLTInvest(x,u,p,t)
    H_impl_k = H_update[1] 
    H_plan_k = H_update[2] 
    
    ######## Store Planned Investments Made During This Time Period & Their Investment Years############
    ##H_plan
    for c in 1:length(H_plan_k)
        xnew[18+c] = H_plan_k[c]
    end
    
    ############Short-Term
    ###Implement ST Investments Immediately
    H_d_t = u_d
    
    ####################Update C_v#######################
    ## C_v
    #SW
    if(p[3][1] == 0)
        xnew[13] = x[13]
    elseif(p[3][1] == 1)
        xnew[13] = x[13]*(1+p[3][2])
    elseif(tf==p[3][3])
        xnew[13] = x[13]*(1+p[3][4])
    else 
        xnew[13]= x[13]
    end
    
    #Gw
    xnew[14] = x[14]*(1+p[3][5])
    
    C_v_new = [ForwardDiff.value(xnew[13]) ForwardDiff.value(xnew[14])]
    
    ####################Update Infrastructure#######################
    Infrast_update = UpdateInfrastructure(x,H_impl_k,C_v_new,p,t)
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
    xnew[17] = x[17] + f_change(x,p,t) 
    #Average Bond Investment 
    xnew[18] = (x[18]*(p[31]-1) + J_b(x,p,t))*(1/p[31])
    
    ####################Store Error#######################
    Error_update = storeError(x,p,t)
    
    for i in 1:length(Error_update)
        xnew[planInvestMemorySize(p) + 18 + i] = Error_update[i]
    end
    
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
#                                  Full Parameter Name Model Variable Name                                                                                                           Definition                   Units Default Value Allowable Range
#   –––––––––––––––––––––––––––––––––––––––––––––––––– ––––––––––––––––––– –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– ––––––––––––––––––––––– ––––––––––––– –––––––––––––––
#                          Max Surface Streamflow Mean             \musmax                                                                       Max mean surface inflow that the city can seek          Bgal/yr OR AFY        200000      [0,\infty)
#                           Max Ground Streamflow Mean             \mugmax                                                             Max mean ground inflow (recharge) that the city can seek          Bgal/yr OR AFY             0      [0,\infty)
#                                  Max Purchased Water             \mupmax                                                                                         Max possible purchased water          Bgal/yr OR AFY             0      [0,\infty)
#                               Max Per-Capita Revenue             \pi_max                                                                        Max revenue that city can extract per citizen       Dollars/person/yr           400     [0, \infty)
#                   Streamflow Variation Scenario Type        \DeltaCvscen                                                       Variation Change Type (0 = no change, 1 = gradual, 2 = sudden)                Unitless             0       \{0,1,2\}
#       Surface Streamflow Percent Change in Variation        \DeltaCvs_pc                                                                   Annual proportional change in streamflow variation                Unitless         0.015    [-1, \infty)
#   Surface Streamflow Time of Sudden Variation Change      \DeltaCvssuddt                                                  Time in model run when sudden change in streamflow variation occurs                      yr            15     [1, \infty)
#            Ground Inflow Percent Change in Variation        \DeltaCvg_pc                                                                Annual proportional change in ground inflow variation                Unitless         0.015    [-1, \infty)
#                  Surface Streamflow Auto-correlation              \rho_s                                                                 1-year lagged auto-correlation in surface streamflow                Unitless           0.2           [0,1]
#                       Ground Inflow Auto-correlation              \rho_g                                                                      1-year lagged auto-correlation in ground inflow                Unitless             0           [0,1]
#                         Endogenous Carrying Capacity             K_endog                                                              Whether the model is using endogenous carrying capacity                 Boolean             1         \{0,1\}
#                               Base Carrying Capacity              K_base                                                                                  If exogenous, set carrying capacity                Unitless       1500000     (0, \infty)
#                                Intrinsic Growth Rate                   r                                                                logistically fit, intrinsic growth rate of population                Unitless           0.1      (0,\infty)
#                               Endogenous Immigration             miendog                                                                    Whether the model is using endogenous immigration                 Boolean             0         \{0,1\}
#                                     Immigration Rate              mirate                                                                  Immigration rate (proportion of current population)                Unitless          0.01     [0, \infty)
#                                Endogenous Emigration             moendog                                                                     Whether the model is using endogenous emigration                 Boolean             0         \{0,1\}
#                                      Emigration Rate              morate                                                                   Emigration rate (proportion of current population)                Unitless          0.02     [0, \infty)
#                         Surface Hedging Policy Ratio                 h_s                           Proportion of existing storage volume + inflow and demand that triggers use cuts (surface)                Unitless             1     [0, \infty)
#                          Ground Hedging Policy Ratio                 h_g                            Proportion of existing storage volume + inflow and demand that triggers use cuts (ground)                Unitless             1     [0, \infty)
#                             Min Surface Storage Fill        \upsilonmins                                                                  Minimum fill proportion of surface storage capacity                Unitless             0          [0, 1]
#                              Min Ground Storage Fill        \upsilonming                                                                   Minimum fill proportion of ground storage capacity                Unitless             0          [0, 1]
#                              Max Delivery Efficiency            \eta_max                                                                                   Max attainable delivery efficiency                Unitless           2.0          [0, 2]
#                         Max Surface Storage Capacity  \bar{\upsilon}maxs                                        Max feasible surface storage capacity (multiple of inflow standard deviation)                Unitless            12     [0, \infty)
#                          Max Ground Storage Capacity  \bar{\upsilon}maxg                                         Max feasible ground storage capacity (multiple of inflow standard deviation)                Unitless            40     [0, \infty)
#                      Max Legal Use of Ground Storage                a_gv                                                               Max legally allowed use of ground storage (proportion)                Unitless             1           [0,1]
#                       Max Legal Use of Ground Inflow                a_gq                                                                Max legally allowed use of ground inflow (proportion)                Unitless             1           [0,1]
#                     Max Legal Use of Surface Storage                a_sv                                                              Max legally allowed use of surface storage (proportion)                Unitless             1           [0,1]
#                      Max Legal Use of Surface Inflow                a_sq                                                               Max legally allowed use of surface inflow (proportion)                Unitless             1           [0,1]
#                                   Min Per Capita Use               d_min                                                                                       Min possible per-capita demand (Bgal/yr OR AFY)/person          0.04     [0, \infty)
#                                  Demand Rebound Rate            \delta_d                                                  Annual decay rate of short-term conservation measures (to baseline)                Unitless           0.5           [0,1]
#                         Background Conservation Rate         \delta_dbar                                                    Annual decay rate of long-term conservation measures (background)                Unitless        0.0116           [0,1]
#                           Surface Storage Decay Rate            \delta_v                                                                        Annual decay rate of surface storage capacity                Unitless         0.001           [0,1]
#                       Delivery Efficiency Decay Rate         \delta_\eta                                                                             Annual decay rate of delivery efficiency                Unitless         0.001           [0,1]
#                        Surface Processing Decay Rate            \deltaws                                                                     Annual decay rate of surface processing capacity                Unitless         0.001           [0,1]
#                         Ground Processing Decay Rate            \deltawg                                                                      Annual decay rate of ground processing capacity                Unitless         0.001           [0,1]
#                   LT Demand Mgmt Implementation Time              \tau_d                                                                        Time to implement long-term demand management                     yrs             3      [1,\infty)
#                        Delivery Efficiency Impl Time           \tau_\eta                                                                  Time to implement delievery efficiency improvements                     yrs             4      [1,\infty)
#                           Storage Capacity Impl Time              \tau_v                                                              Time to implement surface storage capacity improvements                     yrs             5      [1,\infty)
#                         Surface Processing Impl Time              \tauws                                                           Time to implement surface processing capacity improvements                     yrs             3      [1,\infty)
#                          Ground Processing Impl Time              \tauwg                                                            Time to implement ground processing capacity improvements                     yrs             3      [1,\infty)
#                  Surface Flow Augmentation Impl Time            \tau\mus                                                             Time to implement surface flow augmentation improvements                     yrs             4      [1,\infty)
#                   Ground Flow Augmentation Impl Time            \tau\mug                                                              Time to implement ground flow augmentation improvements                     yrs             4      [1,\infty)
#                   Short-Term Goal Supply Sufficiency            \gamma_s                                                                 Goal short-term proportion between supply and demand                Unitless             1      [0,\infty)
#                         Long-Term Goal Supply Buffer            \gamma_l                                                                  Goal long-term proportion between supply and demand                Unitless           1.1      [0,\infty)
#                  Minimum Debt Service Coverage Ratio            \gamma_r                                                              Minimum Debt Service Coverage Ratio Allowed in any Year                unitless             2      [0,\infty)
#                     Short-Term Allocation Proportion              \psi_s                                                           Proportion of investment potential allocated to short-term                Unitless           0.2           [0,1]
#                                    Max Rate Increase              \psi_r                                                                          Max possible proportional increase in rates                Unitless          0.06      [0,\infty)
#                         Delivery Efficiency Priority           \phi_\eta                                                  Proportion of long-term investments directed to delivery efficiency                Unitless           0.2           [0,1]
#                            Storage Capacity Priority              \phi_v                                             Proportion of long-term investments directed to surface storage capacity                Unitless           0.4           [0,1]
#                          Surface Processing Priority              \phiws                                          Proportion of long-term investments directed to surface processing capacity                Unitless           0.3           [0,1]
#                           Ground Processing Priority              \phiwg                                           Proportion of long-term investments directed to ground processing capacity                Unitless             0           [0,1]
#                        Surface Augmentation Priority            \phi\mus                                            Proportion of long-term investments directed to surface flow augmentation                Unitless             0           [0,1]
#                         Ground Augmentation Priority            \phi\mug                                             Proportion of long-term investments directed to ground flow augmentation                Unitless             0           [0,1]
#                                  Cheat Aquifer Limit              \phi_c                                                                  Whether city will follow legal ground pumping limit                 Boolean             1           [0,1]
#                    Short-Term Investment Sensitivity           \lambda_s                                                      Sensitivity (inst. friction component) in short-term investment                Unitless            20     [0, \infty)
#                     Long-Term Investment Sensitivity           \lambda_l                                                       Sensitivity (inst. friction component) in long-term investment                Unitless            20     [0, \infty)
#                             Rate-Setting Sensitivity           \lambda_r                                                               Sensitivity (inst. friction component) in rate-setting                Unitless            20     [0, \infty)
#                   ST Investment Activation Threshold          \epsilon_s                                             Threshold for action (inst. friction component) in short-term investment                Unitless             0     [0, \infty)
#                   LT Investment Activation Threshold          \epsilon_l                                              Threshold for action (inst. friction component) in long-term investment                Unitless             0     [0, \infty)
#                    Rate-Setting Activation Threshold          \epsilon_r                                                      Threshold for action (inst. friction component) in rate-setting                Unitless             0     [0, \infty)
#                  Operating Cost Function Coefficient                 g_o                                                               Coefficient in operating costs function (see equation)   Dollars/(persons*Vol)           2.5     [0, \infty)
#                   LT Dem Mgmt Investment Coefficient              g_dbar                                                                       Coefficient in LT dem mgmt investment function             AFY/Dollars       4.4 E-7      [0,\infty)
#                       Del Eff Investment Coefficient              g_\eta                                                               Coefficient in delivery efficiency investment funciton             AFY/Dollars        0.0033      [0,\infty)
#              Storage Capacity Investment Coefficient              g_vbar                                                                  Coefficient in storage capacity investment function              AF/Dollars         0.003     [0, \infty)
#              SW Proc Capacity Investment Coefficient                g_ws                                                       Coefficient in surface processing capacity investment function             AFY/Dollars       0.00015      [0,\infty)
#              GW Proc Capacity Investment Coefficient                g_wg                                                        Coefficient in ground processing capacity investment function             AFY/Dollars       0.00015      [0,\infty)
#                 SW Inflow Aug Investment Coefficient              g_\mus                                                       Coefficient in surface inflow augmentation investment function             AFY/Dollars        0.0001      [0,\infty)
#                 GW Inflow Aug Investment Coefficient              g_\mug                                                        Coefficient in ground inflow augmentation investment function             AFY/Dollars        0.0001      [0,\infty)
#               Operating Cost Population Scale Factor                z_op                                                                   Population scale factor in operating cost function                unitless         0.563      [0,\infty)
#                   Operating Cost Demand Scale Factor                z_od                                                                       Demand scale factor in operating cost function                unitless         0.831      [0,\infty)
#                  LT Dem Mgmt Investment Scale Factor              z_dbar                                                              Investment scale factor for long-term demand management                unitless             1      [0,\infty)
#                      Del Eff Investment Scale Factor              z_\eta                                                                      Investment scale factor for delivery efficiency                unitless          0.82      [0,\infty)
#             Storage Capacity Investment Scale Factor              z_vbar                                                                         Investment scale factor for storage capacity                unitless             1      [0,\infty)
#             SW Proc Capacity Investment Scale Factor                z_ws                                                                      Investment scale factor for processing capacity                unitless             1      [0,\infty)
#             GW Proc Capacity Investment Scale Factor                z_wg                                                                      Investment scale factor for processing capacity                unitless             1      [0,\infty)
#                SW Inflow Aug Investment Scale Factor              z_\mus                                                                   Investment scale factor for sw inflow augmentation                unitless             1      [0,\infty)
#                GW Inflow Aug Investment Scale Factor              z_\mug                                                                   Investment scale factor for gw inflow augmentation                unitless             1      [0,\infty)
#                                   Case Specification                case                                                            Indicate which case type to follow (0 = default, 1 = PHX)                Unitless             0         \{0,1\}
#                      Mean Surface Inflow Change Type      \Delta\mustype                                            Type of change to surface mean inflow (0 = none, 1 = gradual, 2 = sudden)                Unitless             0       \{0,1,2\}
#                   Mean Surface Inflow Percent Change        \Delta\muspc                                                                                Percent change to surface mean inflow                Unitless             0     [0, \infty)
#                      Mean Sufrace Inflow Change Time         \Delta\must                                          Time that sudden change occurs or time that gradual change will be complete                     yrs             0     [0, \infty)
#                       Mean Ground Inflow Change Type      \Delta\mugtype                                            Type of change to surface mean inflow (0 = none, 1 = gradual, 2 = sudden)                Unitless             0       \{0,1,2\}
#                    Mean Ground Inflow Percent Change        \Delta\mugpc                                                            Percent change to surface mean inflow (sudden or gradual)                Unitless             0     [0, \infty)
#                       Mean Ground Inflow Change Time         \Delta\mugt                             Time that sudden change occurs or time that gradual change will be complete (0 = no end)                     yrs             0     [0, \infty)
#                               Max Surface Processing               wmaxs                                     Max attainable surface processing (proportion of storage capacity + mean inflow)                Unitless             1           [0,1]
#                                Max Ground Processing               wmaxg                                      Max attainable ground processing (proportion of storage capacity + mean inflow)                Unitless           0.5           [0,1]
#                                 Base Groundwater Use            \theta_g                                                           Proportion of Demand that is usually served by groundwater                Unitless           0.0           [0,1]
#                          Forgetting in LT Investment              forget                                                            Whether the model uses forgetting in long-term investment                 Boolean             0         \{0,1\}
#                              Memory in LT Investment              memory                                                                              Years of memory in long-term investment                     yrs            10      [1,\infty)
#                          Projection in LT Investment                proj Whether the moel uses and what type (1 = backwd dif, 2 = omniscent projection) of projection in long-term investment               3 options             0       \{0,1,2\}
#                    Projection Years in LT Investment          proj_years                                                                   Number of projection years in long-term investment                     yrs            10      [1,\infty)
#                                Streamflow Assumption            Q_assume                                      Whether the projection assumes mean streamflow (0) or \mu(1-C_v) streamflow (1)                 Boolean             0         \{0,1\}
#                                Population Assumption            P_assume                                            Whether the projection assumes no water-induced migration (0) or does (1)                 Boolean             0         \{0,1\}
#     ST Dem Mgmt Investment Effectiveness Coefficient              \alpha                                                         Coefficient for effectiveness in ST conservation investments                unitless           0.5      [0,\infty)
#   Absolute Value Error Calculation for Debt Coverage          \gammarabs                                              Whether the error calculation for debt coverage will use absolute value                 Boolean             0           [0,1]
#                                            Bond Life              \tau_b                                                                                                 Life of issued bonds                   years            15      [0,\infty)
#                                   Bond Interest Rate                 i_b                                                                                        Interest Rate of Issued Bonds                unitless          0.04           [0,1]
#               Proportion of Rates from Fixed Charges             \beta_p                                               Proportion of the expected revenue to come from fixed per user charges                unitless           0.5           [0,1]
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
#   \bar{\upsilon}s0 | Surface storage capacity (multiple of inflow standard
#   deviation) | Unitless | 4 | [0, \infty) | | Ground Storage Capacity |
#   \bar{\upsilon}g0 | Ground storage capacity (multiple of inflow standard
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
#   (proportion of withdrawn delivered) | Unitless | 1 | [0, 2] | | Base
#   Per-Capita Demand | \chibar0 | Base per-capita demand, independent of ST
#   conservation (proportion of mean inflow) | Unitless | 8E-7 | [0, \infty) | |
#   Per-Capita Revenue | f0 | Proportion of per-capita revenue to max (\pimax) |
#   Unitless | 0.5 | [0,1] |

#   Setup Function (Default)
#   ––––––––––––––––––––––––––

function Default(;μ_s_max = 200000, μ_g_max = 0.001, μ_purch_max = 0, μ_other_1 = 0, μ_other_2 = 0, π_max = 400, ΔC_v_scen = 0, ΔC_v_s_pc = 0.015, ΔC_v_s_sudd_t = 15, ΔC_v_s_sudd_pc = 2.0, 
        ΔC_v_g_pc = 0.0, ρ_s = 0.2, ρ_g = 0.0, K_endog = 1, K_base = 1500000, r = 0.1, m_i_endog = 0, m_i_rate = 0.01, m_o_endog = 0, m_o_rate = 0.01, h_s = 1.0, h_g = 1.0, υ_min_s = 0.0, 
        υ_min_g = 0.0, η_max = 1.7, υ_bar_max_s = 12.0, υ_bar_max_g = 40.0, a_gv = 1.0, a_gq = 1.0, a_sv = 1.0, a_sq = 1.0, a_q2 = 1.0, a_q3 = 0, a_q4 = 0, a_q5 = 0,
        d_min = 0.04, δ_d = 0.5, δ_dbar = 0.004, δ_v = 0.07, δ_η = 0.07, δ_w_s = 0.07, δ_w_g = 0.07, τ_d = 1.0, τ_v = 1.0, τ_η = 1.0, τ_w_s = 1.0, τ_w_g = 1.0, τ_μ_s = 1.0, τ_μ_g = 1.0, γ_s = 1, 
        γ_l = 1.1, γ_r = 2, ψ_s = 0.0, ψ_r = 0.15, β_η = 0.3, β_v = 0.3, β_w_s = 0.3, β_w_g = 0.0, β_μ_s = 0.0, β_μ_g = 0.0, λ_s = 22.0, λ_l = 22.0, λ_r = 22.0, ϵ_s = 0.0, ϵ_l = 0.0, ϵ_r = 0.0, 
        g_o = 2.5, g_dbar = 5948, g_η = 5873, g_vbar = 6000, g_ws = 4657, g_wg = 4657, g_μs = 4000, g_μg = 4000, z_op = 0.518, z_od = 1.167, z_dbar = 1, z_η = 1, z_vbar = 1, z_ws = 1, z_wg = 1, 
        z_μs = 1, z_μg = 1, case = 0, Δμ_s_type = 0, Δμ_s_pc = 0.0, Δμ_s_t = 0, Δμ_g_type = 0, Δμ_g_pc = 0, Δμ_g_t = 0, Δμ_s_τ = 0, Δμ_g_τ = 0, w_max_s = 1.0, w_max_g = 0.0, θ_g = 0, θ_1 = 0, 
        forget = 1, memory = 10, proj = 0, proj_years = 10, Q_assume = 0, P_assume = 0, α = 0.5, τ_b = 15, i_b = 0.04, γ_r_abs = 0, β_c = 0, ϕ_η = 0.6, ϕ_v = 0.3, ϕ_w_s = 0.3, ϕ_w_g = 0.0, 
        ϕ_μ_s = 0.0, ϕ_μ_g = 0.0, p_0 = 0.65, χ_0 = 8.0E-7, υ_s_0 = 1.0, υ_g_0 = 1.0, υ_bar_s_0 = 4.0, υ_bar_g_0 = 1.0, w_s_0 = 1.0, w_g_0 = 0.0, q_s_0 = 1.0, q_g_0 = 0.0, μ_s_0 = 200000, 
        μ_g_0 = 0.001, C_v_s_0 = 0.1, C_v_g_0 = 0.01, η_0 = 1.0, χbar_0 = 8.0E-7, f_0 = 0.5, J_b_avg_0 = 69000000, β_p = 0.5)
    
    ##########Parameters############
    
    #1:μ
    μbar = [μ_s_max μ_g_max μ_purch_max μ_other_1 μ_other_2]
    
    #2:π_max
    
    #####Other Streamflow
    #3: ΔC_v
    ΔC_v = [ΔC_v_scen ΔC_v_s_pc ΔC_v_s_sudd_t ΔC_v_s_sudd_pc ΔC_v_g_pc]
    
    #4: q_ac
    ρ = [ρ_s, ρ_g]
    
    #24: Select Case
    
    #25: Mean Inflow Change Scneario Definition, Δμ
    Δμ = [Δμ_s_type Δμ_s_pc Δμ_s_t Δμ_g_type Δμ_g_pc Δμ_g_t Δμ_s_τ Δμ_g_τ]
        
    #####Population Growth 
        
    #5: Carrying Capacity, κ
    K_set = [K_endog K_base]
        
    #6: Intrinsic Growth rate, r
        
    #7: Immigration Rate, m_i
    m_i = [m_i_endog, m_i_rate] 
        
    #8: Emigration Rate, m_0
    m_o = [m_o_endog, m_o_rate]
        
    #####Hard Infrastructure Operations
        
    #9: Hedging Policy Coefficient, h
    h = [h_s h_g]
    
    #10: Minimum Allowed Storage Level, υ_min
    
    υ_min = [υ_min_s υ_min_g]
        
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
    δ = [δ_d δ_dbar δ_η δ_v δ_w_s δ_w_g]
        
    #16: Implementation Time, τ_i
    τ_i = [τ_d τ_η τ_v τ_w_s τ_w_g τ_μ_s τ_μ_g]
        
    #####Institutional Context
    #17: Scope Rule, Controller Goals/References
    γ = [γ_s γ_l γ_r]
    
    #18: Choice Rules/Constraints
    ψ = [ψ_s ψ_r]
        
    #19: Investment Allocation Mental Model (β)
    β = [β_η β_v β_w_s β_w_g β_μ_s β_μ_g]
        
    #20,21: Institutional Friction
    λ = [λ_s λ_l λ_r]
    ϵ = [ϵ_s ϵ_l ϵ_r]
        
    #####Costs
    #22: Cost Function Coefficients, g
    g = [g_o g_dbar g_η g_vbar g_ws g_wg g_μs g_μg]
        
    #23: Cost Function Scale Factors, z
    z = [z_op z_od z_dbar z_η z_vbar z_ws z_wg z_μs z_μg]
    
    #30: Short-term investment effectiveness coefficient, α
    
    #31: Bond LIfe, τ_b
    
    #32: Bond Interest Rate, i_b
    
    #33: Whether γ_r error is calculated with absolute value
    
    #34: Willingness to Cheat on Groundwater Allowance, β_c
    
    #35: Investment ($) Allocation Mental Model - just used for g calculation
    ϕ = [ϕ_η ϕ_v ϕ_w_s ϕ_w_g ϕ_μ_s ϕ_μ_g]
    
    #36: Proportion of the Expected Per Capita Revenue from Fixed Charges
    
    
    #####Special Options
    #28: Memory of Long-Term Integral Controller
    τ_m = [forget memory]
        
    #29: Projectioin in Long-Term Controller
    τ_p = [proj proj_years Q_assume P_assume] 
    
    #####Create Parameter Vector
    p = [μbar,π_max,ΔC_v,ρ,K_set,r,m_i,m_o,h,υ_min,η_max,υ_bar_max,a,d_min,δ,τ_i,γ,ψ,β,λ,ϵ,g,z,case,Δμ,w_max,θ,τ_m,τ_p,α,τ_b,i_b,γ_r_abs,β_c,ϕ,β_p]
    
    ##########Initial Conditions############
    e_i0 = 0.0
    
    #####Create Initial Conditions Vector
    x_0 = [p_0, χ_0,υ_s_0,υ_g_0,υ_bar_s_0,υ_bar_g_0,w_s_0,w_g_0,q_s_0,q_g_0,μ_s_0,μ_g_0,C_v_s_0,C_v_g_0,η_0,χbar_0,f_0,J_b_avg_0]
    
    ###Add Memory for Planned Investment Capacity and Time of Implementation 
    for k in 1:length(τ_i)
        if(τ_i[k]>1)
            for i in 1:τ_i[k]-1
                append!(x_0,0.0)
            end
        end
    end
    
    ###Add Memory Variables if Forgetting Setting is On (forget = 1)
    append!(x_0,e_i0) #for e_i0
    
    for i in 1:(memory*2)
        append!(x_0, 0.0)
    end
    
    return [p, x_0]
end;

#   General Parameter Auxilary Variables
#   ––––––––––––––––––––––––––––––––––––––

function planInvestMemorySize(p)
    num = 0
    for k in 1:length(p[16])
        if(p[16][k]>1)
            num += p[16][k]-1
        end
    end
    
    return floor(Int,num)
end;

#   2.2 Phoenix Municipal Area setups
#   ===================================

#   Phoenix Municipal Area Unique Parameters
#   ––––––––––––––––––––––––––––––––––––––––––

#     Full Parameter Name Model Variable Name                                 Definition    Units Default Value Allowable Range
#   ––––––––––––––––––––– ––––––––––––––––––– –––––––––––––––––––––––––––––––––––––––––– –––––––– ––––––––––––– –––––––––––––––
#         Mean SRP Inflow               \mu_1    Mean inflow into PMA through SRP canals      AFY        900000      [0,\infty)
#         Mean CAP Inflow               \mu_2    Mean inflow into PMA through CAP canals      AFY        650491      [0,\infty)
#          SRP Allocation               a_SRP Proportion of SRP inflow allocated to city Unitless         0.244           [0,1]
#          CAP Allocation               a_CAP Proportion of CAP inflow allocated to city Unitless         0.284           [0,1]
#   SRP Demand Proportion            \theta_1  Proportion of Demand that is SRP eligible Unitless           0.5           [0,1]

#   Determine Investment Cost Function Coefficients from Other Parameters
#   –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

function set_g(μ_s_0,μ_g_0,μ_tot,μ_SRP,μ_CAP,a_SRP,a_CAP,a_SRP_NCS,a_gq,a_gv,θ_1,V_g_0, A_l_g_0, V_bar_g_0,
        K_endog,K_base,p_0,χbar_0,χ_0,f_0,π_bar,g_o,z_op,z_od,J_b_avg_0,τ_b,i_b,η_0,w_g_0,w_s_0,δ_η,δ_w_s,
        δ_w_g,z_η, z_w,ϕ_η,ϕ_w,A_SRP_0)
    
    #Determine Maintenance Investment Need ($)
    a_q_0 = (a_SRP*μ_SRP + a_CAP*(μ_s_0-μ_SRP) + a_gq*μ_g_0)/μ_tot
    
    if(K_endog==1)
        K_0 = a_q_0/χbar_0
    else
        K_0 = K_base
    end
    
    P_0 = p_0*K_0
    D_0 = χ_0*P_0*μ_tot
    R_0 = f_0*π_bar*P_0
    C_o_0 = g_o*(P_0^z_op)*(D_0^z_od)
    C_d_0 = (1+τ_b*i_b)*J_b_avg_0
    J_o_0 = C_d_0
    #J_o_0 = R_0 - C_o_0 - C_d_0
    
    #Determine Capacity Need (H) Based on δ
    A_l_s_SRP = A_SRP_0
    A_l_s_CAP = μ_CAP*a_CAP
    A_l_s_0 = A_l_s_SRP + A_l_s_CAP
    A_t_s_0 = w_s_0*μ_s_0
    A_t_g_0 = w_g_0*(μ_g_0+V_bar_g_0)
    
    A_s=min(A_l_s_0,A_t_s_0)
    A_g=min(A_l_g_0,A_t_g_0)
    A_0=A_s+A_g
    
    H_m_η = (η_0)*A_0*δ_η
    H_m_ws = w_s_0*μ_s_0*δ_w_s 
    H_m_wg = w_g_0*(μ_g_0+V_bar_g_0)*δ_w_g
    
    #Determine Available Maintenance Investment (J_m)
    J_m = J_b_avg_0 + J_o_0
    
    J_m_η = ϕ_η*J_m
    J_m_w = ϕ_w*J_m
    
    #Calculate g as a function of H_m, U_m, and z
    g_η = J_m_η/((η_0*H_m_η)^z_η)
    
    g_ws = J_m_w/((H_m_ws^z_w)+(H_m_wg^z_w))
    g_wg = g_ws #treating the g's for processing capacity as the same across surface/ground
    
    return[g_η g_ws g_wg H_m_η H_m_ws H_m_wg]
end;

#   City of Phoenix Setup
#   –––––––––––––––––––––––

function Phoenix(;μ_SRP = 900000,  μ_CAP = 448663, μ_g_max = 690602, π_max = 1000, K_endog = 0, K_base = 1686528, r = 0.0878, m_i_endog = 0, m_i_rate = 0.01, m_o_endog = 0, m_o_rate = 0.01, 
        h_s = 1.0, h_g = 1.0, η_max = 1.5, a_gv = 1, a_gq = 0.05357, a_sv = 0.0, a_SRP = 0.30883, a_CAP = 0.41168, d_min = 0.072, δ_all_set = 1, δ=0.05, δ_d = 0.5, δ_dbar = 0.0003, δ_v = 0.0, 
        δ_η = 0.05, δ_w_s = 0.05, δ_w_g = 0.05, τ_all_set = 1, τ=3.0, τ_d = 1.0, τ_η = 3.0, τ_w_s = 3.0, τ_w_g = 3.0, γ_l = 1.2, γ_r = 2, ψ_s = 0.0, ψ_r = 0.15, β_η = 0.2, β_w_s = 0, 
        β_w_g = 0.7, λ_s = 110, λ_l = 22, λ_r = 22, ϵ_s = 0.0, ϵ_l = 0.0, ϵ_r = 0.0, g_o = 0.1435, g_dbar = 5948, g_vbar = 0, g_μs = 0, g_μg = 0, z_op = 0.5581, z_od = 1.0303, z_dbar = 1, 
        z_η = 1.01266, z_vbar = 1, z_w=1.04019, z_μs = 1, z_μg = 1, Δμ_s_type = 2, Δμ_s_pc = -0.284, Δμ_s_t = 14, Δμ_g_type = 0, Δμ_g_pc = 0, Δμ_g_t = 14, Δμ_s_τ = 0, Δμ_g_τ = 0, θ_g = 0.024, 
        θ_1 = 0.5, forget = 1, memory = 1, proj = 0, proj_years = 5, P_assume = 0, α = 0.5, τ_b = 15, i_b=0.04, γ_r_abs = 0, β_c = 0, ϕ_η = 0.6605, ϕ_w_s = 0.2964, ϕ_w_g = 0.0331, β_p = 0.6322,
        p_0 = 0.86466, χ_0 = 8.96274E-8, υ_bar_s_0 = 7.415E-9, w_s_0 = nothing, w_g_0 = 0.005315, C_v_s_0 = 0.001, C_v_g_0 = 0.001, η_0 = 0.9742, χbar_0 = 8.96274E-8, f_0 = 0.23836, J_b_avg_0 = 69694375, 
        LTSC=240989,CAGRD=0,QC=0,a_CAP_low=0.5324,a_CAP_high=0.3894,a_SRP_NCS=0.06367,A_SRP_0=200275.18)
    
    ####Adjust J_b_avg_0 to τ_b and ι if necessary
    C_d_0 = J_b_avg_0*(1+15*0.04) #with default settings
    J_b_avg_0 = C_d_0/(1+τ_b*i_b) 
    
    ####μ        
    μ_s_0 = μ_SRP + μ_CAP
    μ_g_0 = μ_g_max #μ, KAFY
    μ_tot = μ_s_0 + μ_g_0
    
    ####Set υ_g_0, β_g_0 and β_g_max
    V_g_0 = (100*a_gq*μ_g_0)+LTSC #100 years of safe yield allowance plus LTSC
    V_bar_g_0 = (200*a_gq*μ_g_0) #double the 100 years of safe yield alloance
    υ_bar_g_0 = V_bar_g_0/(C_v_g_0*μ_g_0)
    υ_bar_g_max = υ_bar_g_0
    υ_g_0 = V_g_0/V_bar_g_0
    
    ####W_max
    P_0 = p_0*K_base
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
    if τ_all_set == 1
        τ_η =  τ_w_s = τ_w_g = τ
    end
    
    ####Hard Infrastructure Decay Rates, δ
    if δ_all_set == 1
        δ_v = δ_η = δ_w_s = δ_w_g = δ
    end
    
    ####Investment Cost Function Normalizing Coefficients, g
    z_ws = z_w
    z_wg = z_w
    if(QC==1) #Queen Creek is calibrated to 2016 data to account for H20 aquisition
        p_0_QC=74842/K_base #2016
        χbar_0_QC = 16344.01/(μ_tot*p_0_QC*K_base) #2016
        χ_0_QC=χbar_0_QC
        f_0_QC=23734654/(p_0_QC*K_base*π_max) #2016
        J_b_avg_0_QC=3492917.5 #2016
        
        g=set_g(μ_s_0,μ_g_0,μ_tot,μ_SRP,μ_CAP,a_SRP,a_CAP,a_SRP_NCS,a_gq,a_gv,θ_1,V_g_0, A_l_g_0, V_bar_g_0, K_endog,K_base,p_0_QC,χbar_0_QC,
            χ_0_QC,f_0_QC,π_max,g_o,z_op,z_od,J_b_avg_0_QC,τ_b,i_b,η_0,w_g_0,w_s_0,δ_η,δ_w_s, δ_w_g,z_η, z_w, ϕ_η, (ϕ_w_s+ϕ_w_g),A_SRP_0)
    else
        g=set_g(μ_s_0,μ_g_0,μ_tot,μ_SRP,μ_CAP,a_SRP,a_CAP,a_SRP_NCS,a_gq,a_gv,θ_1,V_g_0, A_l_g_0, V_bar_g_0, K_endog,K_base,p_0,χbar_0,χ_0,
            f_0,π_max,g_o,z_op,z_od,J_b_avg_0,τ_b,i_b,η_0,w_g_0,w_s_0,δ_η,δ_w_s, δ_w_g,z_η, z_w, ϕ_η, (ϕ_w_s+ϕ_w_g),A_SRP_0)
    end
    g_η=g[1]
    g_ws=g[2]
    g_wg=g[3]
    
    
    return Default(μ_s_max = μ_s_0, μ_g_max = μ_g_max, μ_purch_max = 0, μ_other_1=μ_SRP, μ_other_2=μ_CAP, π_max = π_max, ρ_s = 0.2, ρ_g = 0.0, K_endog = K_endog, K_base = K_base, r = r, 
        m_i_endog = m_i_endog, m_i_rate = m_i_rate, m_o_endog = m_o_endog, m_o_rate = m_o_rate, h_s = h_s, h_g = h_g, η_max = η_max, υ_bar_max_s = 6.45E-9, υ_bar_max_g = υ_bar_g_max, a_gv = a_gv, 
        a_gq = a_gq, a_sv = a_sv, a_sq = a_SRP, a_q2 = a_CAP, a_q3=a_CAP_low, a_q4=a_CAP_high, a_q5 = a_SRP_NCS,  d_min = d_min, δ_d = δ_d, δ_dbar = δ_dbar, δ_v = δ_v, δ_η = δ_η, δ_w_s = δ_w_s, 
        δ_w_g = δ_w_g, τ_d = τ_d, τ_η = τ_η, τ_w_s = τ_w_s, τ_w_g = τ_w_g, γ_l = γ_l, γ_r = γ_r, ψ_s = ψ_s, ψ_r = ψ_r, β_η = β_η, β_v = 0, β_w_s = β_w_s, β_w_g = β_w_g, β_μ_s = 0.0, β_μ_g = 0.0, 
        λ_s = λ_s, λ_l = λ_l, λ_r = λ_r, ϵ_s = ϵ_s, ϵ_l = ϵ_l, ϵ_r = ϵ_r, g_o = g_o, g_dbar = g_dbar, g_η = g_η, g_vbar = g_vbar, g_ws = g_ws, g_wg = g_wg, g_μs = g_μs, g_μg = g_μg, z_op = z_op, 
        z_od = z_od, z_dbar = z_dbar, z_η = z_η, z_vbar = z_vbar, z_ws = z_ws, z_wg = z_wg, z_μs = z_μs, z_μg = z_μg, case = 1, Δμ_s_type = Δμ_s_type, Δμ_s_pc = Δμ_s_pc, Δμ_s_t = Δμ_s_t, 
        Δμ_g_type = Δμ_g_type, Δμ_g_pc = Δμ_g_pc, Δμ_g_t = Δμ_g_t, Δμ_s_τ = Δμ_s_τ, Δμ_g_τ = Δμ_g_τ, w_max_s = w_max_s, w_max_g = w_max_g, θ_g = θ_g, θ_1 = θ_1, forget = forget, memory = memory, 
        proj = proj, proj_years = proj_years, Q_assume = 0, P_assume = P_assume, α=α, τ_b = τ_b, i_b=i_b, γ_r_abs = γ_r_abs, ϕ_η = ϕ_η, ϕ_v = 0, ϕ_w_s = ϕ_w_s, ϕ_w_g = ϕ_w_g, ϕ_μ_s = 0, ϕ_μ_g = 0, 
        β_c = β_c, β_p=β_p, p_0 = p_0, χ_0 = χ_0, υ_s_0 = 0.0, υ_g_0 = υ_g_0, υ_bar_s_0 = υ_bar_s_0, υ_bar_g_0 = υ_bar_g_0, w_s_0 = w_s_0, w_g_0 = w_g_0, q_s_0 = 1.0, q_g_0 = 1.0, μ_s_0 = μ_s_0, 
        μ_g_0 = μ_g_0, C_v_s_0 = C_v_s_0, C_v_g_0 = C_v_g_0, η_0 = η_0, χbar_0 = χbar_0, f_0 = f_0, J_b_avg_0=J_b_avg_0) 
end;

#   City of Scottsdale Setup
#   ––––––––––––––––––––––––––

function Scottsdale(;π_max = 1000, K_endog = 0, K_base = 242300, r = 0.143, m_i_endog = 0, m_i_rate = 0.01, m_o_endog = 0, m_o_rate = 0.01, h_s = 1, h_g = 1, η_max = 1.5, a_gv = 1, a_gq = 0.03114, 
        a_sv = 0.0, a_SRP = 0.03067, a_CAP = 0.18076, d_min = 0.074, δ_all_set = 1, δ=0.05, δ_d = 0.5, δ_dbar = 0.0003, τ_all_set = 1, τ=3.0, τ_d = 1.0, τ_η = 3.0, τ_w_s = 3.0, τ_w_g = 3.0, γ_l = 1.2, 
        γ_r = 2, ψ_s = 0.0, ψ_r = 0.15, β_η = 0.2, β_w_s = 0, β_w_g = 0.7, λ_s = 110, λ_l = 22, λ_r = 22, ϵ_s = 0.0, ϵ_l = 0.0, ϵ_r = 0.0, g_o = 0.4316, z_op = 0.5581, z_od = 1.0303, 
        z_dbar = 1, z_η = 1.01266, z_vbar = 1, z_w = 1.04019, z_μs = 1, z_μg = 1, Δμ_s_type = 2, Δμ_s_pc = -0.284, Δμ_s_t = 14, Δμ_g_type = 0, Δμ_g_pc = 0, Δμ_g_t = 14, Δμ_s_τ = 0, Δμ_g_τ = 0, 
        θ_g = 0.075, θ_1 = 0.17, forget = 1, memory = 1, proj = 0, proj_years = 5, P_assume = 0, α = 0.5, τ_b = 15, i_b=0.04, γ_r_abs = 0, β_c = 0, ϕ_η = 0.6605, ϕ_w_s = 0.2409, ϕ_w_g = 0.0886,
        β_p = 0.3586, p_0 = 0.89948, χ_0 = 1.86072E-7, w_g_0 = 0.006977, w_s_0 = nothing, η_0 = 1.0562, χbar_0 = 1.86072E-7, f_0 = 0.39811, J_b_avg_0 = 12500000, LTSC = 93846, CAGRD = 0, 
        a_CAP_low=0.0472, a_CAP_high = 0.2055,a_SRP_NCS=0,A_SRP_0=13632.80)

    
    
    return Phoenix(π_max = π_max, K_endog = K_endog, K_base = K_base, r = r, m_i_endog = m_i_endog, m_i_rate = m_i_rate, m_o_endog = m_o_endog, m_o_rate = m_o_rate, h_s = h_s, h_g = h_g, 
        η_max = η_max, a_gv = a_gv, a_gq = a_gq, a_sv = a_sv, a_SRP = a_SRP, a_CAP = a_CAP, d_min = d_min, δ_all_set = δ_all_set, δ=δ,δ_d = δ_d, δ_dbar = δ_dbar, τ_all_set = τ_all_set, τ=τ, 
        τ_d = τ_d, τ_η = τ_η, τ_w_s = τ_w_s, τ_w_g = τ_w_g, γ_l = γ_l, γ_r = γ_r, ψ_s = ψ_s, ψ_r = ψ_r, β_η = β_η, β_w_s = β_w_s, β_w_g = β_w_g, λ_s = λ_s, λ_l = λ_l, λ_r = λ_r, ϵ_s = ϵ_s, ϵ_l = ϵ_l, 
        ϵ_r = ϵ_r, g_o = g_o, z_op = z_op, z_od = z_od, z_dbar = z_dbar, z_η = z_η, z_vbar = z_vbar, z_w = z_w, z_μs = z_μs, z_μg = z_μg, Δμ_s_type = Δμ_s_type, Δμ_s_pc = Δμ_s_pc, Δμ_s_t = Δμ_s_t, 
        a_CAP_low=a_CAP_low,a_CAP_high=a_CAP_high, a_SRP_NCS=a_SRP_NCS, Δμ_g_type = Δμ_g_type, Δμ_g_pc = Δμ_g_pc, Δμ_g_t = Δμ_g_t, Δμ_s_τ = Δμ_s_τ, Δμ_g_τ = Δμ_g_τ, θ_g = θ_g, θ_1 = θ_1, forget = forget, memory = memory, 
        proj = proj, proj_years = proj_years, P_assume = P_assume, α=α, τ_b = τ_b, i_b=i_b, γ_r_abs = γ_r_abs, β_c = β_c, ϕ_η = ϕ_η, ϕ_w_s = ϕ_w_s, ϕ_w_g = ϕ_w_g, β_p=β_p, p_0 = p_0, χ_0 = χ_0, 
        w_s_0 = w_s_0, w_g_0 = w_g_0, η_0 = η_0, χbar_0 = χbar_0, f_0 = f_0,J_b_avg_0=J_b_avg_0,LTSC=LTSC, CAGRD=CAGRD, QC=0, A_SRP_0=A_SRP_0)
end;

#   Town of Queen Creek Setup
#   –––––––––––––––––––––––––––

function QueenCreek(;π_max = 1000, K_endog = 0, K_base = 101553, r = 0.24, m_i_endog = 0, m_i_rate = 0.01, m_o_endog = 0, m_o_rate = 0.01, h_s = 1, h_g = 1, η_max = 1.5, a_gv = 1, a_gq = 0.02135, 
        a_sv = 0.0, a_SRP = 0.0, a_CAP = 0.01038, d_min = 0.063, δ_all_set = 1, δ=0.05, δ_d = 0.5, δ_dbar = 0.0003, τ_all_set = 1,τ=3.0, τ_d = 1.0, τ_η = 3.0, τ_w_s = 3.0, τ_w_g = 3.0, γ_l = 1.2,
        γ_r = 2, ψ_s = 0.0, ψ_r = 0.15, β_η=0.04, β_w_s=0.04, β_w_g=0.9, λ_s = 110, λ_l = 22, λ_r = 22, ϵ_s = 0.0, ϵ_l = 0.0, ϵ_r = 0.0, g_o = 0.8346, z_op = 0.5581, z_od = 1.0303, 
        z_dbar = 1, z_η = 1.01266, z_vbar = 1, z_w = 1.04019, z_μs = 1, z_μg = 1, Δμ_s_type = 2, Δμ_s_pc = -0.284, Δμ_s_t = 14, Δμ_g_type = 0, Δμ_g_pc = 0, Δμ_g_t = 14, Δμ_s_τ = 0, Δμ_g_τ = 0, 
        θ_g = 0.075, θ_1 = 0.0, forget = 1, memory = 1, proj = 0, proj_years = 5, P_assume = 0, α = 0.5, τ_b = 15, i_b=0.04, γ_r_abs = 0, β_c = 0, ϕ_η = 0.6605, ϕ_w_s = 0.0088, ϕ_w_g = 0.3207, 
        β_p = 0.5801, p_0 = 0.31705, χ_0 = 1.55365E-7, w_s_0 = 0.0005, w_g_0 =0.006737, η_0 = 0.9339, χbar_0 = 1.55365E-7, f_0 = 0.24106, J_b_avg_0 = 1322425.63, LTSC = 289535, CAGRD = 0,
        a_CAP_low=0.0594, a_CAP_high = 0.0013,a_SRP_NCS=0,A_SRP_0=0)
    
    return Phoenix(π_max = π_max, K_endog = K_endog, K_base = K_base, r = r, m_i_endog = m_i_endog, m_i_rate = m_i_rate, m_o_endog = m_o_endog, m_o_rate = m_o_rate, h_s = h_s, h_g = h_g, 
        η_max = η_max, a_gv = a_gv, a_gq = a_gq, a_sv = a_sv, a_SRP = a_SRP, a_CAP = a_CAP, d_min = d_min, δ_all_set = δ_all_set, δ=δ, δ_d = δ_d, δ_dbar = δ_dbar, τ_all_set = τ_all_set, τ=τ, 
        τ_d = τ_d, τ_η = τ_η, τ_w_s = τ_w_s, τ_w_g = τ_w_g, γ_l = γ_l, γ_r=γ_r, ψ_s = ψ_s, ψ_r = ψ_r, β_η = β_η, β_w_s = β_w_s, β_w_g = β_w_g, λ_s = λ_s, λ_l = λ_l, λ_r = λ_r, ϵ_s = ϵ_s, ϵ_l = ϵ_l, 
        ϵ_r = ϵ_r, g_o = g_o, z_op = z_op, z_od = z_od, z_dbar = z_dbar, z_η = z_η, z_vbar = z_vbar, z_w = z_w, z_μs = z_μs, z_μg = z_μg, Δμ_s_type = Δμ_s_type, Δμ_s_pc = Δμ_s_pc, Δμ_s_t = Δμ_s_t, 
        a_CAP_low=a_CAP_low,a_CAP_high=a_CAP_high, a_SRP_NCS=a_SRP_NCS, Δμ_g_type = Δμ_g_type, Δμ_g_pc = Δμ_g_pc, Δμ_g_t = Δμ_g_t, Δμ_s_τ = Δμ_s_τ, Δμ_g_τ = Δμ_g_τ, θ_g = θ_g, θ_1 = θ_1, forget = forget, memory = memory, 
        proj = proj, proj_years = proj_years, P_assume = P_assume, α=α, τ_b = τ_b, i_b=i_b, γ_r_abs = γ_r_abs, β_c = β_c, ϕ_η = ϕ_η, ϕ_w_s = ϕ_w_s, ϕ_w_g = ϕ_w_g, β_p=β_p, p_0 = p_0, χ_0 = χ_0, 
        w_s_0 = w_s_0, w_g_0 = w_g_0, η_0 = η_0, χbar_0 = χbar_0, f_0 = f_0, J_b_avg_0 = J_b_avg_0, LTSC=LTSC, CAGRD=CAGRD,QC=1, A_SRP_0=A_SRP_0)
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
# 
#     4. Phase Plots (vector of 7 plots) (see sub-section c)

#   a.) Generate Dataframe of Trajectory Variables
#   ––––––––––––––––––––––––––––––––––––––––––––––––

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
    O_p_ex = zeros(length(tr.data))
    o_p_ex = zeros(length(tr.data))
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
    Q_p_ex = zeros(length(tr.data))
    q_p_ex = zeros(length(tr.data))
    P_ex = zeros(length(tr.data))
    κbar_ex = zeros(length(tr.data))
    S_ex = zeros(length(tr.data))
    S_proj_ex = zeros(length(tr.data))
    SF_ex = zeros(length(tr.data))
    SF_proj_ex = zeros(length(tr.data))
    DSCR_ex=zeros(length(tr.data)) 
    e_s_ex = zeros(length(tr.data))
    e_l_ex = zeros(length(tr.data))
    e_i_ex = zeros(length(tr.data))
    e_r_ex = zeros(length(tr.data))
    H_d_ex = zeros(length(tr.data))
    R_ex = zeros(length(tr.data))
    R_max_ex = zeros(length(tr.data))
    Y_s_ex = zeros(length(tr.data))
    Y_l_ex = zeros(length(tr.data))
    Y_r_ex = zeros(length(tr.data))
    C_o_ex = zeros(length(tr.data))
    C_d_ex = zeros(length(tr.data))
    J_ex = zeros(length(tr.data)) 
    J_bar_ex = zeros(length(tr.data))
    κ_ex = zeros(length(tr.data))
    J_init_ex = zeros(length(tr.data)) 
    n_i_ex = zeros(length(tr.data)) 
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
    J_m_need_ws_ex = zeros(length(tr.data))
    J_m_need_η_ex = zeros(length(tr.data))
    H_ex = zeros(length(tr.data))
    H_η_ex = zeros(length(tr.data))
    H_w_s_ex = zeros(length(tr.data))
    H_w_g_ex = zeros(length(tr.data))
    H_impl_dbar_ex = zeros(length(tr.data))
    H_impl_η_ex = zeros(length(tr.data))
    H_impl_vbar_ex = zeros(length(tr.data))
    H_impl_w_s_ex = zeros(length(tr.data))
    H_impl_w_g_ex = zeros(length(tr.data))
    H_impl_μ_s_ex = zeros(length(tr.data))
    H_impl_μ_g_ex = zeros(length(tr.data))
    H_m_need_η_ex = zeros(length(tr.data))
    J_m_ex = zeros(length(tr.data))
    J_e_ex = zeros(length(tr.data))
    J_o_ex = zeros(length(tr.data))
    J_b_ex = zeros(length(tr.data))
    H_need_ex = zeros(length(tr.data))
    J_need_ex = zeros(length(tr.data))
    H_e_need_ex = zeros(length(tr.data))
    A_l_g_proj_ex=zeros(length(tr.data))
    A_g_proj_ex=zeros(length(tr.data))
    A_s_proj_ex=zeros(length(tr.data))
    A_l_s_proj_ex=zeros(length(tr.data))
    η_proj_ex = zeros(length(tr.data))
    w_s_proj_ex = zeros(length(tr.data))
    β_η_ex = zeros(length(tr.data))
    β_w_s_ex = zeros(length(tr.data))
    β_w_g_ex = zeros(length(tr.data))
    
    
    #Note Phoenix Specific Variables
    if(p[24]==1)
        O_ex_1 = zeros(length(tr.data))
        O_ex_2 = zeros(length(tr.data))
        o_ex_1 = zeros(length(tr.data))
        o_ex_2 = zeros(length(tr.data))
        Q_1_ex = zeros(length(tr.data))
        Q_2_ex = zeros(length(tr.data))
        Q_b_ex = zeros(length(tr.data))
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
        O_p_ex[i] = O_p(tr.data[i],p,t[i]) 
        o_p_ex[i] = O_p_ex[i]/μ_ex[i]
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
        Q_p_ex[i] = Q_p(tr.data[i],p,t[i])
        q_p_ex[i] = q_p(tr.data[i],p,t[i])
        Q_b_ex[i]=Q_b(tr.data[i],p,t[i])
        κbar_ex[i] = κbar(tr.data[i],p,t[i])
        P_ex[i] = P(tr.data[i],p,t[i])
        S_ex[i] = S(tr.data[i],p,t[i])
        S_proj_ex[i] = S_proj(tr.data[i],p,t[i])
        SF_ex[i] = M_sf(tr.data[i],p,t[i])
        SF_proj_ex[i] = M_sf_l(tr.data[i],p,t[i])
        DSCR_ex[i] = M_dscr(tr.data[i],p,t[i])
        e_s_ex[i] = e_s(tr.data[i],p,t[i])
        e_l_ex[i] = e_l(tr.data[i],p,t[i])
        e_i_ex[i] = e_i(tr.data[i],p,t[i])
        e_r_ex[i] = e_r(tr.data[i],p,t[i])
        H_d_ex[i] = H_d(tr.data[i],p,t[i])
        R_ex[i] = R(tr.data[i],p,t[i])
        R_max_ex[i] = R_max(tr.data[i],p,t[i])
        Y_s_ex[i] = Y_s(tr.data[i],p,t[i])
        Y_l_ex[i] = Y_l(tr.data[i],p,t[i])
        Y_r_ex[i] = Y_r(tr.data[i],p,t[i])
        C_o_ex[i] = C_o(tr.data[i],p,t[i])
        C_d_ex[i] = C_d(tr.data[i],p,t[i])
        J_ex[i] = J(tr.data[i],p,t[i])
        J_bar_ex[i] = J_bar(tr.data[i],p,t[i])
        J_init_ex[i] = J_init(tr.data[i],p,t[i])     
        n_i_ex[i] = n_i(tr.data[i],p,t[i])
        κ_ex[i] = κ(tr.data[i],p,t[i])
        A_s_ex[i] = A_s(tr.data[i],p,t[i])
        A_g_ex[i] = A_g(tr.data[i],p,t[i])
        A_s_proj_ex[i]=A_proj_s(tr.data[i],p,t[i])
        A_g_proj_ex[i]=A_proj_g(tr.data[i],p,t[i])
        A_l_s_ex[i] = A_l_s(tr.data[i],p,t[i])
        A_l_s_proj_ex[i] = A_proj_l_s(tr.data[i],p,t[i])
        A_l_g_proj_ex[i] = A_proj_l_g(tr.data[i],p,t[i])
        A_l_g_ex[i] = A_l_g(tr.data[i],p,t[i])
        A_w_s_ex[i] = A_t_s(tr.data[i],p,t[i])
        A_w_g_ex[i] = A_t_g(tr.data[i],p,t[i])
        A_ex[i] = A(tr.data[i],p,t[i])
        w_max_s_ex[i] = w_max_s(tr.data[i],p,t[i])
        w_max_g_ex[i] = w_max_g(tr.data[i],p,t[i])
        J_m_need_ex[i] = J_m_need(tr.data[i],p,t[i]) 
        J_m_need_ws_ex[i] = J_m_w_s(tr.data[i],p,t[i]) 
        J_m_need_η_ex[i] =J_m_η(tr.data[i],p,t[i]) 
        J_m_ex[i] = J_m(tr.data[i],p,t[i])
        J_e_ex[i] = J_e(tr.data[i],p,t[i])
        J_b_ex[i] = J_b(tr.data[i],p,t[i])
        J_o_ex[i] = J_o(tr.data[i],p,t[i])
        H_need_ex[i] = H_need(tr.data[i],p,t[i])
        H_e_need_ex[i] = H_e_need(tr.data[i],p,t[i])
        J_need_ex[i] = J_need(tr.data[i],p,t[i])
        H_m_need_η_ex[i] = H_m_η(tr.data[i],p,t[i])
        η_proj_ex[i] = η_proj(tr.data[i],p,t[i])
        w_s_proj_ex[i] = w_proj_s(tr.data[i],p,t[i])
        β_k_t = β_k(tr.data[i],p,t[i])
        β_η_ex[i] = β_k_t[1]
        β_w_s_ex[i] = β_k_t[3]
        β_w_g_ex[i] = β_k_t[4]
        
        #Record H 
        H_k_ex = H_k(tr.data[i],p,t[i])
        H_w_s_ex[i] = H_k_ex[4]
        H_w_g_ex[i] = H_k_ex[5]
        H_η_ex[i]=H_k_ex[2]
        H_ex[i] = sum(H_k_ex)
        
        u_t = u(tr.data[i],p,t[i])
        H_impl_k_ex = ImplementOrStoreLTInvest(tr.data[i],u_t,p,t[i])[1]
        H_impl_dbar_ex[i] = H_impl_k_ex[1]
        H_impl_η_ex[i] = H_impl_k_ex[2]
        H_impl_vbar_ex[i] = H_impl_k_ex[3]
        H_impl_w_s_ex[i] = H_impl_k_ex[4]
        H_impl_w_g_ex[i] = H_impl_k_ex[5]
        H_impl_μ_s_ex[i] = H_impl_k_ex[6]
        H_impl_μ_g_ex[i] = H_impl_k_ex[7]
        
        if(p[24]==1)
            O_ex_1[i] = O_1(tr.data[i],p,t[i])
            O_ex_2[i] = O_2(tr.data[i],p,t[i])
            o_ex_1[i] = O_ex_1[i]/(μ_s[i]+μ_g[i])
            o_ex_2[i] = O_ex_2[i]/(μ_s[i]+μ_g[i])
            Q_1_ex[i] = p[1][4]*p[13][4]
            
            Q_CAP = (Q_s(tr.data[i],p,t[i]) - p[1][4])
            CAP_short = p[1][5] - Q_CAP
            NIA_avail = max(0,70022-CAP_short)
            high_avail = Q_CAP - NIA_avail
            Q_2_ex[i] = p[13][6]*NIA_avail + p[13][7]*high_avail
        end
    end
    
    ##Shortage Calculations
    ω_pre = ifelse.(D_bar_ex.> S_ex,(D_bar_ex-S_ex)./D_bar_ex,0)
    ω_post = ifelse.(D_ST_ex.> S_ex,(D_ST_ex-S_ex)./D_ST_ex,0)
    n_ω_pre = ifelse.(ω_pre .> 0, 1, 0)
    n_ω_post = ifelse.(ω_post .> 0, 1, 0)
    
    ##Aggregate Variables into Single Dataframe
    if(p[24]==1)
        vars = DataFrame(t=t,year=year_ex,p=p_ex, χ=χ, υ_s=υ_s, υ_g=υ_g, υ_bar_s=υ_bar_s, υ_bar_g=υ_bar_g, w_s=w_s, w_g=w_g, q_s=q_s_ex, q_g=q_g_ex, Q=Q_ex, q=q_ex, Q_s=Q_s_ex, Q_g=Q_g_ex, Q_p=Q_p_ex, q_p = q_p_ex, μ_s=μ_s, μ_g=μ_g, 
            C_v_s=C_v_s, C_v_g=C_v_g, η=η, χbar=χbar, f=f,  μ=μ_ex, V_s = V_s_ex, V_g = V_g_ex, Vbar_s = Vbar_s_ex, Vbar_g = Vbar_g_ex, D=D_ex, D_proj=D_proj_ex,υ_proj_g=υ_proj_g_ex, O=O_ex, o=o_ex, O_d = O_d_ex, o_d = o_d_ex, 
            O_s=O_s_ex, o_s=o_s_ex, O_g=O_g_ex, o_g=o_g_ex, O_p=O_p_ex, o_p=o_p_ex, O_f=O_f_ex, o_f=o_f_ex, O_1=O_ex_1, o_1=o_ex_1, O_2=O_ex_2, o_2=o_ex_2, κbar = κbar_ex, P=P_ex, S=S_ex, S_proj=S_proj_ex,Q_a = Q_a_ex, q_a = q_a_ex, 
            Q_a_s = Q_a_s_ex, Q_a_g = Q_a_g_ex, SF=SF_ex, SF_proj=SF_proj_ex,DSCR=DSCR_ex, e_s=e_s_ex, e_l = e_l_ex, e_r = e_r_ex, e_i = e_i_ex, H_d=H_d_ex, R=R_ex, R_max = R_max_ex, Y_s=Y_s_ex, Y_l=Y_l_ex, Y_r = Y_r_ex, C_o = C_o_ex, 
            C_d = C_d_ex, J=J_ex, ω_pre=ω_pre,ω_post=ω_post,n_ω_pre=n_ω_pre,n_ω_post=n_ω_post,κ = κ_ex, Q_1=Q_1_ex, Q_2=Q_2_ex, n_i=n_i_ex, J_init=J_init_ex, A_s=A_s_ex,A_g=A_g_ex,A=A_ex, A_l_s=A_l_s_ex, A_l_g=A_l_g_ex, A_w_s=A_w_s_ex,
            A_w_g=A_w_g_ex, w_max_s = w_max_s_ex, w_max_g = w_max_g_ex, J_m_need=J_m_need_ex, J_m_need_ws= J_m_need_ws_ex,J_m_need_η=J_m_need_η_ex,H_w_s=H_w_s_ex,H_w_g=H_w_g_ex,H_impl_dbar =H_impl_dbar_ex,H_impl_η=H_impl_η_ex,
            H_impl_vbar =H_impl_vbar_ex, H_impl_w_s =H_impl_w_s_ex,H_impl_w_g =H_impl_w_g_ex,H_impl_μ_s =H_impl_μ_s_ex,H_impl_μ_g =H_impl_μ_g_ex, J_m=J_m_ex, J_e=J_e_ex, H_η=H_η_ex,J_b_avg=J_b_avg,J_bar=J_bar_ex,J_b=J_b_ex,J_o=J_o_ex, 
            D_bar=D_bar_ex,D_ST = D_ST_ex, d_ST = d_ST_ex,H=H_ex,H_need=H_need_ex,J_need=J_need_ex,H_m_need_η=H_m_need_η_ex,H_e_need=H_e_need_ex,A_l_g_proj=A_l_g_proj_ex,A_g_proj=A_g_proj_ex, A_s_proj=A_s_proj_ex,η_proj=η_proj_ex,
            A_l_s_proj=A_l_s_proj_ex, β_η=β_η_ex,β_w_s=β_w_s_ex,β_w_g=β_w_g_ex,Q_b=Q_b_ex,w_s_proj=w_s_proj_ex)
    else 
        vars = DataFrame(t=t,year=year_ex,p=p_ex, χ=χ, υ_s=υ_s, υ_g=υ_g, υ_bar_s=υ_bar_s, υ_bar_g=υ_bar_g, w_s=w_s, w_g=w_g, q_s=q_s_ex, q_g=q_g_ex, Q=Q_ex, q=q_ex, Q_s=Q_s_ex, Q_g=Q_g_ex, Q_p=Q_p_ex, q_p = q_p_ex, μ_s=μ_s, μ_g=μ_g, 
            C_v_s=C_v_s, C_v_g=C_v_g, η=η, χbar=χbar, f=f,  μ=μ_ex, V_s = V_s_ex, V_g = V_g_ex, Vbar_s = Vbar_s_ex,υ_proj_g=υ_proj_g_ex,Vbar_g = Vbar_g_ex, D=D_ex, D_proj=D_proj_ex, O=O_ex, o=o_ex, O_d = O_d_ex, o_d = o_d_ex, O_s=O_s_ex, 
            o_s=o_s_ex, O_g=O_g_ex, o_g=o_g_ex, O_p=O_p_ex, o_p=o_p_ex, O_f=O_f_ex, o_f=o_f_ex, κbar = κbar_ex, P=P_ex, S=S_ex, S_proj=S_proj_ex, Q_a = Q_a_ex, q_a = q_a_ex, Q_a_s = Q_a_s_ex, Q_a_g = Q_a_g_ex, SF=SF_ex, SF_proj=SF_proj_ex, 
            DSCR=DSCR_ex, e_s=e_s_ex, e_l = e_l_ex, e_r = e_r_ex, e_i = e_i_ex, H_d=H_d_ex, R=R_ex, R_max = R_max_ex, Y_s=Y_s_ex, Y_l=Y_l_ex, Y_r = Y_r_ex, C_o = C_o_ex, C_d = C_d_ex, J=J_ex, ω_pre=ω_pre,ω_post=ω_post,n_ω_pre=n_ω_pre,
            n_ω_post=n_ω_post, κ = κ_ex, n_i=n_i_ex, J_init=J_init_ex, A_s=A_s_ex,A_g=A_g_ex,A=A_ex,A_l_s=A_l_s_ex, A_l_g=A_l_g_ex, A_w_s=A_w_s_ex,A_w_g=A_w_g_ex,w_max_s = w_max_s_ex, w_max_g = w_max_g_ex,J_m_need=J_m_need_ex,
            J_m_need_ws= J_m_need_ws_ex,J_m_need_η=J_m_need_η_ex,H_w_s=H_w_s_ex,H_w_g=H_w_g_ex,H_impl_dbar =H_impl_dbar_ex,H_impl_η =H_impl_η_ex,H_impl_vbar =H_impl_vbar_ex, H_impl_w_s =H_impl_w_s_ex,H_impl_w_g =H_impl_w_g_ex,
            H_impl_μ_s =H_impl_μ_s_ex,H_impl_μ_g =H_impl_μ_g_ex, J_m=J_m_ex, J_e=J_e_ex,H_η=H_η_ex,J_b_avg=J_b_avg,J_bar=J_bar_ex,J_b=J_b_ex,J_o=J_o_ex, D_bar=D_bar_ex,D_ST = D_ST_ex, d_ST = d_ST_ex,H=H_ex,H_need=H_need_ex,J_need=J_need_ex,
            H_m_need_η=H_m_need_η_ex,H_e_need=H_e_need_ex,A_l_g_proj=A_l_g_proj_ex,A_g_proj=A_g_proj_ex,A_s_proj=A_s_proj_ex,η_proj=η_proj_ex,A_l_s_proj=A_l_s_proj_ex, β_η=β_η_ex,β_w_s=β_w_s_ex,β_w_g=β_w_g_ex)
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
    a_q_ex = a_q(x_0,p,1)
    
    ##Note Key Times
    if(((1/x_0[1])-1)/((x_0[2]*vars.κ[1])-1) <= 0)
        t_μ = 0
    else
        t_μ = (1/p[6])*log(((1/x_0[1])-1)/((x_0[2]*vars.κ[1])-1)) #time when the demand is projected to reach the mean streamflow
    end
    
    
    ##Plot Test Plots
    #Shortage, Demand, Supply, Population
    plt_short=plot(vars.t, [vars.S./vars.μ vars.κ vars.χbar.*vars.P vars.ω_pre vars.ω_post], labels = ["Supply" "Demand" "Demand_base" "Shortage_PreCons" "Shortage_w/Cons"], xlabel = "Year", ylabel = "Flow/μ", 
        linecolor = [:green :blue :indigo :pink :red], title = "Shortage",legend=:outerright)
    vline!([t_μ],labels="t_μ",linestyle=:dash, linecolor = :magenta)
    
    #Flows
    if(p[24]==1)
        plt_flows = plot(vars.t, [vars.q_a vars.o_d], labels=["In_avail" "Use"], ylims = (0,vars.q_a[1]*1.5), xlabel = "Year", ylabel="Flow/μ", title = "Inflows & Use", 
        legend = :outerright, linecolor = [:darkgreen :darkgreen], linestyle=[:solid :dot])
    else
        plt_flows = plot(vars.t, [vars.q_a vars.Q_a_s./vars.μ p[13][4].*(1 .+ vars.C_v_s) p[13][4].*(1 .- vars.C_v_s) vars.Q_a_g./vars.μ vars.o_d vars.o_s vars.o_g vars.o_p vars.o_f], 
            labels= ["In_avail" "In_s" "+σ_s" "-σ_s" "In_g" "Use_all" "Use_s" "Use_g" "Use_p" "Use_f"], ylims = (0,1.5*a_q_ex), xlabel = "Year", ylabel="Flow/μ", title = "Inflows & Use", legend = :outerright, 
            linecolor = [:black :green :green :green :turquoise :black :green :turquoise :blue :brown], linestyle=[:solid :solid :dash :dash :solid :dot :dot :dot :dot :dot])
    end
    vline!([t_μ],labels="t_μ",linestyle=:dash, linecolor = :magenta)
    
    #Storage Volume & Capacity
    plt_stor = plot(vars.t, [vars.υ_s.*vars.υ_bar_s vars.υ_bar_s], labels = ["Stor_Vol" "Stor_Capac"], ylims = (0,p[12][1].*1.5), xlabel = "Year", ylabel="Vol/(μ*C_v)", title = "Reservoir Storage", 
        legend = :outerright, linecolor = [:grey :black])
    vline!([t_μ],labels="t_μ",linestyle=:dash, linecolor = :magenta)
    
    #Error
    plt_e = plot(vars.t,[vars.e_s vars.e_i vars.e_r], labels = ["Short-Term" "Long-Term" "Rates"], xlabel = "Year", ylabel = "Error", title = "Error", 
        legend = :outerright, ylims = (-3,3), linecolor = [:blue :green :brown])
    vline!([t_μ],labels="t_μ",linestyle=:dash, linecolor = :magenta)
    
    #Attention
    plt_ρ = plot(vars.t,[vars.Y_s vars.Y_l vars.Y_r],labels = ["Short-Term" "Long-Term" "Rates"], xlabel="Year", ylabel = "Attention", 
        title = "Attention", legend=:outerright, linecolor = [:blue :green :brown])
    vline!([t_μ],labels="t_μ",linestyle=:dash, linecolor = :magenta)

    #Financial Flows
    plt_fin = plot(vars.t,[vars.R.*0.000001 vars.C_o.*0.000001 vars.C_d.*0.000001 vars.J_init.*0.000001./p[16][2]], labels = ["Rev" "Op Costs" "Debt Serv" "Inv Long"], 
        xlabel="Year", ylabel = "Dollars (M)", title = "Financial Flows", legend=:outerright, ylims=(0,(p[2]*vars.κ[1].*0.000001)), 
        linecolor = [:brown :pink :blue :green])
    vline!([t_μ],labels="t_μ",linestyle=:dash, linecolor = :magenta)
    
    #Rate-Setting
    plt_f = plot(vars.t,vars.f, labels = "f", xlabel="Time(years)", ylabel = "Per-Capita Rates/Max Per-Capita Rates", title = "Rates - Max Rates Ratio", legend=:outerright, 
        ylims=(0,1.2), linecolor = :brown)
    vline!([t_μ],labels="t_μ",linestyle=:dash, linecolor = :magenta)
    
    #Delivery Efficiency
    plt_η = plot(vars.t, [vars.η vars.H_impl_η./vars.A], xlabel = "Year",labels = ["η" "Invest_η"],ylabel="Del. Eff (%Outflow)", title = "Delivery Efficiency", 
        linecolor = :purple, linestyle = [:solid :dot], legend=:outerright, ylims = (0,2.2))
    hline!([p[11]], labels="η_max", linestyle=:dash, linecolor=:purple)
    vline!([t_μ],labels="t_μ",linestyle=:dash, linecolor = :magenta)
    
    #Per-Capita Demand
    plt_d = plot(vars.t, [vars.χ vars.χbar (vars.H_d./vars.μ).*(-1) (vars.H_impl_dbar./vars.μ).*(-1) p[14]./vars.μ], xlabel = "Year", ylabel="Flows/cap/μ",
        labels=["Dem_PC" "BaseDem_PC" "Conserv_ST" "Conserv_LT" "Dem_PC_min"], title="Per-Capita Demand", ylims = (0,1.1*vars.χbar[1]), 
        legend=:outerright, linecolor = [:blue :indigo :blue :indigo :indigo], linestyle = [:solid :solid :dot :dot :dash])
    vline!([t_μ],labels="t_μ",linestyle=:dash, linecolor = :magenta)
    
    #Pumping Capacity
    plt_w = plot(vars.t, [vars.w_s vars.w_g vars.H_impl_w_s./(vars.Vbar_s.+vars.Q_a_s) vars.H_impl_w_g./(vars.Vbar_s.+vars.Q_a_s) vars.w_max_s vars.w_max_g], xlabel = "Year",
        labels = ["w_s" "Invest_w_s" "w_g" "Invest_w_g" "w_s_max" "w_g_max"], ylabel="Surface Pump Capac (%Storage)", title = "Processing Capacity", 
        linecolor = [:green :turquoise :green :turquoise :green :turquoise], linestyle = [:solid :solid :dot :dot :dash :dash], legend=:outerright, ylims = (0,1.2))
    vline!([t_μ],labels="t_μ",linestyle=:dash, linecolor = :magenta)
    
    #Mean Inflow
    plt_μ = plot(vars.t, [vars.μ_s./vars.μ[1] vars.μ_g./vars.μ[1]], xlabel = "Year", labels = ["μ_s" "μ_g"], ylabel = "Mean Inflow/Total Initial Mean Inflow (μ_0)", title = "Mean River Inflow", 
        linecolor = [:green :turquoise], linestyle = :solid, legend = :outerright, ylims = (0,1.5))
    hline!([p[1][1]/vars.μ[1]], labels="μ_s_max", linestyle=:dash, linecolor=:green)
    hline!([p[1][2]/vars.μ[1]], labels="μ_g_max", linestyle=:dash, linecolor=:turquoise)
    vline!([t_μ],labels="t_μ",linestyle=:dash, linecolor = :magenta)
    
    #Population
    plt_P = plot(vars.t, [vars.P.*0.001 vars.κ.*0.001 vars.κbar.*0.001], ylabel="1000 Persons", xlabel = "Year", title = "Population", labels=["Pop" "κ" "κbar"], ylims=(0,vars.κ[1]*1.2*0.001), 
        linecolor = :orange, linestyle = [:solid :dash :dot], legend = :outerright)
    vline!([t_μ],labels="t_μ",linestyle=:dash, linecolor = :magenta)
    
    #Aggregate Plots
    plt = plot(plt_short, plt_flows, plt_stor, plt_e, plt_ρ, plt_fin, plt_f, plt_η, plt_d, plt_w, plt_P, plt_μ, size = (1800,1500), layout = (4,3))
    
    #Save in Plot List
    plt_list = [plt_short plt_flows plt_stor plt_e plt_ρ plt_fin plt_f plt_η plt_d plt_w plt_P plt_μ plt]
    
    return plt_list
end;

#   Dimensional Time Series
#   -------------------------

function timeSeriesPlot_dim(vars,p,x_0,units)
    a_q_ex = a_q(x_0,p,1)
    
    ##Note Key Times; If Phoenix scenario, note CAP shock time. If general, not time that demand reaches μ
    if(p[24]==1)
        t_CAP = p[25][3]
        year_CAP = vars.year[1] + t_CAP
    else
        if(((1/x_0[1])-1)/((x_0[2]*vars.κ[1])-1) <= 0)
            t_μ = 0
        else
            t_μ = (1/p[6])*log(((1/x_0[1])-1)/((x_0[2]*vars.κ[1])-1)) #time when the demand is projected to reach the mean streamflow
        end
        year_μ = vars.year[1] + t_μ
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
    else
        vline!([year_μ],labels="t_μ",linestyle=:dash, linecolor = :magenta)
    end
    
    #Flows
    if(p[24] == 1)
        plt_flows = plot(vars.year, [vars.Q_a.*0.001 vars.Q_1.*0.001 vars.Q_2.*0.001 vars.Q_a_g.*0.001 vars.O_d.*0.001 vars.O_1.*0.001 vars.O_2.*0.001 vars.O_g.*0.001], 
            labels=["In_all" "In_SRP" "In_CAP" "In_GW" "Use_all" "Use_SRP" "Use_CAP" "Use_GW"], ylims = (0,1.2*0.001*vars.Q_a[1]), xlabel = "Year", ylabel = "KAFY", title = "Inflows & Use", legend = :outerright, 
            linecolor = [:black :darkgreen :purple :turquoise :black :darkgreen :purple :turquoise], linestyle=[:solid :solid :solid :solid :dot :dot :dot :dot], legendtitle="Inflows & Use")
        vline!([year_CAP],labels="t_CAP",linestyle=:dash, linecolor = :magenta)
    else
        if(units == "AF")
            plt_flows = plot(vars.year, [vars.Q_a.*0.001 vars.Q_a_s.*0.001 (1 .+ vars.C_v_s).*vars.μ_s.*p[13][4].*0.001 (1 .- vars.C_v_s).*vars.μ_s.*p[13][4].*0.001 vars.Q_a_g.*0.001 (1 .+ vars.C_v_g).*vars.μ_g.*p[13][2].*0.001 (1 .- vars.C_v_g).*vars.μ_g.*p[13][2].*0.001 vars.O_d.*0.001 vars.O_s.*0.001 vars.O_g.*0.001 vars.O_p.*0.001 vars.O_f.*0.001], 
                labels=["In_all" "In_s" "+σ_s" "-σ_s" "In_g" "+σ_g" "-σ_g" "Use_all" "Use_s" "Use_g" "Use_p" "Flood"], ylims = (0,1.5*0.001*vars.μ[1]), xlabel = "Year", title = "Inflows & Use", ylabel = "KAFY", legend = :outerright, 
                linecolor = [:black :green :green :green :turquoise :turquoise :turquoise :black :green :turquoise :blue :brown], linestyle=[:solid :solid :dash :dash :solid :dash :dash :dot :dot :dot :dot :dot], legendtitle="Inflows & Use")
        else
            plt_flows = plot(vars.year, [vars.Q_a.*0.000000001 vars.Q_a_s.*0.000000001 (1 .+ vars.C_v_s).*vars.μ_s.*p[13][4].*0.000000001 (1 .- vars.C_v_s).*vars.μ_s.*p[13][4].*0.000000001 vars.Q_a_g.*0.000000001 (1 .+ vars.C_v_g).*vars.μ_g.*p[13][2].*0.000000001 (1 .- vars.C_v_g).*vars.μ_g.*p[13][2].*0.000000001 vars.O_d.*0.000000001 vars.O_s.*0.000000001 vars.O_g.*0.000000001 vars.O_p.*0.000000001 vars.O_f.*0.000000001], 
                labels=["In_all" "In_s" "+σ_s" "-σ_s" "In_g" "+σ_g" "-σ_g" "Use_all" "Use_s" "Use_g" "Use_p" "Flood"], ylims = (0,1.5*0.000000001*vars.μ[1]*a_q_ex), xlabel = "Year", title = "Inflows & Use", ylabel = "Bgal.yr", 
                legend = :outerright, linecolor = [:black :green :green :green :turquoise :turquoise :turquoise :black :green :turquoise :blue :brown], 
                linestyle=[:solid :solid :dash :dash :solid :dash :dash :dot :dot :dot :dot :dot], legendtitle="Inflows & Use")
        end
        vline!([year_μ],labels="t_μ",linestyle=:dash, linecolor = :magenta)
    end
    
    #Storage Volume & Capacity
    SY = p[13][2]*p[1][2]*100
    
    if(p[24] ==1)
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
        vline!([year_μ],labels="t_μ",linestyle=:dash, linecolor = :magenta) 
    end
    
    #Error
    plt_e = plot(vars.year,[vars.e_s vars.e_l vars.e_i vars.e_r], labels = ["Short-Term" "Long-Term,p" "Long-Term,i" "Rates"], xlabel = "Year", ylabel = "Error", title = "Error", 
        legend = :outerright, ylims = (-3,3), linecolor = [:blue :grey :green :brown], legendtitle="Error")
    if(p[24]==1)
        vline!([year_CAP],labels="t_CAP",linestyle=:dash, linecolor = :magenta)
    else
        vline!([year_μ],labels="t_μ",linestyle=:dash, linecolor = :magenta)
    end
    
    
    #Attention
    plt_Y = plot(vars.year,[vars.Y_s vars.Y_l vars.Y_r],labels = ["Short-Term" "Long-Term" "Rates"], xlabel="Year", ylabel = "Attention", 
        title = "Attention", legend=:outerright, linecolor = [:blue :green :brown], legendtitle="Attention")
    if(p[24]==1)
        vline!([year_CAP],labels="t_CAP",linestyle=:dash, linecolor = :magenta)
    else
        vline!([year_μ],labels="t_μ",linestyle=:dash, linecolor = :magenta)
    end

    #Financial Flows
    plt_fin = plot(vars.year,[vars.R.*0.000001 vars.C_o.*0.000001 vars.C_d.*0.000001 vars.J_bar.*0.000001 vars.J.*0.000001 vars.J_m_need.*0.000001], 
        labels = ["Rev" "Op Costs" "Debt Serv" "MaxInvest" "Inv Long" "Inv_Maint"], 
        xlabel="Year", ylabel = "Dollars (M)", title = "Financial Flows", legend=:outerright, ylims=(0,(maximum(vars.R_max)*0.000001*2)), legendtitle="Financial Flows", 
        linecolor = [:purple :red :orange :black :green :brown])
    if(p[24]==1)
        vline!([year_CAP],labels="t_CAP",linestyle=:dash, linecolor = :magenta)
    else
        vline!([year_μ],labels="t_μ",linestyle=:dash, linecolor = :magenta)
    end
    
    #Rate-Setting
    plt_f = plot(vars.year,vars.f.*p[2], labels = "f", xlabel="Time(years)", ylabel = "Per-Capita Rates (Dollars/yr)", title = "Rates - Max Rates Ratio", legend=:outerright, 
        ylims=(0,1.2*maximum(vars.f.*p[2])), linecolor = :brown, legendtitle="Rates")
    hline!([p[2]], labels = "Rate_PC_max", linestyle=:dash, linecolor=:brown) 
    if(p[24]==1)
        vline!([year_CAP],labels="t_CAP",linestyle=:dash, linecolor = :magenta)
    else
        vline!([year_μ],labels="t_μ",linestyle=:dash, linecolor = :magenta)
    end
    
    #Delivery Efficiency
    plt_η = plot(vars.year, [vars.η vars.H_impl_η./vars.A], xlabel = "Year",labels = ["η" "Invest_η"],ylabel="Del. Eff (%Outflow)", title = "Delivery Efficiency", 
        linecolor = :purple, linestyle = [:solid :dot], legend=:outerright, ylims = (0,p[11]*1.2), legendtitle="Delivery Efficiency")
    hline!([p[11]], labels="η_max", linestyle=:dash, linecolor=:purple)
    if(p[24]==1)
        vline!([year_CAP],labels="t_CAP",linestyle=:dash, linecolor = :magenta)
    else
        vline!([year_μ],labels="t_μ",linestyle=:dash, linecolor = :magenta)
    end
    
    #Per-Capita Demand
    if(units == "AF")
        plt_d = plot(vars.year, [vars.d_ST.*892.15 vars.χbar.*vars.μ.*892.15 vars.H_d.*(-1).*vars.μ.*892.15 vars.H_impl_dbar.*vars.μ.*(-1).*892.15], xlabel = "Year", ylabel="PC Dem (GPCD)",
            labels=["Dem_PC" "BaseDem_PC" "Conserv_ST" "Conserv_LT"], title="Per-Capita Demand", ylims = (0,1.1*vars.χbar[1]*892.15*vars.μ[1]), legend=:outerright, linecolor = [:blue :indigo :blue :indigo], 
            linestyle = [:solid :solid :dot :dot], legendtitle="Per-Capita Demand")
        hline!([p[14]*892.15], labels="Dem_PC_min", linestyle=:dash, linecolor=:indigo)
    else
       plt_d = plot(vars.year, [vars.d_ST./365 vars.χbar.*vars.μ./365 vars.H_d.*(-1)./365 vars.H_impl_plan_dbar.*(-1)./365], xlabel = "Year", ylabel="PC Dem (GPCD)", 
            labels=["Dem_PC" "BaseDem_PC" "Conserv_ST" "Conserv_LT"], title="Per-Capita Demand", ylims = (0,1.1*vars.χbar[1]/365*vars.μ[1]), legend=:outerright, linecolor = [:blue :indigo :blue :indigo], 
            linestyle = [:solid :solid :dot :dot], legendtitle="Per-Capita Demand")
        hline!([p[14]/365], labels="Dem_PC_min", linestyle=:dash, linecolor=:indigo)
    end
    if(p[24]==1)
        vline!([year_CAP],labels="t_CAP",linestyle=:dash, linecolor = :magenta)
    else
        vline!([year_μ],labels="t_μ",linestyle=:dash, linecolor = :magenta)
    end
    
    #Pumping Capacity
    if(units == "AF")
        plt_w = plot(vars.year, [vars.w_s.*(vars.υ_bar_s.*vars.C_v_s .+ 1).*vars.μ_s.*0.001 vars.w_g.*(vars.υ_bar_g.*vars.C_v_g .+ 1).*vars.μ_g.*0.001 vars.H_impl_w_s.*0.001 vars.H_impl_w_g*0.001 vars.w_max_s.*(vars.υ_bar_s.*vars.C_v_s .+ 1).*vars.μ_s.*0.001 vars.w_max_g.*(vars.υ_bar_g.*vars.C_v_g .+ 1).*vars.μ_g.*0.001], 
            xlabel = "Year", ylabel = "KAF/yr", labels = ["w_s" "w_g" "Invest_w_s" "Invest_w_g" "w_s_max" "w_g_max"], title = "Processing Capacity", linecolor = [:green :turquoise :green :turquoise :green :turquoise], 
            linestyle = [:solid :solid :dot :dot :dash :dash], legend=:outerright, ylims = (0, max(maximum(vars.A_l_s*1.2*0.001),maximum(vars.A_l_g.*1.2*0.001))), legendtitle="Processing Capacity")
    else
        plt_w = plot(vars.year, [vars.w_s.*(vars.υ_bar_s.*vars.C_v_s .+ 1).*vars.μ_s.*0.000000001 vars.w_g.*(vars.υ_bar_g.*vars.C_v_g .+ 1).*vars.μ_g.*0.000000001 vars.H_impl_w_s.*0.000000001 vars.H_impl_w_g*0.000000001 vars.w_max_s.*(vars.υ_bar_s.*vars.C_v_s .+ 1).*vars.μ_s.*0.000000001 vars.w_max_g.*(vars.υ_bar_g.*vars.C_v_g .+ 1).*vars.μ_g.*0.000000001], 
            xlabel = "Year", ylabel = "Bgal/yr", labels = ["w_s" "w_g" "Invest_w_s" "Invest_w_g" "w_s_max" "w_g_max"], title = "Processing Capacity", linecolor = [:green :turquoise :green :turquoise :green :turquoise], 
            linestyle = [:solid :solid :dot :dot :dash :dash], legend=:outerright, ylims = (0, max(maximum(vars.w_max_s.*(vars.V_s .+ vars.μ_s)),maximum(vars.w_max_g.*(vars.Vbar_g .+ vars.μ_g)))*0.000000001*1.2), legendtitle="Processing Capacity")
    end
    if(p[24]==1)
        vline!([year_CAP],labels="t_CAP",linestyle=:dash, linecolor = :magenta)
    else
        vline!([year_μ],labels="t_μ",linestyle=:dash, linecolor = :magenta)
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
    else
        vline!([year_μ],labels="t_μ",linestyle=:dash, linecolor = :magenta)
    end
    
    #Population
    plt_P = plot(vars.year, [vars.P.*0.001 vars.κ.*0.001 vars.κbar.*0.001], ylabel="Persons (K)", xlabel = "Year", title = "Population", labels=["Pop" "κ" "κbar"], 
        ylims=(0,vars.κ[1]*1.2*0.001), linecolor = :orange, linestyle = [:solid :dash :dot], legend = :outerright, legendtitle="Population")
    if(p[24]==1)
        vline!([year_CAP],labels="t_CAP",linestyle=:dash, linecolor = :magenta)
    else
        vline!([year_μ],labels="t_μ",linestyle=:dash, linecolor = :magenta)
    end
    
    #Aggregate Plots
    plt = plot(plt_short, plt_flows, plt_stor, plt_e, plt_Y, plt_fin, plt_f, plt_η, plt_d, plt_w, plt_P, plt_μ, size = (1800,1500), layout = (4,3))
    
    #Save in Plot List
    plt_list = [plt_short plt_flows plt_stor plt_e plt_Y plt_fin plt_f plt_η plt_d plt_w plt_P plt_μ plt]
    
    return plt_list
end;

#   c.) Generate Phase Plots
#   ––––––––––––––––––––––––––

#   Description of Plots
#   ----------------------

#   The phase plot function returns a vector of 7 plots
# 
#     1. Population vs. Per-Capita Demand
# 
#     2. Storage Capacity vs. Per-Capita Demand
# 
#     3. Total Infrastructure vs. Population
# 
#     4. Total Infrastructure vs. Demand (normalized by mean inflow)
# 
#     5. Total Infrastructure vs. Mean Inflow
# 
#     6. Total Infrastructure vs. Coefficient of Variation
# 
#     7. All Plots

#   Institutional Friction Test
#   -----------------------------

function run_instfricttest(x_0, p, t_run, year_0; controllertype = "L", sensORactivthresh = "ϵ", f=f)
    p_orig = copy(p)
    x_0_orig = copy(x_0)
    
    if(sensORactivthresh == "ϵ")
        index_1 = 21
        var_range = collect(-1:0.1:1.0)
    else
        index_1 = 20
        var_range = collect(1:1:100)
    end
    
    var_def = p_orig[index_1]
    
    index_2 = ifelse(controllertype == "L", 2, ifelse(controllertype == "S", 1, 3))
    
    p_set = Array{Any}(nothing,length(var_range),length(p_orig))
    
    for i in 1:length(var_range)
        for j in 1:length(p_orig)
            if j == index_1
                temp = copy(var_def)
                temp[index_2] = var_range[i]
                p_set[i,j] = temp
            else
                p_set[i,j] = p_orig[j]
            end
        end
    end
    
    varsDF_IF = Array{Any}(nothing, length(var_range))
    
    for i in 1:length(var_range)
        model_i =  DiscreteDynamicalSystem(f, x_0_orig, p_set[i,:])
        varsDF_IF[i] = createVarsDF(model_i, p_set[i,:], t_run, year_0)
    end
    
    return varsDF_IF
        
end;

#   Phase Plots
#   -------------

function phaseSpacePlot(vars,x_0,p,instfricttest,controllertype = "L", sensORactivthresh="ϵ")
    denom_h = 2 + ifelse(p[26][1]==0, 0, 1) + ifelse(p[26][2]==0, 0, 1) 
    denom_μ = 1 + ifelse(p[1][2]==0, 0, 1)
    
    ## No Streamflow Infrastructure in Phoenix
    if(p[24]==1)
        denom_μ=0
        denom_h=3
    end
    
        t_CAP = trunc(Int, p[25][3])
        
        ### Pop & PC Demand
        plt_Pdp = plot(vars.p, vars.χ, title="Population-PC Demand", ylabel="Inflow Normalized PC Demand (Vol/cap)", xlabel="Capacity Noramlized Population", 
        legend=:outerright, ylims = (0,vars.χbar[1]*1.5), xlims=(0,1.2), color = :black, labels = "trajectory")
        scatter!([vars.p[1]], [vars.χ[1]], color = :green, labels = "t_0")
        scatter!([vars.p[end]], [vars.χ[end]], color = :brown, labels = "t_end")
        if(p[24] == 1)
            scatter!([vars.p[t_CAP]], [vars.χ[t_CAP]], color = :red, labels = "CAP Shock")
        end
        
        ### Storage & PC Demand
        plt_VmaxD = plot(vars.υ_bar_s, vars.χ, title = "Res Stor Capacity-PC Demand", ylabel="Inflow Normalized Per-Capita Demand", xlabel="Inflow Std. Dev. Normalized Storage Capacity", labels = "trajectory", linecolor = :black,
        legend=:outerright, xlims = (0,p[12][1]*1.5), ylims=(0,vars.χbar[1]*1.5))
        scatter!([vars.υ_bar_s[1]], [vars.χ[1]], color = :green, labels = "t_0")
        scatter!([vars.υ_bar_s[end]], [vars.χ[end]], color = :brown, labels = "t_end")
        if(p[24] == 1)
            scatter!([vars.υ_bar_s[t_CAP]], [vars.χ[t_CAP]], color = :red, labels = "CAP Shock")
        end
        
        ########Infrastructure Plots##############
        
        Infrast_h = (vars.η./p[11] .+ ifelse(p[24][1]==1,0,vars.υ_bar_s./p[12][1]) .+ ifelse.(p[26][1]==0,0,vars.w_s./vars.w_max_s) .+ ifelse.(p[26][2]==0,0,vars.w_g./vars.w_max_g))/denom_h
        Infrast_d = p[14]./(vars.χbar.*vars.μ)
        if(p[24]==1)
            Infrast_μ = ifelse(p[24]==1,0,(vars.μ_s./p[1][1] + ifelse.(p[1][2]==0, 0, vars.μ_g./p[1][2]))/denom_μ)
        end
        
        Infrast = (denom_h.*Infrast_h .+ Infrast_d .+ denom_μ.*Infrast_μ)./(1+denom_h+denom_μ)
    
        ### Infrastructure & Population 
        plt_PInf = plot([vars.P vars.κ], Infrast, labels = ["trajectory" "K"], xlabel = "Population", ylabel = "Infrastructure",legend=:outerright, ylims = (0,1), linestyle = [:solid :dash], xlims=(0,vars.κ[1]*1.5), 
            linecolor = :black, title = "Infrastructure-Population")
        scatter!([vars.P[1]], [Infrast[1]], color = :green, labels = "t_0")
        scatter!([vars.P[end]], [Infrast[end]], color = :brown, labels = "t_end")
        if(p[24] == 1)
            scatter!([vars.P[t_CAP]], [Infrast[t_CAP]], color = :red, labels = "CAP Shock")
        end
        
        ### Infrastructure & Demand
        a_q_ex = a_q(x_0,p,1)
        plt_DInf = plot((vars.D./vars.μ)./a_q_ex, Infrast, xlabel = "Available Inflow Normalized Demand", labels = "trajectory", ylabel = "Infrastructure",legend=:outerright,ylims = (0,1),xlims=(0,2.2),
            linecolor=:black, title = "Infrastructure-Demand")
        scatter!([(vars.D./vars.μ)[1]/a_q_ex], [Infrast[1]], color = :green, labels = "t_0")
        scatter!([(vars.D./vars.μ)[end]/a_q_ex], [Infrast[end]], color = :brown, labels = "t_end")
        if(p[24] == 1)
            scatter!([(vars.D./vars.μ)[t_CAP]/a_q_ex], [Infrast[t_CAP]], color = :red, labels = "CAP Shock")
        end
    
        ### Infrastructure & Mean Inflow
        plt_μInf = plot(vars.μ_s, Infrast, xlabel = "μ_s", ylabel = "Infrastructure", labels = "trajectory", linecolor = :black, legend=:outerright, ylims = (0,1),xlims=(0,p[1][1]*1.2), 
            title = "Infrastructure - Inflows")
        scatter!([vars.μ_s[1]], [Infrast[1]], color = :green, labels = "t_0")
        scatter!([vars.μ_s[end]], [Infrast[end]], color = :brown, labels = "t_end")
        if(p[24] == 1)
            scatter!([vars.μ_s[t_CAP]], [Infrast[t_CAP]], color = :red, labels = "CAP Shock")
        end
    
        ### Infrastructure & Inflow Variation
        plt_CvInf = plot(vars.C_v_s, Infrast, xlabel = "C_v_s", ylabel = "Infrastructure",label = "trajectory", linecolor = :black, legend=:outerright,ylims = (0,1),xlims=(0,1), 
            title = "Infrastructure - St Dev")
        scatter!([vars.C_v_s[1]], [Infrast[1]], color = :green, labels = "t_0")
        scatter!([vars.C_v_s[end]], [Infrast[end]], color = :brown, labels = "t_end")
        if(p[24] == 1)
            scatter!([vars.C_v_s[t_CAP]], [Infrast[t_CAP]], color = :red, labels = "CAP Shock")
        end
    
        ###Aggregate Plots
        plt_phase = plot(plt_Pdp, plt_VmaxD, plt_PInf, plt_DInf, plt_μInf, plt_CvInf, size = (1700,1000),layout=(2,3))
        
        return [plt_Pdp, plt_VmaxD, plt_PInf, plt_DInf, plt_μInf, plt_CvInf, plt_phase]
end;

#   d.) Aggregate Output
#   ––––––––––––––––––––––

function UWIIM_output(model; setup=Default(), t_run=100, year_0=2010, units="AF", instfricttest = 0, controllertype = "L", sensORactivthresh = "ϵ")
    p=setup[1]
    x_0=Float64.(setup[2])
    
    varsDF = createVarsDF(model,p,t_run,year_0)
    
    plt_ts_nd = timeSeriesPlot_nondim(varsDF,p,x_0)
    
    plt_ts_dim = timeSeriesPlot_dim(varsDF,p,x_0,units)
    
    plt_phase=phaseSpacePlot(varsDF, x_0, p, instfricttest)
    
    output = [varsDF, plt_ts_nd, plt_ts_dim, plt_phase]
    
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
# 
#     4. f (default = f): model definition. there's really no reason to
#        change this...
# 
#     5. instfricttest (default = 0): run an institutional friction test? 0
#        = don't run.

function run_UWIIM(setup::Any=Default(); t_run=100, year_0=2010, units="AF", f=f, instfricttest = 0, controllertype = "L", sensORactivthresh = "ϵ")
    model = create_UWIIM(setup;f=f)
    
    output = UWIIM_output(model; setup=setup, t_run=t_run, year_0=year_0, units=units, instfricttest = instfricttest, controllertype = controllertype, sensORactivthresh = sensORactivthresh)
    
    return output
end;

#   4. Example Ouptputs
#   ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

#   With Phoenix Setup
#   ====================

#Generate Output (based on 1 trajectory with Phoenix setup) 
#output_test_PHX = run_UWIIM(Phoenix(); t_run=50,year_0=2010,units="AF");
#output_test_PHX[3][end]

#   With Scottsdale Setup
#   =======================

#Generate Output (based on 1 trajectory with Phoenix setup)
#output_test_Sc = run_UWIIM(Scottsdale(); t_run=50,year_0=2010,units="AF")
#output_test_Sc[3][end]

#   With Queen Creek Setup
#   ========================

#Generate Output
#output_test_QC_01 = run_UWIIM(QueenCreek(Δμ_s_pc=-0.01); t_run=50,year_0=2010,units="AF")
#output_test_QC_01[3][end]

#output_test_QC = run_UWIIM(QueenCreek(Δμ_s_pc=-0); t_run=50,year_0=2010,units="AF")
#output_test_QC[3][end]