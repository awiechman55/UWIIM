#   Set Up Needed Code and Environment
#   ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

using Distributed

num_cores = parse(Int, ENV["SLURM_CPUS_PER_TASK"])
#num_cores = 16
addprocs(num_cores)

println("Number of cores: ", nprocs())
println("Number of workers: ", nworkers())

t1 = time()
@everywhere include("UWIIM_PHX_sensanalysis.jl")
elapsed_time = time() - t1
println("Time to Load Model: ", elapsed_time, " seconds")

#   Specify Parameters for Parameter Space Run
#   ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

#write_p_info(50,1)

#   run1d = true run2d = false run_3d = true;

#   1-D Runs
#   ========

for i in 1:length(runs_1d_index)
    print(p_names[runs_1d_saindex[i]]*" run time: ")
    
    @time runSA(runs_1d_index[i],runs_1d_matrix[i],runs_1d_matrixindex[i],runs_1d_saindex[i],nothing,nothing,
        nothing,nothing,["PHX"],[PHX_x_0],[PHX_p],file_preamble,num_t,num_it,false,1)
end

#   2-D Runs
#   ========

#   3-D Runs
#   ========

for i in 2:length(runs_3d_index_1)
    print(p_names[runs_3d_saindex_1[i]]*"-"*p_names[runs_3d_saindex_2[i]]*" 3D full run time: ")
    @time runSA(runs_3d_index_1[i],runs_3d_matrix_1[i],runs_3d_matrixindex_1[i],runs_3d_saindex_1[i],runs_3d_index_2[i],
      runs_3d_matrix_2[i],runs_3d_matrixindex_2[i],runs_3d_saindex_2[i],["PHX"],[PHX_x_0],[PHX_p],file_preamble,num_t,num_it,true,1)
end;

####### 3D Experiments - gradual range
#for i in 1:4#length(runs_3d_index_1)
#    print(p_names[runs_3d_saindex_1[i]]*"-"*p_names[runs_3d_saindex_2[i]]*" 3D gradual run time: ")
#    @time runSA(runs_3d_index_1[i],runs_3d_matrix_1[i],runs_3d_matrixindex_1[i],runs_3d_saindex_1[i],runs_3d_index_2[i],
#        runs_3d_matrix_2[i],runs_3d_matrixindex_2[i],runs_3d_saindex_2[i],["PHX"],[PHX_x_0],[PHX_p],file_preamble,num_t,num_it,true,3)
#end;