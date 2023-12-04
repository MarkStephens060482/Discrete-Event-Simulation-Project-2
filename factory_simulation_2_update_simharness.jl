include("factory_simulation_2_update.jl")

### Simulation Harness
#number of realisations
M = 150
rerun = true # rerun all simulations
output = false #produce output for simulations
if output == false
    rerun = true
end
Time_limit =350_000.0

# Define path for Parameter_set data file
parameter_path = pwd() # directory name
#read the csv file 
parameter_set = read_parameters(parameter_path*"/parameter_set.csv")
perf_times = zeros(M,1) # initialise a vector to hold performance run times of the simulation
#Initialise the DataFrame of simulated Data 
state_df =  create_state_df()
entity_df = create_entity_df()
#Generate M realisations of the simulation
i=0
Total_T1 = @elapsed begin
    for seed = 1:M
    # for the given parameters in the set of parameters 
        for parameter in eachrow(parameter_set)
            i+=1
            # define the output files        
            (file_entities,file_state) = output_file(seed, parameter,output)
            # Run the simulation based on the flag or missing data files for particular seeds and parameters
            if rerun || !(isfile(file_entities) && isfile(file_state))
               # run_factory_simulation_2(seed::Int64, Time_limit::Float64,state_df::DataFrame,entity_df::DataFrame, parameter::DataFrameRow, output::Bool, file_state::String, file_entities::String)
               perf_times[seed] = @elapsed begin
                    (state_df,entity_df) = run_factory_simulation_2(seed, Time_limit,state_df,entity_df, parameter, output, file_state, file_entities)
                   state_df = vcat(state_df)
                   entity_df = vcat(entity_df)
                end
            end
        end
        println("The run time for $(seed) realisations of sim time of $(Time_limit) is $(round(sum(perf_times),digits = 2)) secs")
        println("The average run time per simulation is $(round((sum(perf_times)/i),digits = 2)) secs")
    end   
    state_df[!,:total_orders] =  state_df.length_queue1 + state_df.length_queue2 + state_df.in_service1 + state_df.in_service2
end
println("The total run time for $(M) realisations of sim time of $(Time_limit) is $(round(sum(perf_times),digits = 2)) secs")
println("The average run time per simulation is $(round(mean(perf_times),digits = 2)) secs")
println("Total time to run the code is $(round(Total_T1,digits = 2)) secs")

#write Parameter and  M realisations State and Entity Data to CSV
filename_path = pwd()*"/data"
# make the path if it does not exist
mkpath(dir)
# output function
output_DataFrame_to_CSV(state_df, filename_path)
output_DataFrame_to_CSV(entity_df, filename_path)