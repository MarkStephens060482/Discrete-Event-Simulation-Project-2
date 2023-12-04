using DataStructures
using Distributions
using StableRNGs
using Printf
using Dates
using CSV
using DataFrames

### use one global variable
const n_servers = 2

### Entity data structure for each order
mutable struct Order
    id::Int64
    arrival_time::Float64                   # time when the order arrives at the factory
    start_service_times::Array{Float64,1}   # array of times when the order starts construction at machine i
    finish_service_times::Array{Float64,1}  # array of times when the order finishes construction at machine i
    queue_wait_times::Array{Float64,1}      # array of times for the duration of waiting in the queue 
    completion_time::Float64 # time when the order is complete
end
# generate a newly arrived order (where paint_time and completion_time are unknown)
Order(id::Int64, arrival_time::Float64 ) = Order(id, arrival_time, Array{Float64,1}(undef,2), Array{Float64,1}(undef,2), Array{Float64,1}(undef,2), Inf)

### Events
abstract type Event end 

struct Arrival <: Event # order arrives
    id::Int64         # a unique event id
    time::Float64     # the time of the event 
end

mutable struct Finish <: Event # an order finishes processing at machine i
    id::Int64         # a unique event id
    time::Float64     # the time of the event
    server::Int64     # ID of the server that is finishing
end

struct Null <: Event # order arrives
    id::Int64    
end

### parameter structure
struct Parameters
    seed::Int64
    T::Float64
    mean_interarrival::Float64
    mean_machine_times::Array{Float64,1}
    stdev_machine1_time::Float64
    max_queue::Int64         # space available in each queue
    time_units::String
end

### State
mutable struct SystemState
    time::Float64                               # the system time (simulation time)
    stall_time_start::Float64                   # The starting time of the stall event
    stall_duration::Float64                     # The cumulative duraction of stalling M1.
    n_entities::Int64                           # the number of entities to have been served
    n_events::Int64                             # tracks the number of events to have occur + queued
    event_queue::PriorityQueue{Event,Float64}   # to keep track of future arravals/services
    order_queues::Array{Queue{Order},1}         # the system queues (1 is the arrival queue)
    in_service::Array{Union{Order,Nothing},1}   # the order currently in service at machine i if there is one
end
function SystemState( P::Parameters ) # create an initial (empty) state
    init_time = 0.0
    init_stall_time_start = 0.0
    init_stall_duration = 0.0
    init_n_entities = 0
    init_n_events = 0
    init_event_queue = PriorityQueue{Event,Float64}()
    init_order_queues = Array{Queue{Order},1}(undef,n_servers)
    init_order_queues[1] = Queue{Order}()
    init_order_queues[2] = Queue{Order}()
    init_in_service = Array{Union{Order,Nothing},1}(undef,n_servers)
    for i=1:n_servers
        init_in_service[i] = nothing
    end
    return SystemState( init_time,
                        init_stall_time_start,
                        init_stall_duration,
                        init_n_entities,
                        init_n_events,
                        init_event_queue,
                        init_order_queues,
                        init_in_service)
end

# setup random number generators
struct RandomNGs
    rng::StableRNGs.LehmerRNG
    interarrival_time::Function
    machine_times::Array{Function,1}
end
# constructor function to create all the pieces required
function RandomNGs( P::Parameters )
    rng = StableRNG( P.seed ) # create a new RNG with seed set to that required
    interarrival_time() = rand(rng, Exponential(P.mean_interarrival))
    #interarrival_time() = P.mean_interarrival
    machine_times = Array{Function,1}(undef,n_servers)
    machine_times[1] = () -> rand(rng, Normal(P.mean_machine_times[1],P.stdev_machine1_time))   # create this as a function to be consistent
    machine_times[2] = () -> rand(rng, Exponential(P.mean_machine_times[2]))
        
    return RandomNGs( rng, interarrival_time,  machine_times )
end

# initialisation function for the simulation
function initialise( P::Parameters )
    # construct random number generators and system state
    R = RandomNGs( P )
    system = SystemState( P )

    # add an arrival at time 0.0
    t0 = 0.0
    system.n_events += 1
    enqueue!( system.event_queue, Arrival(0,t0),t0)

    return (system, R)
end

### output functions (I am using formatted output, but it could use just println)

function write_parameters( output::IO, P::Parameters ) # function to writeout parameters
    T = typeof(P)
    for name in fieldnames(T)
        println( output, "# parameter: $name = $(getfield(P,name))" )
    end
end
write_parameters( P::Parameters ) = write_parameters( stdout, P )
function write_metadata( output::IO ) # function to writeout extra metadata
    (path, prog) = splitdir( @__FILE__ )
    println( output, "# file created by code in $(prog)" )
    t = now()
    println( output, "# file created on $(Dates.format(t, "yyyy-mm-dd at HH:MM:SS"))" )
end

  # function to write out headings state and events data
function write_state_header(event_file::IO)  
    println(event_file,"time,stall_time_start,stall_duration,event_id,event_type,timing,length_event_queue,length_queue1,length_queue2,in_service1,in_service2")
end

# produce an empty State dataframe
#"mean_interarrival,mean_machine1_time,stdev_machine1_time,mean_machine2_time,max_queue"
create_state_df() = hcat(DataFrame( [name => Float64[] for name in ["mean_interarrival","mean_machine1_time","mean_machine2_time","Time","stall_time_start","stall_duration"]]),
    DataFrame( [name => Int64[] for name in ["seed","max_queue","length_queue1","length_queue2","in_service1","in_service2"]]) )
    
function write_state( event_file::IO, system::SystemState, event::Event, timing::AbstractString; debug_level::Int=0)
    if typeof(event) <: Finish
        type_of_event = "Finish($(event.server))"
    else
        type_of_event = typeof(event)
    end
     
    @printf(event_file,
            "%12.3f,%12.3f,%12.3f,%6d,%9s,%6s,%4d,%4d,%4d,%4d,%4d\n",
            system.time,
            system.stall_time_start,
            system.stall_duration,
            event.id,
            type_of_event,
            timing,
            length(system.event_queue),
            length(system.order_queues[1] ),
            length(system.order_queues[2] ),
            system.in_service[1] ==nothing ? 0 : 1, 
            system.in_service[2] ==nothing ? 0 : 1 
            )
end

function write_state_df(state_df::DataFrame,system::SystemState,P::Parameters)
    data = [
        P.mean_interarrival,
        P.mean_machine_times[1],
        P.mean_machine_times[2],
        system.time,
        system.stall_time_start,
        system.stall_duration,
        P.seed,
        P.max_queue,
        length(system.order_queues[1] ),
        length(system.order_queues[2] ),
        system.in_service[1] ==nothing ? 0 : 1, 
        system.in_service[2] ==nothing ? 0 : 1]
    push!(state_df, data)
    return state_df
end

# function to write out entity data headings and data from entity attributes
function write_entity_header( entity_file::IO ) 
    println(entity_file,"id,arrival_time,M1_queue_wait_time,M1_start_service_time,M1_finish_service_time,M2_queue_wait_time,M2_start_service_time,M2_finish_service_time,completion_time")
end
#"mean_interarrival,mean_machine1_time,stdev_machine1_time,mean_machine2_time,max_queue"
create_entity_df() = DataFrame( [name => Float64[] for name in ["seed",
                                                                "mean_interarrival",
                                                                "mean_machine1_time",
                                                                "mean_machine2_time",
                                                                "max_queue",
                                                                "arrival_time",
                                                                "M1_queue_wait_time",
                                                                "M1_start_service_time",
                                                                "M1_finish_service_time",
                                                                "M2_queue_wait_time",
                                                                "M2_start_service_time",
                                                                "M2_finish_service_time",
                                                                "completion_time"]])

function write_entity( entity_file::IO, entity::Order; debug_level::Int=0)
    
     @printf(entity_file,
            "%6d,%12.3f,%12.3f,%12.3f,%12.3f,%12.3f,%12.3f,%12.3f,%12.3f\n",
            entity.id,
            entity.arrival_time,
            entity.queue_wait_times[1],
            entity.start_service_times[1],
            entity.finish_service_times[1],
            entity.queue_wait_times[2],
            entity.start_service_times[2],
            entity.finish_service_times[2],
            entity.completion_time            
            )
end

function write_entity_df(entity_df::DataFrame,entity::Order,P::Parameters)
    data = [P.seed,
            P.mean_interarrival,
            P.mean_machine_times[1],
            P.mean_machine_times[2],
            P.max_queue,
            entity.arrival_time,
            entity.queue_wait_times[1],
            entity.start_service_times[1],
            entity.finish_service_times[1],
            entity.queue_wait_times[2],
            entity.start_service_times[2],
            entity.finish_service_times[2],
            entity.completion_time]
    push!(entity_df, data)
    return entity_df
end


# A function to contruct the output file names
function output_file(seed::Int64, parameter::DataFrameRow,output::Bool)
    M1 = round(parameter.mean_machine1_time, digits = 2)
    M2 = round(parameter.mean_machine2_time, digits = 2)
    maxQ = parameter.max_queue
    A = round(parameter.mean_interarrival, digits = 2)
    # file directory and name; * concatenates strings.
    dir = pwd()*"/data/"*"/M1_"*string(M1)*"_M2_"*string(M2)*"_maxQ_"*string(maxQ)*"_A_"*string(A)*"/seed_"*string(seed) # directory name
    mkpath(dir)                          # this creates the directory 
    file_entities = dir*"/entities.csv"  # the name of the data file 
    file_state = dir*"/state.csv"        # the name of the data file 
    return (file_entities,file_state)
end

#Write Dataframe to CSV file
macro Name(arg)
   string(arg)
end

function output_DataFrame_to_CSV(dataframe::DataFrame, filename_path::String)
    # file directory and name; * concatenates strings.
    name = @Name(dataframe)
    CSV.write(filename_path*string(name)*".csv", dataframe)   
end

# A function to write the output data of the simulation to the files    
function write_output_file(file_entities::String,file_state::String,P::Parameters) 
    fid_entities = open(file_entities, "w") # open the file for writing
    fid_state = open(file_state, "w")       # open the file for writing
    
    # loop through the output files and write the metadata and parameters to set up the files
    for fid in [fid_entities,fid_state]
        write_metadata( fid )
        write_parameters( fid, P )
    end
    
    #write out the headers to the output files
    write_entity_header(fid_entities)
    write_state_header(fid_state)
    
    return (fid_entities,fid_state)
end

#A function to read the input parameter csv file
function read_parameters(file_parameter::String)
    parameter_set = CSV.read(file_parameter, DataFrame, comment="#" )
    return parameter_set    
end

### Update functions
function update!( system::SystemState, P::Parameters, R::RandomNGs, e::Event )
    throw( DomainError("invalid event type" ) )
end

function improved_update!( system::SystemState, P::Parameters, R::RandomNGs, e::Event )
    throw( DomainError("invalid event type" ) )
end

function move_to_server!( system::SystemState, R::RandomNGs, server::Integer )
    # move the order order from a queue into construction
    system.in_service[server] = dequeue!(system.order_queues[server])
    system.in_service[server].start_service_times[server] = system.time # start service 'now'   
    
    if server < n_servers
        #Calculate queueing time
        system.in_service[server].queue_wait_times[server] = system.time - system.in_service[server].arrival_time
    else
       #Calculate queueing time
       system.in_service[server].queue_wait_times[server] = system.time - system.in_service[server].finish_service_times[server-1]
    end
    completion_time = system.time + R.machine_times[server]() # best current guess at service time
    # create a finish event for the machine current constructing the item
    system.n_events += 1
    finish_event = Finish( system.n_events, completion_time, server )
    enqueue!( system.event_queue, finish_event, completion_time )
    return nothing
end

function update!( system::SystemState, P::Parameters, R::RandomNGs, event::Arrival )
    # create an arriving order and add it to the 1st queue
    system.n_entities += 1    # new entity will enter the system
    new_order = Order( system.n_entities, event.time )
    enqueue!(system.order_queues[1], new_order)
    
    # generate next arrival and add it to the event queue
    future_arrival = Arrival(system.n_events, system.time + R.interarrival_time())
    enqueue!(system.event_queue, future_arrival, future_arrival.time)

    # if the construction machine is available, the order goes to service
    if system.in_service[1] == nothing
        move_to_server!( system, R, 1 )
    end
    return nothing
end

function improved_update!( system::SystemState, P::Parameters, R::RandomNGs, event::Arrival )
    # create an arriving order and add it to the 1st queue
    system.n_entities += 1    # new entity will enter the system
    new_order = Order( system.n_entities, event.time )
    enqueue!(system.order_queues[1], new_order)
        
    # generate next arrival and add it to the event queue
    future_arrival = Arrival(system.n_events, system.time + R.interarrival_time())
    enqueue!(system.event_queue, future_arrival, future_arrival.time)

    # if  machine 1 is available and queue 2 is less than capacity, the order goes to service
    if system.in_service[1] == nothing && length(system.order_queues[2]) < P.max_queue
        
        move_to_server!( system, R, 1 )
    # queue 2 is at capacity, stall the order
    elseif length(system.order_queues[2]) == P.max_queue
        
        #record start of stall time if it hadn't been stated yet
        if system.stall_time_start == 0
            system.stall_time_start = system.time
        end
    end
    return nothing
end

function stall_event!( system::SystemState, event::Event )
    # defer an event until after the next event in the list
    next_event_time = peek( system.event_queue )[2]
    event.time = next_event_time + eps() # add eps() so that this event occurs just after the next event
    enqueue!(system.event_queue, event, event.time)
    return nothing
end

function update!( system::SystemState, P::Parameters, R::RandomNGs, event::Finish )
    server = event.server
    if server < n_servers && length(system.order_queues[server+1]) >= P.max_queue
        # if the server finishes, but there are too many people in the next queue,
        # then defer the event until the queue has space, i.e, the next finish event
        # but finding the next event is easy, and next finish is hard, so we stall by one
        stall_event!( system, event )
    else
        # otherwise treat this as normal finish of service
        departing_order = deepcopy( system.in_service[server] )
        system.in_service[server] = nothing
        
        if !isempty(system.order_queues[server]) # if someone is waiting, move them to service
            move_to_server!( system, R, server )
        end

        if server < n_servers
            # move the customer to the next queue
            enqueue!(system.order_queues[server+1], departing_order)
            if system.in_service[server+1] === nothing
                move_to_server!( system, R, server+1 )
            end
        else
            # or return the entity when it is leaving the system for good
            departing_order.completion_time = system.time
            return departing_order 
        end
    end
    return nothing
end

function improved_update!( system::SystemState, P::Parameters, R::RandomNGs, event::Finish )
    server = event.server
    #stall_event!( system, event )
    # lamp order will depart machine 
    departing_order = deepcopy( system.in_service[server] )
    departing_order.finish_service_times[server] = system.time
    system.in_service[server] = nothing
    
    if server < n_servers  # Machine 1 finish
        # move the lamp to the next queue, Machine 2 queue
        enqueue!(system.order_queues[server+1], departing_order)
        # if  Machine 2 is available, move to service straight away.
        if system.in_service[server+1] === nothing
            move_to_server!( system, R, server+1 )
        end
        # if lamp order is waiting in the Machine 1 queue and machine 2 queue is less than capacity, move lamp order to service of machine 1
        if length(system.order_queues[server+1]) < P.max_queue && !isempty(system.order_queues[server])
           #then move to machine 1 
           move_to_server!( system, R, server )
           # No Stall
           #Start_stall_time is reset
           system.stall_time_start = 0
           #Stall duration is zero
           system.stall_duration = 0 
        # if next lamp order is waiting in machine 1 queue and machine 2 queue is at capacity, next order is stalled and record time.
        elseif length(system.order_queues[server+1]) >= P.max_queue && !isempty(system.order_queues[server])
           #record start of stall time in system state
           system.stall_time_start = system.time
        end
    else     # Machine 2 finish
        # if Machine 2 queue is full and machine queue 1 is not empty
        if length(system.order_queues[server]) >= P.max_queue && !isempty(system.order_queues[server-1])
            # Record Stall duration
            system.stall_duration = (system.time - system.stall_time_start)
            #reset stall start time 
            system.stall_time_start = 0
            #move the next lamp order into service of machine 2.
            move_to_server!( system, R, server )
            #move the next lamp order into service of machine 1
            move_to_server!( system, R, server-1 )
        # if Machine 2 queue is not empty     
        elseif !isempty(system.order_queues[server]) 
            #move the next lamp order into service of machine 2.
            move_to_server!( system, R, server )  
        end
        # return the lamp when it is leaving the system for good
        departing_order.completion_time = system.time
        return departing_order 
    end
    return nothing
end

function run!( system::SystemState, P::Parameters, R::RandomNGs,fid_state::IO, fid_entities::IO,state_df::DataFrame,entity_df::DataFrame, output_level::Integer)
    # main simulation loop
    while system.time < P.T
        if P.seed ==1 && system.time <= 200.0
            # debug information for first few events whenb seed = 1
            println("system time : $(round(system.time,digits=2)), Q1 = $(length(system.order_queues[1])), S1 = $(system.in_service[1] ==nothing ? 0 : 1), Q2 = $(length(system.order_queues[2])), S2 = $(system.in_service[2] ==nothing ? 0 : 1)...") 
        end

        # grab the next event from the event queue
        (event, time) = dequeue_pair!(system.event_queue)
        system.time = time  # advance system time to the new arrival
        system.n_events += 1      # increase the event counter
        
        # write out event and state data before event
        if output_level>=2
            write_state( fid_state, system, event, "before")
            state_df = write_state_df(state_df,system,P)
        elseif output_level==1 
            # add state to dataframe
            state_df = write_state_df(state_df,system,P)
        end
        #resent the Stall time counter
        if system.stall_duration > 0
            system.stall_duration = 0
        end
        
        # update the system based on the next event, and spawn new events. 
        # return arrived/departed customer.
        departure = improved_update!( system, P, R, event )
         
        # write out event and state data after event for debugging
        if output_level>=2
            write_state( fid_state, system, event, "after")
        end
        
        # write out entity data if it was a departure from the system
        if departure !== nothing && output_level>=2
            write_entity( fid_entities, departure )
            entity_df = write_entity_df(entity_df,departure,P)
        elseif output_level == 1 && departure !== nothing
            # add state to dataframe
            entity_df = write_entity_df(entity_df,departure,P)
        end
    end
    return (state_df,entity_df)
end

### main worker function to run the simulation
function run_factory_simulation_2(seed::Int64, Time_limit::Float64,state_df::DataFrame,entity_df::DataFrame, parameter::DataFrameRow, output::Bool, file_state::String, file_entities::String)
    if seed == 1
        println("running simulation...parameters:")
        println("Arrival = $(round(parameter.mean_interarrival,digits = 2))")
        println("machine 1mean = $(round(parameter.mean_machine1_time,digits = 2))")
        println("machine 2mean = $(round(parameter.mean_machine2_time,digits = 2))")
        println("machine 1std = $(round(parameter.stdev_machine1_time,digits = 2))")
        println("machine 2 queuemax = $(round(parameter.max_queue,digits = 2))")
    end
    # define the parameter object
    time_units = "minutes"
    P = Parameters(seed,
                   Time_limit,
                   parameter.mean_interarrival,
                   [parameter.mean_machine1_time,parameter.mean_machine2_time],
                   parameter.stdev_machine1_time,
                   parameter.max_queue,
                   time_units) 
    
    #Setup and write the output files
    if output
        output_level = 2
        (fid_entities,fid_state) = write_output_file(file_entities::String,file_state::String, P::Parameters) 
    else
        output_level = 1
        (fid_entities,fid_state) = (IOBuffer(),IOBuffer())
    end
    
    #Initialisae and  run the actual simulation
    (system,rngs) = initialise( P ) 
    (state_df1,entity_df1) = run!( system, P, rngs, fid_state, fid_entities,state_df,entity_df,output_level)

    # to close the files
    if output
        close( fid_entities )
        close( fid_state )
    end
    return (state_df1,entity_df1)
end

    
