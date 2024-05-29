using JuMP # building models
using DataStructures # using dictionaries with a default value
using HiGHS # solver for the JuMP model
using CSV # readin of CSV files
using DataFrames # data tables
using Statistics # mean function
using Plots  # generate graphs
using Plots.Measures
using StatsPlots # additional features for plots
include(joinpath(@__DIR__, "colors.jl")) # colors for the plots

### some helper functions ###
# read the csv files
readcsv(x; dir=@__DIR__) = CSV.read(joinpath(dir, x), DataFrame, stringtype=String)
# readin function for parameters; this makes handling easier
readin(x::AbstractDataFrame; default=0,dims=1) = DefaultDict(default,Dict((dims > 1 ? Tuple(row[y] for y in 1:dims) : row[1]) => row[dims+1] for row in eachrow(x)))
readin(x::AbstractString; dir=@__DIR__, kwargs...) = readin(readcsv(x, dir=dir); kwargs...)

### Read in of parameters ###
# We define our sets from the csv files
technologies = readcsv("technologies.csv").technology
fuels = readcsv("fuels.csv").fuel
hour = 1:120
n_hour = length(hour)

# in addition, we now also read a set for the storages
storages = readcsv("storages.csv").storage

# Also, we read our input parameters via csv files
Demand = readin("demand.csv", dims=1)
OutputRatio = readin("outputratio.csv", dims=2)
InputRatio = readin("inputratio.csv", dims=2)
VariableCost = readin("variablecost.csv", dims=1)
InvestmentCost = readin("investmentcost.csv", dims=1)
EmissionRatio = readin("emissionratio.csv", dims=1)
DemandProfile = readin("demandprofile.csv", default=1/n_hour, dims=2)
CapacityFactor = readin("capacity_factors_2018.csv",default=0, dims=2)
TagDispatchableTechnology = readin("tagdispatchabletechnology.csv",default=1,dims=1)

# we need to ensure that all non-variable technologies do have a CapacityFactor of 1 at all times
for t in technologies
    if TagDispatchableTechnology[t] > 0
        for h in hour
            CapacityFactor[h,t] = 1
        end
    end
end

# we can test if solar does still produce during the night
CapacityFactor[2,"SolarPV"]

### Also, we need to read in our additional storage parameters
InvestmentCostStorage = readin("investmentcoststorage.csv",dims=1)
E2PRatio = readin("e2pratio.csv",dims=1)
StorageChargeEfficiency = readin("storagechargeefficiency.csv",dims=2)
StorageDischargeEfficiency = readin("storagedischargeefficiency.csv",dims=2)
StorageLosses = readin("storagelosses.csv",dims=1)

# our emission limit
EmissionLimit = 5000

# define the dictionary for max capacities with specific default value
MaxCapacity = readin("maxcapacity.csv", default=999, dims=1)

# define the dictionary for storage max capacities with specific default value
MaxStorageCapacity = readin("maxstoragecapacity.csv", default=50, dims=1)

### building the model ###
# instantiate a model with an optimizer
ESM = Model(HiGHS.Optimizer)

# this creates our variables
@variable(ESM, TotalCost[technologies] >= 0)
@variable(ESM, FuelProductionByTechnology[hour,technologies, fuels] >= 0)
@variable(ESM, Capacity[technologies] >=0)
@variable(ESM, FuelUseByTechnology[hour,technologies, fuels] >=0)
@variable(ESM, TechnologyEmissions[technologies] >=0)
@variable(ESM,Curtailment[hour,fuels] >=0)

### And we also need to add our new variables for storages
@variable(ESM,StorageEnergyCapacity[s=storages,f=fuels; StorageDischargeEfficiency[s,f]>0]>=0)
@variable(ESM,StorageCharge[s=storages, hour, f=fuels; StorageDischargeEfficiency[s,f]>0]>=0)
@variable(ESM,StorageDischarge[s=storages, hour, f=fuels; StorageDischargeEfficiency[s,f]>0]>=0)
@variable(ESM,StorageLevel[s=storages, hour, f=fuels; StorageDischargeEfficiency[s,f]>0]>=0)
@variable(ESM,TotalStorageCost[storages] >= 0)

## constraints ##
# Generation must meet demand
@constraint(ESM, EnergyBalance[h in hour,f in fuels],
    sum(FuelProductionByTechnology[h,t,f] for t in technologies) + sum(StorageDischarge[s,h,f] for s in storages if StorageDischargeEfficiency[s,f]>0) == 
    Demand[f]*DemandProfile[h,f] + sum(FuelUseByTechnology[h,t,f] for t in technologies) + Curtailment[h,f] + sum(StorageCharge[s,h,f] for s in storages if StorageChargeEfficiency[s,f] > 0)
)

# calculate the total cost
@constraint(ESM, ProductionCost[t in technologies],
    sum(FuelProductionByTechnology[h,t,f] for f in fuels, h in hour) * VariableCost[t] + Capacity[t] * InvestmentCost[t] == TotalCost[t]
)

# limit the production by the installed capacity
@constraint(ESM, ProductionFuntion_disp[h in hour, t in technologies, f in fuels;TagDispatchableTechnology[t]>0],
    OutputRatio[t,f] * Capacity[t] * CapacityFactor[h,t] >= FuelProductionByTechnology[h,t,f]
)

# for variable renewables, the production needs to be always at maximum
@constraint(ESM, ProductionFuntion_res[h in hour, t in technologies, f in fuels;TagDispatchableTechnology[t]==0],
    OutputRatio[t,f] * Capacity[t] * CapacityFactor[h,t] == FuelProductionByTechnology[h,t,f]
)

# define the use by the production
@constraint(ESM, UseFunction[h in hour,t in technologies, f in fuels],
    InputRatio[t,f] * sum(FuelProductionByTechnology[h,t,ff] for ff in fuels) == FuelUseByTechnology[h,t,f]
)

# define the emissions
@constraint(ESM, TechnologyEmissionFunction[t in technologies],
    sum(FuelProductionByTechnology[h,t,f] for f in fuels, h in hour) * EmissionRatio[t] == TechnologyEmissions[t]
)

# limit the emissions
@constraint(ESM, TotalEmissionsFunction,
    sum(TechnologyEmissions[t] for t in technologies) <= EmissionLimit
)

# installed capacity is limited by the maximum capacity
@constraint(ESM, MaxCapacityFunction[t in technologies],
     Capacity[t] <= MaxCapacity[t]
)

### Add your storage constraints here
@constraint(ESM, StorageChargeFunction[s in storages, h in hour, f in fuels; StorageDischargeEfficiency[s,f]>0], 
    StorageCharge[s,h,f] <= StorageEnergyCapacity[s,f]/E2PRatio[s]
)

@constraint(ESM, StorageDischargeFunction[s in storages, h in hour, f in fuels; StorageDischargeEfficiency[s,f]>0], 
    StorageDischarge[s,h,f] <= StorageEnergyCapacity[s,f]/E2PRatio[s]
)

@constraint(ESM, StorageLevelFunction[s in storages, h in hour, f in fuels; h>1 && StorageDischargeEfficiency[s,f]>0], 
    StorageLevel[s,h,f] == StorageLosses[s]*StorageLevel[s,h-1,f] + StorageCharge[s,h,f]*StorageChargeEfficiency[s,f] - StorageDischarge[s,h,f]/StorageDischargeEfficiency[s,f]
)

@constraint(ESM, StorageLevelStartFunction[s in storages, h in hour, f in fuels; h==1 && StorageDischargeEfficiency[s,f]>0], 
    StorageLevel[s,h,f] == 0.5*StorageLosses[s]*StorageEnergyCapacity[s,f]+ StorageCharge[s,h,f]*StorageChargeEfficiency[s,f] - StorageDischarge[s,h,f]/StorageDischargeEfficiency[s,f]
)

#new constraint, replaces StorageAnnualBalanceFunction
@constraint(ESM, StorageLevelEndFunction[s in storages, h in hour, f in fuels; h==n_hour && StorageDischargeEfficiency[s,f]>0], 
    StorageLevel[s,h,f] == 0.5*StorageEnergyCapacity[s,f]
)

@constraint(ESM, MaxStorageLevelFunction[s in storages, h in hour, f in fuels; StorageDischargeEfficiency[s,f]>0], 
    StorageLevel[s,h,f] <= StorageEnergyCapacity[s,f]
)

@constraint(ESM, StorageCostFunction[s in storages], 
    TotalStorageCost[s] == sum(StorageEnergyCapacity[s,f]*InvestmentCostStorage[s] for f in fuels if StorageDischargeEfficiency[s,f]>0)
)

#new constraint: installed storage capacity is limited by the maximum storage capacity
@constraint(ESM, MaxStorageCapacityFunction[s in storages],
     sum(StorageEnergyCapacity[s,f] for f in fuels if StorageDischargeEfficiency[s,f]>0) <= MaxStorageCapacity[s]
)


# the objective function
# total costs should be minimized
@objective(ESM, Min, sum(TotalCost[t] for t in technologies) + sum(TotalStorageCost[s] for s in storages))

# this starts the optimization
# the assigned solver (here HiGHS) will takes care of the solution algorithm
optimize!(ESM)
# reading our objective value
objective_value(ESM)

# some result analysis
value.(FuelProductionByTechnology)
value.(Capacity)
value.(StorageEnergyCapacity)
value.(StorageDischarge)
value.(StorageLevel)
value.(StorageCharge)

df_res_production = DataFrame(Containers.rowtable(value,FuelProductionByTechnology; header = [:Hour, :Technology, :Fuel, :value]))
df_res_capacity = DataFrame(Containers.rowtable(value,Capacity; header = [:Technology, :value]))

df_storage_production = DataFrame(Containers.rowtable(value,StorageDischarge; header = [:Technology, :Hour, :Fuel, :value]))
df_storage_charge = DataFrame(Containers.rowtable(value,StorageCharge; header = [:Technology, :Hour, :Fuel, :value]))
df_storage_level = DataFrame(Containers.rowtable(value,StorageLevel; header = [:Technology, :Hour, :Fuel, :value]))

append!(df_res_production, df_storage_production)

transform!(df_res_production, "Technology" => ByRow(x-> colors[x]) => "Color")
transform!(df_res_capacity, "Technology" => ByRow(x-> colors[x]) => "Color")

# and some plots
bar(
    df_res_capacity.Technology,
    df_res_capacity.value,
    title="Installed Capacity by Technology",
    color=df_res_capacity.Color,
    linewidth=0,
    rotation=90
)

gdf_production_by_fuel = groupby(df_res_production, :Fuel)
sto_charge = Dict((row.Technology, row.Fuel, row.Hour) => row.value for row in eachrow(df_storage_charge))

n_fuels = length(gdf_production_by_fuel)
plts = map(enumerate(pairs(gdf_production_by_fuel))) do (i,(k,v))
    p = groupedbar(
        v.Hour,
        v.value,
        group=v.Technology,
        bar_position=:stack,
        title="$(k[1])",
        linewidth=0,
        color=v.Color,
        legend=i == n_fuels ? (0.1,-0.5) : false,
        bottom_margin=i == n_fuels ? 20mm : 2mm,
        legend_column=5
    )

    d = [Demand[k[1]]*DemandProfile[h,k[1]] for h in hour]
    u = sum(value.(FuelUseByTechnology)[:,t, k[1]] for t in technologies)
    c = [sum(get(sto_charge, (s, k[1], h), 0) for s in storages) for h in hour]
    du = d .+ u.data
    dus = d .+ u.data .+ c
    plot!(p, hour, d, color=:black, linewidth=2, label="Demand")
    plot!(p, hour, du, color=:black, linestyle=:dash, linewidth=2, label="Demand + Use")
    plot!(p, hour, dus, color=:black, linestyle=:dot, linewidth=2, label="Demand + Use + Storage")

    return p
end

plot(plts..., layout=(n_fuels,1), size=(1600,1200))


sto_prod = Dict((row.Technology, row.Fuel, row.Hour) => row.value for row in eachrow(df_storage_production))
sto_lvl = Dict((row.Technology, row.Fuel, row.Hour) => row.value for row in eachrow(df_storage_level))

plt_storage_lvl = map(storages) do s
    
   val = [[get(sto_prod, (s, f, h), 0) for h in hour ] for f in fuels if StorageDischargeEfficiency[s,f]>0]

   color = colors[s]
    p = bar(
        hour,
        val,
        title="$(s)",
        label="Production",
        color=:green,
        ylabel="GW",
        linewidth=0
    )

    val_charg = [[get(sto_charge, (s, f, h), 0) for h in hour ] for f in fuels if StorageDischargeEfficiency[s,f]>0]

    bar!(p,
        hour,
        -val_charg,
        title="$(s)",
        label="Charge",
        color=:red,
        linewidth=0,
        legend=:bottomleft
    )


    val_lvl = [[get(sto_lvl, (s, f, h), 0) for h in hour ] for f in fuels if StorageDischargeEfficiency[s,f]>0]


    p2 = twinx(p)

    plot!(
        p2,
        hour,
        val_lvl,
        label = "Level",
        color=:black,
        linewidth=3,
        ylabel="GWh",
        legend=:bottomright
    )

    

    return p
end

plot(plt_storage_lvl...)