## ------------------------------------------------------------------
@time begin
    using MetXBase
    using MetXOptim
    using MetXNetHub
    using MetXCultureHub
    using Gurobi
    using ContextDBs
end

# ------------------------------------------------------------------
include("1_setup_sim.jl")

## ------------------------------------------------------------------
# NET0
## ------------------------------------------------------------------
@tempcontext ["NET0" => v"0.1.0"] let

    # base model
    net0 = pull_net("SysBioChalmers_Human_GEM")

    # I'm opening all intakes in Ham medium so they are more than enough for a growth of 0.1 h^{-1}.
    bounds!(net0, "MAR09034", -1000.0, 1000.0) # glucose (metid: MAM01965e)
    bounds!(net0, "MAR09063", -1000.0, 1000.0) # glutamine (metid: MAM01975e)
    bounds!(net0, "MAR09046", -1000.0, 1000.0) # valine (metid: MAM03135e)
    bounds!(net0, "MAR09065", -1000.0, 1000.0) # cysteine (metid: MAM01628e)
    bounds!(net0, "MAR09068", -1000.0, 1000.0) # proline (metid: MAM02770e)
    bounds!(net0, "MAR09062", -1000.0, 1000.0) # asparagine (metid: MAM01369e)
    bounds!(net0, "MAR09039", -1000.0, 1000.0) # isoleucine (metid: MAM02184e)
    bounds!(net0, "MAR09040", -1000.0, 1000.0) # leucine (metid: MAM02360e)
    bounds!(net0, "MAR09066", -1000.0, 1000.0) # arginine (metid: MAM01365e)
    bounds!(net0, "MAR09038", -1000.0, 1000.0) # histidine (metid: MAM02125e)
    bounds!(net0, "MAR09041", -1000.0, 1000.0) # lysine (metid: MAM02426e)
    bounds!(net0, "MAR09042", -1000.0, 1000.0) # methionine (metid: MAM02471e)
    bounds!(net0, "MAR09043", -1000.0, 1000.0) # phenylalanine (metid: MAM02724e)
    bounds!(net0, "MAR09045", -1000.0, 1000.0) # tryptophan (metid: MAM03089e)
    bounds!(net0, "MAR09064", -1000.0, 1000.0) # tyrosine (metid: MAM03101e)
    bounds!(net0, "MAR09061", -1000.0, 1000.0) # alanine (metid: MAM01307e)
    bounds!(net0, "MAR09067", -1000.0, 1000.0) # glycine (metid: MAM01986e)
    bounds!(net0, "MAR09069", -1000.0, 1000.0) # serine (metid: MAM02896e)
    bounds!(net0, "MAR09044", -1000.0, 1000.0) # threonine (metid: MAM02993e)
    bounds!(net0, "MAR09070", -1000.0, 1000.0) # aspartate (metid: )
    bounds!(net0, "MAR09071", -1000.0, 1000.0) # glutamate (metid: MAM01974e)
    bounds!(net0, "MAR09109", -1000.0, 1000.0) # biotin (metid: MAM01401e)
    bounds!(net0, "MAR09145", -1000.0, 1000.0) # pantothenate (metid: MAM02680e)
    bounds!(net0, "MAR09167", -1000.0, 1000.0) # lipoic acid (metid: MAM02394e)
    bounds!(net0, "MAR09074", -1000.0, 1000.0) # sulfate (metid: MAM02946e)
    bounds!(net0, "MAR09404", -1000.0, 1000.0) # retinoate (metid: MAM02833e)
    bounds!(net0, "MAR09035", -1000.0, 1000.0) # linoleate (metid: MAM02387e)
    bounds!(net0, "MAR09036", -1000.0, 1000.0) # linolenate (metid: MAM02389e)
    bounds!(net0, "MAR09083", -1000.0, 1000.0) # choline (metid: MAM01513e)
    bounds!(net0, "MAR09361", -1000.0, 1000.0) # inositol (metid: MAM02171e)
    bounds!(net0, "MAR09143", -1000.0, 1000.0) # riboflavin (metid: MAM02842e)
    bounds!(net0, "MAR09378", -1000.0, 1000.0) # nicotinamide (metid: MAM02583e)
    bounds!(net0, "MAR09144", -1000.0, 1000.0) # pyridoxine (metid: MAM02817e)
    bounds!(net0, "MAR09423", -1000.0, 1000.0) # thymidine (metid: MAM02996e)
    bounds!(net0, "MAR09151", -1000.0, 1000.0) # alpha-tocopherol (metid: MAM01327e)
    bounds!(net0, "MAR09146", -1000.0, 1000.0) # folate (metid: MAM01830e)
    bounds!(net0, "MAR09358", -1000.0, 1000.0) # hypoxanthine (metid: MAM02159e)
    bounds!(net0, "MAR09269", -1000.0, 1000.0) # aquacob(III)alamin (metid: MAM01361e)
    bounds!(net0, "MAR09153", -1000.0, 1000.0) # gamma-tocopherol (metid: MAM01935e)
    bounds!(net0, "MAR09159", -1000.0, 1000.0) # thiamin (metid: MAM02982e)
    bounds!(net0, "MAR09048", -1000.0, 1000.0) # O2 (metid: MAM02630e)
    bounds!(net0, "MAR09047", -1000.0, 1000.0) # H2O (metid: MAM02040e)
    bounds!(net0, "MAR09072", -1000.0, 1000.0) # Pi (metid: MAM02751e)
    bounds!(net0, "MAR09076", -1000.0, 1000.0) # Fe2+ (metid: MAM01821e)

    @stage! "net0" => CacheRef(net0)

    nothing
end

## ------------------------------------------------------------------
# Essential reactions in the medium
@tempcontext ["NET0.ESSENTIALS.EXCHANGES" => v"0.1.0"] let

    # Load net0
    @stage! net0_db = query(["ROOT", "NET0"])
    net0 = net0_db["net0"][]

    # globals
    @stage! globals_db = query(["ROOT", "GLOBALS"])
    
    # opm
    solver = globals_db["LP_SOLVER"]
    opm = FBAOpModel(net0, solver)
    
    # essentials
    @stage! biom_th = 1e-2
    @stage! non_essentials = String[]
    @stage! essentials = String[]
    for ri in eachindex(net0.rxns)
        net0.subSystems[ri] == "Exchange/demand reactions" || continue
        
        rxn = net0.rxns[ri]
        _lb, _ub = bounds(opm, rxn)
        _lb < 0 || continue # only intakes
        bounds!(opm, rxn, 0.0, 0.0)
        
        @time optimize!(opm)
        biom = solution(opm, "MAR13082")
        if biom < biom_th
            println("."^60)
            println("Essential! ")
            summary(net0, rxn)
            push!(essentials, rxn)
        else
            push!(non_essentials, rxn)
        end
        
        bounds!(opm, rxn, _lb, _ub)
    end

end

## ------------------------------------------------------------------
# NET1
# I will reduce the essential intakes so they are at the minimum for a given maximum biomass
@tempcontext ["NET1" => v"0.1.0"] let
    
    # Load net0
    @stage! net0_db = query(["ROOT", "NET0"])
    net0 = net0_db["net0"][]

    # globals
    @stage! globals_db = query(["ROOT", "GLOBALS"])
    
    # essentials
    @stage! essentials_db = query(["ROOT", "NET0.ESSENTIALS.EXCHANGES"])
    essentials = essentials_db["essentials"]

    
    # opm
    solver = globals_db["LP_SOLVER"]
    opm = FBAOpModel(net0, solver)

    # minimize essentials
    @stage! biom_max = 0.1 # h^{-1}
    bounds!(opm, "MAR13082", biom_max, biom_max)
    set_linear_obj!(opm, essentials, MAX_SENSE)
    optimize!(opm)
    lb_ = solution(opm, essentials)
    @assert all(lb_ .< 0)

    # set bounds into network
    lb!(net0, essentials, lb_)

    # test biomass (note biomass is unbounded)
    opm = FBAOpModel(net0, solver)
    set_linear_obj!(opm, "MAR13082", MAX_SENSE)
    optimize!(opm)
    @show solution(opm, "MAR13082")

    @stage! "net1" => CacheRef(net0)

    nothing
end

## ------------------------------------------------------------------
## ------------------------------------------------------------------
## ------------------------------------------------------------------

# ## ------------------------------------------------------------------
# # NET2
# # The boxed version of NET1
# ## ------------------------------------------------------------------
# @tempcontext ["NET2" => v"0.1.0"] let

#     # Load net0
#     net0_db = query(["ROOT", "WIP", "NET0"])
#     net0 = net0_db["net0"][]

#     # globals
#     globals_db = query(["ROOT", "WIP", "GLOBALS"])

#     # boxing    
#     @stage! cid = (:NET1, hash(net0))
#     @stage! solver = globals_db["LP_SOLVER"]
#     @stage! nths = globals_db["NTHREADS"]
#     cfile, net1 = withcachedat(PROJ, :set!, cid) do
#         # TODO: check no intakes are blocked
#         _net1 = box(net0, solver; nths, verbose = true)
#         return CacheRef(_net1) # return just the ref
#     end

#     @stage! cfile 
#     @stage! net1

#     nothing

# end



## ------------------------------------------------------------------
# To compute the shadow price of exch i
# 1. find max intake of i given the biomass.
# 2. give the biomass if the intake is reduced by a fraction. 
# NOTE: Note that this is more strict than just perturving around any bound. 
# If a baund is bigger than the maximum the shadow price will always be 0. 
# Now, perturbing around the maximum, the shadow price will be zero only if the network have the flexibility to move the load (biomass) around.
# Essential mets should always have non zero shadow price.
