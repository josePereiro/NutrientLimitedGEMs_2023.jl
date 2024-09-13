## ------------------------------------------------------------
@time begin
    using Random
    using MetXGEMs
    using MetXBase
    using ProjFlows
    using Statistics
    using BlobBatches
    using Base.Threads
    using Combinatorics
    using Base.Threads: Atomic
    using NutrientLimitedGEMs_2023
end

# TODO: Add Graphs.jl kind of functionality for getting basic stuff
# ------------------------------------------------------------
include("1_setup.jl")
include("2_utils.jl")

## ------------------------------------------------------------
# delete downset duplicates (this is not exaustive, but it must eliminate the extremes)
let
    n0 = 0 # init
    n1 = Inf # non-ignored count
    lk = ReentrantLock()
    _downset_hashes_ = Set{UInt64}()
    _th_readdir(;n0, n1, nthrs = 10) do bbi, bb

        haskey(bb["meta"], "core_koma.ver") || return :ignore
        haskey(bb["meta"], "core_feasets.ver") || return :ignore
        # haskey(bb["meta"], "core_biomass.ver") || return :ignore
        # haskey(bb["meta"], "core_nut_sp.ver") || return :ignore

        # lock(bb) do
            # load frame
            feasets_frame = bb["core_feasets"]
            core_strip_db = bb["core_strip"]

            do_save = false
            for (feasets_blob0, strip_blob0) in zip(feasets_frame, core_strip_db)
                strip_koset = strip_blob0["koset"]
                for _fealen in keys(feasets_blob0)
                    feasets_blob1 = feasets_blob0[_fealen]
                    _feahash = hash(sort(strip_koset[1:_fealen])::Vector{Int16})
                    if _feahash in _downset_hashes_
                        @assert haskey(feasets_blob0, _fealen)
                        delete!(feasets_blob0, _fealen)
                        @assert !haskey(feasets_blob0, _fealen)
                        do_save = true
                    else
                        lock(lk) do
                            push!(_downset_hashes_, _feahash)
                            # reset limit
                            if length(_downset_hashes_) > Inf
                                empty!(_downset_hashes_)
                            end
                        end
                    end
                end
            end
            do_save && println("do_save! ", do_save)
            do_save && serialize(bb, "core_feasets")
            return :continue
        # end
    end # _th_readdir
end