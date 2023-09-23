## ------------------------------------------------------------
#  Utils
function _run_cmd(cmd::String; 
        ignorestatus = false, 
        detach = false, 
        dir = pwd()
    )
    src = tempname()
    try
        # create script
        rm(src; force = true)
        open(src, "w") do io
            print(io, cmd)
        end
        run(Cmd(["chmod", "777", src]))
        # run
        cmd = Cmd(Cmd(["bash", "-c", src]); ignorestatus, detach, dir)
        run(cmd; wait = true)
    finally
        rm(src; force = true)
    end
    return nothing
end

## ------------------------------------------------------------
let
    projtoml = Base.active_project()
    @sync let
        @async _run_cmd("""pwd"""; dir = @__DIR__)
        for p in 1:3
            @async while true
                _run_cmd(
                    """\
                    julia -t1 --project=$(projtoml) ./3_core_koma.jl -- \
                    "SIMVER:ECOLI-CORE-BEG2007-PHASE_I-0.1.0"
                    """; dir = @__DIR__
                )
            end
            @async while true 
                _run_cmd(
                    """\
                    julia -t1 --project=$(projtoml) ./3_core_koma.jl -- \
                    "SIMVER:ECOLI-CORE-BEG2007-PHASE_II-0.1.0"
                    """; dir = @__DIR__
                )
            end
            @async while true
                 _run_cmd(
                    """\
                    julia -t1 --project=$(projtoml) ./3_core_koma.jl -- \
                    "SIMVER:ECOLI-CORE-BEG2007-PHASE_III-0.1.0"
                    """; dir = @__DIR__
                )
            end
        end
    end
    nothing
end
