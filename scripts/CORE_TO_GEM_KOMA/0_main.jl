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

## ------------------------------------------------------------
# (julia -t5 --project scripts/CORE_TO_GEM_KOMA/6_fva.jl -- "SIMVER:ECOLI-CORE-BEG2007-PHASE_III-0.1.0") &
@sync let
    projtoml = Base.active_project()
    @async _run_cmd("""\
        while true; do
            julia -t5 --project=$(projtoml) ./4_koma_strip.jl -- "SIMVER:ECOLI-CORE-BEG2007-PHASE_I-0.1.0";\
            julia -t5 --project=$(projtoml) ./5_core_feasets.jl -- "SIMVER:ECOLI-CORE-BEG2007-PHASE_I-0.1.0";\
            julia -t5 --project=$(projtoml) ./5.1_enumerate.jl -- "SIMVER:ECOLI-CORE-BEG2007-PHASE_I-0.1.0";\
            julia -t5 --project=$(projtoml) ./6_fva.jl -- "SIMVER:ECOLI-CORE-BEG2007-PHASE_I-0.1.0";\
            sleep 100
        done
        """; dir = @__DIR__
    )
    @async _run_cmd("""\
        while true; do
            julia -t5 --project=$(projtoml) ./4_koma_strip.jl -- "SIMVER:ECOLI-CORE-BEG2007-PHASE_II-0.1.0";\
            julia -t5 --project=$(projtoml) ./5_core_feasets.jl -- "SIMVER:ECOLI-CORE-BEG2007-PHASE_II-0.1.0";\
            julia -t5 --project=$(projtoml) ./5.1_enumerate.jl -- "SIMVER:ECOLI-CORE-BEG2007-PHASE_II-0.1.0";\
            julia -t5 --project=$(projtoml) ./6_fva.jl -- "SIMVER:ECOLI-CORE-BEG2007-PHASE_II-0.1.0";\
            sleep 100
        done
        """; dir = @__DIR__
    )
    @async _run_cmd("""\
        while true; do
            julia -t5 --project=$(projtoml) ./4_koma_strip.jl -- "SIMVER:ECOLI-CORE-BEG2007-PHASE_III-0.1.0";\
            julia -t5 --project=$(projtoml) ./5_core_feasets.jl -- "SIMVER:ECOLI-CORE-BEG2007-PHASE_III-0.1.0";\
            julia -t5 --project=$(projtoml) ./5.1_enumerate.jl -- "SIMVER:ECOLI-CORE-BEG2007-PHASE_III-0.1.0";\
            julia -t5 --project=$(projtoml) ./6_fva.jl -- "SIMVER:ECOLI-CORE-BEG2007-PHASE_III-0.1.0";\
            sleep 100
        done
        """; dir = @__DIR__
    )
end

## ------------------------------------------------------------
# let
#     projtoml = Base.active_project()

#     # @sync 
#     let
#         # _run_cmd("""pwd"""; dir = @__DIR__)
#         for p in 1:1
#             while true
#                 for srcscript in [
#                     "./4_koma_strip.jl", 
#                 ]
#                     for arg in [
#                             "SIMVER:ECOLI-CORE-BEG2007-PHASE_I-0.1.0",
#                             "SIMVER:ECOLI-CORE-BEG2007-PHASE_II-0.1.0",
#                             "SIMVER:ECOLI-CORE-BEG2007-PHASE_III-0.1.0",
#                         ]
#                         @async let
#                             """\
#                             julia -t1 --project=$(projtoml) $(srcscript) -- $(arg)
#                             """ |> println
#                             for i in 1:10
#                                 pritln("computing")
#                                 sleep(3)
#                             end
#                         end
#                         # @async _run_cmd(
#                         #     """\
#                         #     julia -t1 --project=$(projtoml) $(srcscript) -- $(arg)
#                         #     """; dir = @__DIR__
#                         # )
#                     end
#                 end
#             end # for arg 
#         end
#     end
#     nothing
# end
