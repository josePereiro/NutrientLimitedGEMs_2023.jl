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

# ## ------------------------------------------------------------
# # 3_core_koma
# let
#     projtoml = Base.active_project()
#     @sync let
#         @async _run_cmd("""pwd"""; dir = @__DIR__)
#         for p in 1:1
#             @async while true
#                 _run_cmd(
#                     """\
#                     julia -t1 --project=$(projtoml) ./3_core_koma.jl -- \
#                     "SIMVER:ECOLI-CORE-BEG2007-PHASE_1"
#                     """; dir = @__DIR__
#                 )
#             end
#             @async while true 
#                 _run_cmd(
#                     """\
#                     julia -t1 --project=$(projtoml) ./3_core_koma.jl -- \
#                     "SIMVER:ECOLI-CORE-BEG2007-PHASE_2"
#                     """; dir = @__DIR__
#                 )
#             end
#             @async while true
#                  _run_cmd(
#                     """\
#                     julia -t1 --project=$(projtoml) ./3_core_koma.jl -- \
#                     "SIMVER:ECOLI-CORE-BEG2007-PHASE_3"
#                     """; dir = @__DIR__
#                 )
#             end
#         end
#     end
#     nothing
# end

## ------------------------------------------------------------
# 4_koma_strip
let
projtoml = Base.active_project()
@sync let
    @async _run_cmd("""pwd"""; dir = @__DIR__)
    for p in 1:1
        @async while true
            _run_cmd(
                """\
                julia -t1 --project=$(projtoml) ./4_koma_strip.jl -- \
                "SIMVER:ECOLI-CORE-BEG2007-PHASE_1"
                """; dir = @__DIR__
            )
        end
        @async while true 
            _run_cmd(
                """\
                julia -t1 --project=$(projtoml) ./4_koma_strip.jl -- \
                "SIMVER:ECOLI-CORE-BEG2007-PHASE_2"
                """; dir = @__DIR__
            )
        end
        @async while true
             _run_cmd(
                """\
                julia -t1 --project=$(projtoml) ./4_koma_strip.jl -- \
                "SIMVER:ECOLI-CORE-BEG2007-PHASE_3"
                """; dir = @__DIR__
            )
        end
    end
end
nothing
end
