# TODO: package this
using Base.Threads

## ------------------------------------------------------------------------
t = ReentrantLock()

_info(msg; kwargs...) = lock(t) do
    @info(msg, kwargs...)
end

progress = ["" for t in 1:nthreads()]

## ------------------------------------------------------------------------
# show progress
@async while true
    for (tid, file) in enumerate(progress)
        isempty(file) && continue
        @info("Doing ", 
            tid,
            file = basename(file), 
            size = filesize(file)
        )
    end
    sleep(3.0)
end

## ------------------------------------------------------------------------
# convert
# dir = "/Users/Pereiro/Downloads/Revolutions Podcast - S4_ Haitian Revolution"
dir = "/Users/Pereiro/Downloads/Revolutions Podcast - S5_ Simón Bolívar"
@threads for src in readdir(dir; join = true)

    tid = threadid()

    endswith(src, ".mp4") || continue
    dest = replace(src, ".mp4" => ".aac")
    
    rm(dest; force = true)
    proc = run(`ffmpeg -i $(src) $(dest)`; wait = false)
    progress[tid] = dest
    wait(proc)
    
    progress[tid] = ""

end

# bash
# for s in *.mp4; do
#     rm -f "${s%.mp4}.aac"
#     ffmpeg -threads 3 -i "${s}" "${s%.mp4}.aac"
#     exit
# done