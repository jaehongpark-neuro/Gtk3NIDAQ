using Pkg

Pkg.activate(pwd())

# using GLMakie
using Gtk, GtkObservables
# using DataStructures
using JLD
using CircularArrayBuffers
using NIDAQ
using FFTW, DSP
using LinearAlgebra
using Base.Threads



b =  GtkBuilder(filename="realtimeusv.glade")

win = b["win"]
grid = b["wholegrid"]


airbox = b["airbox"]
sngbox = b["sngbox"]

recordstartcb = b["recordstart"]
recordendcb = b["recordend"]


sngcnvs = canvas(UserUnit)
push!(sngbox,sngcnvs )

sngimg = rand(257,3903)


@spawn showall(win)



mutable struct stft_param
    nfft::Int64
    skip::Int64
    nout::Int64
    window::Vector{Float64}
    fin::Vector{Float64}
    tmp::Vector{ComplexF64}
    plan::FFTW.rFFTWPlan{Float64, -1, false, 1, Tuple{Int64}}

    function stft_param()
        nfft = 512
        skip = 128
        nout = (nfft >> 1)+1
        window = hanning(nfft)
        fin = zeros(Float64,nfft)
        tmp = zeros(ComplexF64,nout)
        plan = plan_rfft(fin)

        sp = new(nfft,skip,nout,window,fin,tmp,plan)
        return sp
    end
end

# apply hanning window filter to smoothe out the ends
function sig!(data::Vector{Float64},f::Vector{Float64},ui::UnitRange{Int64},wa::Vector{Float64})
    # data = whole audio 1 to end
    # f = filtered data mutated
    # ui = unit range 1:512 513:1024 ....
    # wa = hanning window
    of =1
    for i = ui
        @inbounds f[of] = data[i]*wa[of]
        of+=1
    end
    nothing
end

# convert into real number
function fout!(fout::Matrix{Float64},tmp::Vector{ComplexF64},offset::Int64)
    for i = 1:length(tmp)
        @inbounds fout[offset+i] = abs2(tmp[i])
    end
    nothing
end

function Stft(data::Vector{Float64},
                        nfft::Int64,skip::Int64,nout::Int64,
                        window::Vector{Float64},
                        fin::Vector{Float64},tmp::Vector{ComplexF64},
                        plan::FFTW.rFFTWPlan{Float64, -1, false, 1, Tuple{Int64}})
    npts = length(data) # for spectrogram not fixed length inputs
    nblocks = Int(floor((npts-nfft)/skip))+1
    # println("$(nblocks)")
    fout = zeros(Float64,nout,nblocks)
    offset = 0
    @inbounds begin for k = 1:nblocks
            period = (1:nfft) .+ (k-1)*skip
            sig!(data,fin,period,window)
            mul!(tmp,plan,fin)
            fout!(fout,tmp,offset)
            offset += nout
        end
    end
    return fout
end

function Stft(data::Vector{Float64},sp::stft_param)
    Stft(data,sp.nfft,sp.skip,sp.nout,sp.window,sp.fin,sp.tmp,sp.plan)
end

function nblock(fs,s)
    sam = Int(floor(fs*s))
    nblocks = Int(floor((sam-512)/128))+1  
end

# const NIDAQ_TID = Ref{Int64}(1)
# if Threads.nthreads() > 1
#     NIDAQ_TID[] = 2
# end

mutable struct Recording
    task::DAQTask{AI}
    channels::Dict{Int64, String}
    fs::Int64
    refresh::Float64

    # samples
    history_seconds::Float64
    history_samples::Int64

    air::CircularArrayBuffer{Float64}
    audio::CircularArrayBuffer{Float64}

    airind::StepRange{Int64, Int64} 
    audioind::StepRange{Int64, Int64}
    
    sng::CircularArrayBuffer{Float64}
    chan::Channel

    function Recording(discnvs::GtkObservables.Canvas{UserUnit},
                    sngimage::Matrix{Float64};
                        channels = Dict(1=>"air",5=>"audio"),
                        fs=250_000,refresh=20,
                        history_seconds = 2)

        # for display e.g. 2s
        history_samples = history_seconds * fs
        
        air = CircularArrayBuffer(zeros(history_samples))
        audio = CircularArrayBuffer(zeros(history_samples))
        sng = CircularArrayBuffer(zeros(257,nblock(fs,history_seconds)))

        # for each channel
        airind = 1:2:Int(fs/refresh*2)
        audioind = 2:2:Int(fs/refresh*2)
        
        chan = Channel(500)
        dev = DefaultDev()
        task = DAQTask{AI}()

        push!(task, dev.channels[AI][1], alias = "air", tcfg = DAQmx.Diff, range = (0.0,3.0))
        push!(task, dev.channels[AI][5], alias = "audio", tcfg = DAQmx.Diff, range = (-2.0,2.0))

        rec = new(task,channels,fs,refresh,history_seconds,history_samples,
        air,audio,airind,audioind,sng,chan)

        pr = stft_param()

        signal_updates(rec,pr,discnvs,sngimage)
        finalizer(rec) do rec
            close(rec.chan)
        end
        return rec
    end
end

function update_disbuffer!(wholesng::CircularArrayBuffer{Float64},sng::Matrix{Float64})
    highcut =0.01
    for i in 1:size(sng,2)
        temp = sng[:,i]
        reverse!(temp)
        temp = clamp.(temp,0.0,highcut)./highcut
        push!(wholesng,temp)
    end
end

function buffer_to_display!(wholesng::CircularArrayBuffer{Float64},
                            display::Matrix{Float64})
    display .= wholesng
    nothing
end

function update_display!(cnvs::GtkObservables.Canvas{UserUnit},
                        wholesng::CircularArrayBuffer{Float64},
                        display::Matrix{Float64})
    
    buffer_to_display!(wholesng,display)

    @guarded draw(cnvs) do widget
        set_coordinates(cnvs,BoundingBox(0, 1, 0, 1))
        ctx = getgc(cnvs)
        copy!(ctx,display)
    end
end

function signal_updates(rec::Recording,sp::stft_param,
    discnvs::GtkObservables.Canvas{UserUnit},
    sngimage::Matrix{Float64},leftover::Int64=0)
    @async begin
        for x in rec.chan
            wholedata = vec(x)
            
            air = wholedata[rec.airind]
            audio = wholedata[rec.audioind]

            
            append!(rec.air, air)
            append!(rec.audio,audio)
            # println("length of whole channle is : $(length(audio))")
            
            currentdata = leftover + length(audio)
            
            ff = Stft(rec.audio[end-currentdata:end],stft_param())
            used = size(ff,2)*sp.skip
            leftover = currentdata - used
            println("size of blocks is : $(size(ff,2))")

            update_disbuffer!(rec.sng,ff)
            update_display!(discnvs,rec.sng,sngimage)

        end
    end
end


function (rec::Recording)()
    @spawn record2!(rec.task, $(rec.fs), $(rec.refresh); feed = rec.chan)
end

function recordstart(w)
    wait(REC())
end

id = signal_connect(recordstart, recordstartcb, "clicked")

function recordend(w)
    NIDAQ.stop(REC.task)
end

id = signal_connect(recordend,recordendcb,"clicked")

# initialize

hs = 2
fs = 250_000
sngimg = rand(257,nblock(fs,hs))
REC = Recording(sngcnvs,sngimg;fs=fs,refresh=20,history_seconds=hs)

# actual recording start
wait(REC())

# stop
NIDAQ.stop(REC.task)

REC.chan
put!(REC.chan,ones(25000))


@async for i in 1:1000
    put!(REC.chan,clamp.(rand(25000),0.0,0.1)./0.1)
    sleep(0.05)
end

put!(REC.chan,ones(25000))


