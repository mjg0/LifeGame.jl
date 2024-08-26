module LifeGameMPI

using LifeGame, MPI

export LifeGridMPI, step!



# Determine the sub-width of the nth slice (of N, starting at 0) given a total width
function slicewidth(width, N, n=0)
    gridwidth = cld(width, LifeGame.CELLS_PER_CLUSTER)
    subwidth = cld(gridwidth, N) * LifeGame.CELLS_PER_CLUSTER
    return min((n+1)*subwidth, width) - n*subwidth
end



function sliceindex(width, N, j)
    return divrem(j-1, slicewidth(width, N)).+(0, 1) # Add 1 for 1-based j index
end



struct LifeGridMPI <: AbstractMatrix{Bool}
    lifegrid::LifeGrid
    width::Int64
    commsize::Int64
    commrank::Int64
    leftsendbuf ::MPI.Buffer # MPI buffer
    leftrecvbuf ::MPI.Buffer # MPI buffer
    rightsendbuf::MPI.Buffer # MPI buffer
    rightrecvbuf::MPI.Buffer # MPI buffer
    halorequests::MPI.MultiRequest # Halo exchange async requests

    function LifeGridMPI(m::Integer, n::Integer)
        if !MPI.Initialized()
            MPI.Init()
        end

        commsize = MPI.Comm_size(MPI.COMM_WORLD)
        commrank = MPI.Comm_rank(MPI.COMM_WORLD)

        thisrankwidth = slicewidth(n, commsize, commrank)

        lg = LifeGrid(m, thisrankwidth)

        leftsendbuf  = MPI.Buffer(@view lg.grid[:,begin+1])
        leftrecvbuf  = MPI.Buffer(@view lg.grid[:,begin])
        rightsendbuf = MPI.Buffer(@view lg.grid[:,end-1])
        rightrecvbuf = MPI.Buffer(@view lg.grid[:,end])

        halorequests = MPI.MultiRequest(4)

        return new(lg, n, commsize, commrank,
                   leftsendbuf, leftrecvbuf, rightsendbuf, rightrecvbuf, halorequests)
    end

    function LifeGridMPI(grid::Union{BitMatrix,AbstractMatrix{Bool}})
        lg = LifeGridMPI(size(grid)...)
        lg .= grid
    end

    LifeGridMPI(grid::AbstractMatrix{<:Number}) = LifeGridMPI(grid.!=0)
end



# AbstractArray interface for LifeGridMPI
function Base.size(lg::LifeGridMPI)
    return size(lg.lifegrid, 1), lg.width
end

Base.@propagate_inbounds function Base.getindex(lg::LifeGridMPI, i::Integer, j::Integer)
    process, jdx = sliceindex(lg.width, lg.commsize, j)
    jdx = min(jdx, size(lg.lifegrid, 2))
    value = lg.lifegrid[i,jdx] # dummy for all but `process`
    return MPI.Bcast(value, process, MPI.COMM_WORLD)
end

Base.@propagate_inbounds function Base.setindex!(lg::LifeGridMPI, value,
                                                 i::Integer, j::Integer)
    process, jdx = sliceindex(lg.width, lg.commsize, j)
    if process == lg.commrank
        lg.lifegrid[i,jdx] = value
    end
    @boundscheck if process >= lg.commsize
        BoundsError(lg, (i, j))
    end
    return value
end



function step!(lg::LifeGridMPI; kw...)
    # Exchange halos before running the update
    reqs = lg.halorequests
    @inbounds if lg.commrank > 0             # Left halo
        MPI.Isend( lg.leftsendbuf,  MPI.COMM_WORLD, reqs[1]; dest  =lg.commrank-1, tag=0)
        MPI.Irecv!(lg.leftrecvbuf,  MPI.COMM_WORLD, reqs[2]; source=lg.commrank-1, tag=1)
    end
    @inbounds if lg.commrank < lg.commsize-1 # Right halo
        MPI.Isend( lg.rightsendbuf, MPI.COMM_WORLD, reqs[3]; dest  =lg.commrank+1, tag=1)
        MPI.Irecv!(lg.rightrecvbuf, MPI.COMM_WORLD, reqs[4]; source=lg.commrank+1, tag=0)
    end
    MPI.Waitall(reqs)

    # Update and return
    @inbounds LifeGame.step!(lg.lifegrid; kw...)
    return lg
end



end # module
