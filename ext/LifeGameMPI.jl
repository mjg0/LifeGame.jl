module LifeGameMPI

using MPI, LifeGame



const COMM = MPI.COMM_WORLD

commsize() = MPI.Comm_size(COMM)

commrank() = MPI.Comm_rank(COMM)



function thisrankrows(gridheight::Integer)
    rowsperrank = cld(gridheight, commsize())

    firstrow = min(rowsperrank * commrank() + 1, gridheight)
    lastrow = min((rowsperrank + 1) * commrank(), gridheight)

    if firstrow == gridheight
        firstrow += 1 # indicate that this process is in charge of no columns
    end

    return firstrow, lastrow
end



struct LifeGridMPI{LG<:LifeGame.LifeGrid} <: AbstractMatrix{Bool}
    lifegrid::LG
    height::Int64
    aboverecvbuf::MPI.Buffer
    belowrecvbuf::MPI.Buffer
    abovesendbuf::MPI.Buffer
    belowsendbuf::MPI.Buffer
    requests::MPI.MultiRequest # Halo exchange async requests
end



function LifeGame.MPILifeGrid(m::Integer, n::Integer)
    # Ensure MPI initialization
    if !MPI.Initialized()
        MPI.Init()
    end

    # Which columns will this process be in charge of?
    firstrow, lastrow = thisrankrows(m)

    lg = LifeGrid(lastrow - firstrow + 3, n) # Add two for buffer rows above and below

    aboverecvbuf = MPI.buffer(@view lg.grid[begin, :])
    belowrecvbuf = MPI.buffer(@view lg.grid[end, :])
    abovesendbuf = MPI.buffer(@view lg.grid[begin+1, :])
    belowsendbuf = MPI.buffer(@view lg.grid[end-1, :])

    return LifeGridMPI{typeof(lg)}(
        lg,
        m,
        aboverecvbuf,
        belowrecvbuf,
        abovesendbuf,
        belowsendbuf,
        MPI.MultiRequest(4),
    )
end

function LifeGame.MPILifeGrid(grid::AbstractMatrix{T}) where {T<:Number}
    lg = MPILifeGrid(size(grid)...)
    lg .= grid .!= zero(T)
end



# AbstractArray interface for MPILifeGrid
function Base.size(lg::LifeGridMPI)
    return lg.height, size(lg, 2)
end

Base.@propagate_inbounds function Base.getindex(lg::LifeGridMPI, i::Integer, j::Integer)
    firstcol, lastcol = thisrankcolumns(lg.width)

    value = if j ∈ firstcol:lastcol
        lg.grid[i+1, j-firstcol+1]
    else
        false
    end

    return MPI.Allreduce(value, |, COMM)
end

Base.@propagate_inbounds function Base.setindex!(
    lg::LifeGridMPI,
    value,
    i::Integer,
    j::Integer,
)
    firstcol, lastcol = thisrankcolumns(lg.width)

    if j ∈ firstcol:lastcol
        lg.grid[i+1, j-firstcol+1] = value
    end

    return lg
end



function LifeGame.step!(lg::LifeGridMPI)
    firstrow, lastrow = thisrankrows(lg.height)
    reqs = lg.requests

    # Exchange halos with the process above
    if commrank() > 0 && firstrow ≤ lg.height
        MPI.Isend(lg.abovesendbuf, COMM, reqs[1]; dest = commrank() - 1, tag = 0)
        MPI.Irecv(lg.aboverecvbuf, COMM, reqs[2]; source = commrank() - 1, tag = 1)
    end

    # Exchange halos with the process below
    if commrank() < commsize() && lastrow < lg.height
        MPI.Isend(lg.belowsendbuf, COMM, reqs[3]; dest = commrank() + 1, tag = 1)
        MPI.Irecv(lg.belowrecvbuf, COMM, reqs[4]; source = commrank() + 1, tag = 0)
    end

    # Wait for syncs to finish
    MPI.Waitall(reqs)

    # Update this section
    step!(lg.lifegrid)

    return lg
end



end # module
