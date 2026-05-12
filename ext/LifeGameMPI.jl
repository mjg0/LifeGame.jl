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



function LifeGame.MPILifeGrid(m::Integer, n::Integer; rule = "B3/S23")
    # Ensure MPI initialization
    if !MPI.Initialized()
        MPI.Init()
    end

    # Which columns will this process be in charge of?
    firstrow, lastrow = thisrankrows(m)
    M = lastrow - firstrow + 1

    # Internal LifeGrid, with an extra rows at the top and bottom
    lg = LifeGrid(M + 2, n; rule = rule)

    aboverecvbuf = MPI.Buffer(@view lg.grid[begin, :])
    belowrecvbuf = MPI.Buffer(@view lg.grid[end, :])
    abovesendbuf = MPI.Buffer(@view lg.grid[begin+1, :])
    belowsendbuf = MPI.Buffer(@view lg.grid[end-1, :])

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

function LifeGame.MPILifeGrid(grid::AbstractMatrix{T}; kw...) where {T<:Number}
    lg = MPILifeGrid(size(grid)...; kw...)
    lg .= grid .!= zero(T)
end



# AbstractArray interface for MPILifeGrid
function Base.size(lg::LifeGridMPI)
    return lg.height, size(lg.lifegrid, 2)
end

Base.@propagate_inbounds function Base.getindex(lg::LifeGridMPI, i::Integer, j::Integer)
    firstrow, lastrow = thisrankrows(lg.height)

    value = if i ∈ firstrow:lastrow
        lg.lifegrid[i-firstrow+2, j]
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
    firstrow, lastrow = thisrankrows(lg.height)

    if i ∈ firstrow:lastrow
        lg.lifegrid[i-firstrow+1, j] = value
    end

    return lg
end



function LifeGame.step!(lg::LifeGridMPI)
    firstrow, lastrow = thisrankrows(lg.height)
    reqs = lg.requests

    # Exchange halos with the process above
    if commrank() > 0 && firstrow ≤ lg.height
        MPI.Isend(lg.abovesendbuf, COMM, reqs[1]; dest = commrank() - 1, tag = 0)
        MPI.Irecv!(lg.aboverecvbuf, COMM, reqs[2]; source = commrank() - 1, tag = 1)
    end

    # Exchange halos with the process below
    if commrank() < commsize()-1 && lastrow < lg.height
        MPI.Isend(lg.belowsendbuf, COMM, reqs[3]; dest = commrank() + 1, tag = 1)
        MPI.Irecv!(lg.belowrecvbuf, COMM, reqs[4]; source = commrank() + 1, tag = 0)
    end

    # Wait for syncs to finish
    MPI.Waitall(reqs)

    # Update this section
    step!(lg.lifegrid)

    return lg
end



end # module
