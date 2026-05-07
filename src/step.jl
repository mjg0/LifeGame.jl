export step!



Base.@propagate_inbounds @inline function
updatedchunkhalos(chunk::AbstractVector{C}, ::Type{H}) where {C, H}
    Cbits = 8*sizeof(C)
    Hbits = 8*sizeof(H)

    lhalo = zero(H)
    rhalo = zero(H)

    lshift = Cbits-2
    rshift = 1

    for k in 1:Hbits
        cell = chunk[k]
        lhalo |= H((cell >> lshift) & one(C)) << (Hbits-k)
        rhalo |= H((cell >> rshift) & one(C)) << (Hbits-k)
    end

    return lhalo, rhalo
end



#Base.@propagate_inbounds @inline function
#updatedchunkhalos(lg::LifeGrid{R, C, H}, i, j) where {R, C, H}
#    return updatedchunkhalos(view(lg.grid, i:i+8*sizeof(H)-1, j), H)
#end



"""
    step!(lg::LifeGrid; parallel::Bool)

Update `lg` one generation according to the [`rule`](@ref) associated with it.

All cells outside of the grid boundary are fixed at zero.

By default, `parallel` is set to true if the grid backing `lg` has multiple columns.

A generic algorithm for updating each cluster in the grid is used by default. The compiler
does a decent job of optimizing for most rules, but hand-tuning the cluster update function
can improve performance by 10% for some rules. See the extended help for
[`LifeGame.updatedcluster`](@ref) for instructions on specializing the cluster update.
Specializations are provided for commonly used rules (`B3/S23`, `B36/S23`, and `B2/s`).
"""
function step!(lg::LifeGrid; parallel=size(lg.grid, 2)>3)
    if parallel
        return stepraw!(lg)
    else
        disable_polyester_threads() do
            return stepraw!(lg)
        end
    end
end



"""

`out` is 2 longer than `in` and doesn't get its first or last cells updated
"""
Base.@propagate_inbounds @inline function
updatehalos!(out::AbstractVector{T}, in::AbstractVector{T}, lhalo::H, rhalo::H) where {T, H}
    Hbits = 8*sizeof(H)
    Tbits = 8*sizeof(T)
    centermask = (typemax(T) << 2) >> 1

    @simd for i in 1:Hbits
        lbit = T(lhalo >> (Hbits-i) & one(H)) << (Tbits-1)
        rbit = T(rhalo >> (Hbits-i) & one(H))

        out[i+1] = (in[i] & centermask) | lbit | rbit
    end
end



#Base.@propagate_inbounds @inline function updategridchunk!(lg::LifeGrid{R,C,H}, i, j, cur::AbstractVector{C}) where {R,C,H}
#    N = 8*sizeof(H)
#    @simd for k in 1:N
#        lg.grid[i+k-1,j] = updatedcluster(cur[k], cur[k+1], cur[k+2], R)
#    end
#end



"""

`in` should be two elements longer than `out`
"""
Base.@propagate_inbounds @inline function
updategridchunk!(out::AbstractVector, in::AbstractVector, rule::R, ::Type{H}) where {R, H}
    @simd for i in 1:8*sizeof(H)
        out[i] = updatedcluster(in[i], in[i+1], in[i+2], rule)
    end
end



Base.@propagate_inbounds @inline function
gridchunk(lg::LifeGrid{R, C, H}, I, J) where {R, C, H}
    Hbits = 8*sizeof(H)
    i = (I-1)*Hbits+2
    return view(lg.grid, i:i+Hbits-1, J)
end



Base.@propagate_inbounds @inline function
zerotrailing!(cells::AbstractVector, shift, ::Type{H}) where H
    @simd for i in 1:8*sizeof(H)
        cells[i] = (cells[i] >> shift) << shift
    end
end



function stepraw!(lg::LifeGrid{R, C, H}) where {R, C, H}
    Cbits = 8*sizeof(C)
    Hbits = 8*sizeof(H)

    grid = lg.grid

    J1 = firstindex(lg.grid, 2)+1
    J2 = lastindex(lg.grid, 2)-1

    bufflen = Hbits+2

    endshift = mod(-size(lg, 2), Cbits-2)+1

    @inbounds @batch for J in J1:J2
        current = lg.colbuffers1[Threads.threadid()]
        next = lg.colbuffers2[Threads.threadid()]

        inhalosleft   = view(lg.righthalos[1], :, J-1)
        inhalosright  = view(lg.lefthalos[ 1], :, J+1)
        outhalosleft  = view(lg.lefthalos[ 2], :, J)
        outhalosright = view(lg.righthalos[2], :, J)

        above = zero(C)
        updatehalos!(current, gridchunk(lg, 1, J), inhalosleft[1], inhalosright[1])

        # Update this row
        for I in firstindex(inhalosleft):lastindex(inhalosleft)-1
            updatehalos!(next, gridchunk(lg, I+1, J), inhalosleft[I+1], inhalosright[I+1])

            current[1] = above
            current[bufflen] = next[2]

            updategridchunk!(gridchunk(lg, I, J), current, R, H)

            # Zero trailing cells
            if J==J2
                zerotrailing!(gridchunk(lg, I, J), endshift, H)
            end

            outhalosleft[I], outhalosright[I] = updatedchunkhalos(gridchunk(lg, I, J), H)

            above = current[bufflen-1]

            current, next = next, current
        end

        I = lastindex(inhalosleft)

        current[1] = above
        current[bufflen] = zero(C)

        updategridchunk!(gridchunk(lg, I, J), current, R, H)

        if J==J2
            zerotrailing!(gridchunk(lg, I, J), endshift, H)
        end

        # Zero the last+1 cell
        grid[lg.height+2,J] = zero(C)

        outhalosleft[I], outhalosright[I] = updatedchunkhalos(gridchunk(lg, I, J), H)
    end

    lg.lefthalos[ 1], lg.lefthalos[ 2] = lg.lefthalos[ 2], lg.lefthalos[ 1]
    lg.righthalos[1], lg.righthalos[2] = lg.righthalos[2], lg.righthalos[1]

    return lg
end