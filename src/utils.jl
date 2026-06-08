highbit(T::Type{<:Unsigned}) = one(T) << (8*sizeof(T)-1)
highbit(::T) where {T} = highbit(T)

lowbit(T::Type{<:Unsigned}) = one(T)
lowbit(::T) where {T} = lowbit(T)

nbits(T::Type{<:Unsigned}) = 8*sizeof(T)
nbits(::T) where {T} = nbits(T)
