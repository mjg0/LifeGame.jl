"""
    BitWriter{I<:IO}

A buffered writer for writing values a few bits at a time to an `IO`.

Bits are packed from high to low within each byte. Any trailing partial byte is written when
[`flushbits!`](@ref) is called or when [`bitwriter`](@ref) exits.

---

    BitWriter(io::IO; bufsize=1024^2)

Return a `BitWriter` wrapping `io`, with a buffer defaulting to 1 MiB.
"""
mutable struct BitWriter{I<:IO}
    io::I
    bytes::Vector{UInt8}
    idx::Int64 # 0-based for bitshifting convenience

    function BitWriter(io::I; bufsize=1024^2) where {I<:IO}
        return new{I}(io, Vector{UInt8}(undef, bufsize), 0)
    end
end



"""
    flushbits!(bw::BitWriter)

Write all buffered bits in `bw` to its underlying `IO`.

If the buffered bit count is not a multiple of 8, the final byte is written with its unused
rightmost bits set to zero.
"""
function flushbits!(bw::BitWriter)
    if bw.idx > 0
        write(bw.io, @view bw.bytes[1:cld(bw.idx, 8)])
        bw.idx = 0
    end
end



"""
    writebits!(bw::BitWriter, bits::Unsigned, n::Integer)

Write the leftmost `n` bits of `bits` to `bw` and return `n`.

`n` must be no larger than the number of bits in `eltype(bits)`.
"""
Base.@propagate_inbounds function writebits!(
    bw::BitWriter,
    bits::T,
    n::Integer,
) where T<:Unsigned
    @boundscheck if n > nbits(T)
        throw(ArgumentError("not enough bits in input"))
    end

    # Return the number of bits written
    ret = n

    # Iterate until all requested bits have been written
    while n > 0
        # Zero the upcoming bit if we're starting in on it
        i = bw.idx ÷ 8 + 1
        if bw.idx % 8 == 0
            bw.bytes[i] = 0x00
        end

        # Get the bits to write this iteration
        writecount = min(8 - (bw.idx % 8), n)
        mask = typemax(T) << (nbits(T) - writecount)
        newbits = bits & mask

        # Write them to bw.bytes
        shift = nbits(T) - 8 + bw.idx % 8
        bw.bytes[i] |= UInt8(newbits >> shift)

        # Flush the buffer if we've reached the end
        bw.idx += writecount
        if bw.idx == length(bw.bytes) * 8
            flushbits!(bw)
        end

        # Prepare for the next iteration
        n -= writecount
        bits = bits << writecount
    end

    return ret
end



"""
    bitwriter(f, io::IO; kw...)

Call `f` with a temporary [`BitWriter`](@ref) wrapping `io`.

`kw` are passed to the [`BitWriter`](@ref) constructor.

Any buffered bits are flushed before `bitwriter` returns, even if `f` throws.

# Examples

Write padded unsigned integers to a file as packed bits, discarding one bit of padding on
each side of each element:

```julia
f = open("grid.dat", "w")
bitwriter(f) do bw
    for elem in mypaddedgrid
        writebits!(bw, elem << 1, 8*sizeof(elem)-2)
    end
end
```
"""
Base.@propagate_inbounds function bitwriter(f, io::IO; kw...)
    bw = BitWriter(io; kw...)
    try
        return f(bw)
    finally
        flushbits!(bw)
    end
end



"""
    BitReader{I<:IO}

A buffered reader for reading unsigned values a few bits at a time from an `IO`.

Bits are read from left to right within each byte.

---

    BitReader(io::IO; bufsize=1024^2, maxbits=typemax(Int64))

Return a `BitReader` wrapping `io`, with a buffer defaulting to 1 MiB.

The returned `BitReader` will read at most `cld(maxbits, 8)` bytes from `io`.
"""
mutable struct BitReader{I<:IO}
    io::I
    bytes::Vector{UInt8}
    idx::Int64
    maxbits::Int64
    readbits::Int64

    function BitReader(io::I; bufsize=1024^2, maxbits=typemax(Int64)) where {I<:IO}
        br = new{I}(io, Vector{UInt8}(undef, bufsize), 0, maxbits, 0)
        refillbits!(br)
        return br
    end
end



"""
    refillbits!(br::BitReader)

Refill `br`'s byte buffer from its underlying `IO`.

At most `br.maxbits` in total are read across all calls to `refillbits!`.
"""
function refillbits!(br::BitReader)
    br.idx = 0

    readcount = min(length(br.bytes), cld(br.maxbits - br.readbits, 8))
    return readbytes!(br.io, br.bytes, readcount)
end



"""
    readbits!(br::BitReader, ::Type{T}, n::Integer) where T<:Unsigned

Read `n` bits from `br` and return them as a `T`.

The bits read are placed in the leftmost `n` bits of the returned value. `n` must be no
more than the number of bits in `T`. An `EOFError` is thrown if less than `n` bits remain in
`br`.
"""
Base.@propagate_inbounds function readbits!(
    br::BitReader,
    ::Type{T},
    n::Integer,
) where T<:Unsigned
    @boundscheck if n > nbits(T)
        throw(ArgumentError("not enough bits in T"))
    end

    bits = zero(T)
    endshift = nbits(T) - n

    # Iterate until all requested bits have been read
    while n > 0
        # Make sure we haven't hit the read limit
        readcount = min(8 - (br.idx % 8), n)
        if br.readbits + readcount > br.maxbits
            throw(EOFError())
        end

        # Move the return bits over to make room for this iteration
        n -= readcount
        bits = bits << readcount

        # Which bits should be read this iteration?
        i = br.idx ÷ 8 + 1
        mask = 0xFF >> (br.idx % 8) # bits to remove from left side
        shift = 8 - readcount - (br.idx % 8) # bits to remove from right side
        newbits = (br.bytes[i] & mask) >> shift

        # Write newbits to bits
        bits |= newbits

        # Update indices
        br.idx += readcount
        br.readbits += readcount

        # Refill the buffer if we've reached the end
        if br.idx == length(br.bytes) * 8
            refillbits!(br)
        end
    end

    # Move the return value to the high side of bits and return it
    return bits << endshift
end



"""
    bitreader(f, io::IO; kw...)

Call `f` with a [`BitReader`](@ref) wrapping `io`.

`kw` are passed to the [`BitReader`](@ref) constructor.

If `maxbits` is specified, no more than that many bits are read from `io`.

# Examples

Read padded unsigned integers from a file of packed bits, adding one bit of padding on each
side of each element:

```julia
T = eltype(mypaddedgrid)
f = open("grid.dat")
bitreader(f) do br
    for i in eachindex(mypaddedgrid)
        mypaddedgrid = readbits!(bw, T, 8*sizeof(T)-2) >> 1
    end
end
```
"""
Base.@propagate_inbounds bitreader(f, io::IO; kw...) = f(BitReader(io; kw...))
