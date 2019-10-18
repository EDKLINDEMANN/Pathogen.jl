struct EventExtents{M <: ILM}
  exposure::Union{Nothing, Float64}
  infection::Union{Nothing, Float64}
  removal::Union{Nothing, Float64}

  function EventExtents{M}(ex::F,
                           in::F,
                           re::F) where {
                           F <: Float64,
                           S <: SEIR,
                           M <: ILM{S}}
    return new{M}(ex, in, re)
  end

  function EventExtents{M}(ex::F,
                           in::F) where {
                           F <: Float64,
                           S <: SEI,
                           M <: ILM{S}}
    return new{M}(ex, in, nothing)
  end

  function EventExtents{M}(in::F,
                           re::F) where {
                           F <: Float64,
                           S <: SIR,
                           M <: ILM{S}}
    return new{M}(nothing, in, re)
  end

  function EventExtents{M}(in::F) where {
                           F <: Float64,
                           S <: SI,
                           M <: ILM{S}}
    return new{M}(nothing, in, nothing)
  end
end

function Base.show(io::IO, x::EventExtents{}) where {M <: ILM}
  print(io, "$M event extents")
end
