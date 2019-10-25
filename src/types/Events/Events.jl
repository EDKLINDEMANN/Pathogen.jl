struct Events{M <: ILM}
  exposure::Union{Nothing, Vector{Float64}}
  infection::Union{Nothing, Vector{Float64}}
  removal::Union{Nothing, Vector{Float64}}
  individuals::Int64
  start_time::Float64

  function Events{M}(e::V,
                     i::V,
                     r::V;
                     start_time::F=0.0) where {
                     S <: SEIR,
                     M <: ILM{S},
                     F <: Float64
                     V <: Vector{F}}
    if length(unique((length.([e; i; r])))) != 1
      throw(ErrorException("Length of event time vectors must be equal"))
    end
    return new{M}(e, i, r, length(i), start_time)
  end

  function Events{M}(n::Int64;
                     start_time::Float64=0.0) where {
                     S <: SEIR,
                     M <: ILM{S}}
    return new{M}(fill(NaN, n), fill(NaN, n), fill(NaN, n), n, start_time)
  end

  function Events{M}(e::V,
                     i::V;
                     start_time::F=0.0) where {
                     S <: SEI,
                     M <: ILM{S},
                     F <: Float64
                     V <: Vector{F}}
    if length(unique((length.([e; i])))) != 1
      throw(ErrorException("Length of event time vectors must be equal"))
    end
    return new{M}(e, i, nothing, length(i), start_time)
  end

  function Events{M}(n::Int64;
                     start_time::Float64=0.0) where {
                     S <: SEI,
                     M <: ILM{S}}
    return new{M}(fill(NaN, n), fill(NaN, n), nothing, n, start_time)
  end

  function Events{M}(i::V,
                     r::V;
                     start_time::F=0.0) where {
                     S <: SIR,
                     M <: ILM{S},
                     F <: Float64
                     V <: Vector{F}}
    if length(unique((length.([i; r])))) != 1
      throw(ErrorException("Length of event time vectors must be equal"))
    end
    return new{M}(nothing, i, r, length(i), start_time)
  end

  function Events{M}(n::Int64;
                     start_time::Float64=0.0) where {
                     S <: SIR,
                     M <: ILM{S}}
    return new{M}(nothing, fill(NaN, n), fill(NaN, n), n, start_time)
  end

  function Events{M}(i::V;
                     start_time::F=0.0) where {
                     S <: SI,
                     M <: ILM{S},
                     F <: Float64
                     V <: Vector{F}}
    if length(unique((length.([i; r])))) != 1
      throw(ErrorException("Length of event time vectors must be equal"))
    end
    return new{M}(nothing, i, nothing, length(i), start_time)
  end

  function Events{M}(n::Int64;
                     start_time::Float64=0.0) where {
                     S <: SI,
                     M <: ILM{S}}
    return new{M}(nothing, fill(NaN, n), nothing, n, start_time)
  end
end

function Base.show(io::IO, x::Events{M}) where M <: ILM
  return print(io, "$T event times (n=$(x.individuals))")
end

function Base.getindex(x::Events{M},
                       state::DiseaseState) where {
                       S <: DiseaseStateSequence,
                       M <: ILM{S}}
  if state == State_E
    return x.exposure
  elseif state == State_I
    return x.infection
  elseif state == State_R
    return x.removal
  else
    @error "Unrecognized indexing disease state"
  end
end

function Base.getindex(x::Events{M},
                       states::Vector{DiseaseState}) where {
                       S <: DiseaseStateSequence,
                       M <: ILM{S}}
  y = x[states[1]]
  for i = 2:length(states)
    y = hcat(y, x[states[i]])
  end
  return y
end

function Events{M}(a::Array{Float64,2};
                   start_time::Float64=0.0) where {
                   S <: DiseaseStateSequence,
                   M <: ILM{S}}
  if size(a, 2) == 4
    return Events{M}(a[:,1], a[:,2], a[:,3], a[:,4], start_time = start_time)
  elseif size(a, 2) == 3
    return Events{M}(a[:,1], a[:,2], a[:,3], start_time = start_time)
  elseif size(a, 2) == 2
    return Events{M}(a[:,1], a[:,2], start_time = start_time)
  else
    throw(ErrorException("Invalid array size for construction of $(Events{T})"))
  end
end

function Events{M}(x::Vector{DiseaseState};
                   start_time::Float64=0.0) where {
                   S <: DiseaseStateSequence,
                   M <: ILM{S}}
  events = Events{T}(length(x), start_time = start_time)
  for i = 1:length(x)
    if x[i] in _state_progressions[S]
      for j = _state_progressions[S][2:findfirst(Ref(x[i]) .== _state_progressions[S])]
        events[j][i] = -Inf
      end
    else
      @error "Invalid initial state for individual $i"
    end
  end
  return events
end

function Base.convert(::Type{Array{Float64, 2}},
                      x::Events{M}) where {
                      S <: DiseaseStateSequence,
                      M <: ILM{S}}
  return x[_state_progressions[S][2:end]]
end

function Base.convert(::Type{Vector{Float64}},
                      x::Events{M}) where {
                      S <: DiseaseStateSequence,
                      M <: ILM{S}}
  return x[_state_progressions[S][2:end]][:]
end

function Base.convert(::Type{Array{Float64, 2}},
                      x::Vector{Events{M}}) where {
                      M <: ILM}
  y = convert(Vector{Float64}, x[1])'
  for i = 2:length(x)
    y = vcat(y, convert(Vector{Float64}, x[i])')
  end
  return y
end

function Base.minimum(x::Events{M}) where M <: ILM
  y = convert(Array{Float64, 2}, x)
  return minimum(y[.!isnan.(y)])
end

function Base.maximum(x::Events{M}) where M <: ILM
  y = convert(Array{Float64, 2}, x)
  return maximum(y[.!isnan.(y)])
end

function Statistics.mean(x::Vector{Events{M}}) where M <: ILM
  return Events{M}(mean([convert(Array{Float64, 2}, i) for i in x]))
end
