struct Event{T <: ILM}
  time::Float64
  individual::Union{Nothing, Int64}
  new_state::Union{Nothing, DiseaseState}

  function Event{M}(time::Float64) where M <: ILM
    if time !== Inf
      throw(ErrorException("This initialization method is intended for generated event times of Inf, i.e. when an epidemic is complete and no further events are possible."))
    end
    return new{T}(Inf, nothing, nothing)
  end

  function Event{M}(time::Float64, id::Int64, new_state::DiseaseState) where M <: ILM
    if !(new_state in _state_progressions[T])
      throw(ErrorException("Invalid disease state provided for $(T) Epidemic Models.")
    end
    return new{T}(time, id, new_state)
  end

  function Event(new_time::Float64, event::Event{M}) where M <: ILM
    return new{T}(new_time, event.individual, event.new_state)
  end
end

function Base.copy(x::Event{M}) where M <: ILM
  return Event{M}(copy(x.time), copy(x.id), copy(x.new_state))
end

function _new_transmission(e::Event{M}) where {S <: Union{SEIR, SEI}, M <: ILM{S}}
  return e.new_state == State_E
end

function _new_transmission(e::Event{M}) where {S <: Union{SIR, SI}, M <: ILM{S}}
  return e.new_state == State_I
end

function _time(x::Event{M}) where M <: ILM
  return x.time
end

function _individual(x::Event{M}) where M <: ILM
  return x.individual
end

function _new_state(x::Event{M}) where M <: ILM
  return x.new_state
end

function Base.show(io::IO, x::Event{M}) where M <: ILM
  return print(io, "Transition of individual $(x.individual) into $(x.new_state) state at t = $(x.time)")
end
