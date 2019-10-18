function generate(::Type{Event},
                  rates::EventRates{M},
                  time::Float64) where {
                  S <: DiseaseStateSequence,
                  M <: ILM{S}}
  totals = Weights([sum(rates[state]) for state in _state_progressions[T][2:end]])
  if sum(totals) == Inf
    new_state = sample(_state_progressions[T][2:end], totals)
    id = sample(1:rates.individuals, Weights(rates[new_state]))
    return Event{M}(time, id, new_state)
  elseif sum(totals) == 0.0
    return Event{M}(Inf)
  else
    # Generate new state
    new_state = sample(_state_progressions[T][2:end], totals)
    # Generate event individual
    id = sample(1:rates.individuals, Weights(rates[new_state]))
    return Event{M}(time + rand(Exponential(1.0 / sum(totals))), id, new_state)
  end
end

function generate(::Type{Event},
                  last_event::Event{M},
                  σ::Float64,
                  extents::EventExtents{M},
                  obs::EventObservations,
                  events::Events{M}) where {
                  S <: DiseaseStateSequence,
                  M <: ILM{S}}
  lowerbound, upperbound = _bounds(last_event, extents, obs, events)
  time = rand(TruncatedNormal(last_event.time,
                              σ,
                              lowerbound,
                              upperbound))
  return Event(time, last_event)
end

function generate(::Type{Event},
                  last_event::Event{M},
                  σ::Float64,
                  extents::EventExtents{M},
                  obs::EventObservations,
                  events::Events{M},
                  network::TransmissionNetwork) where {
                  S <: DiseaseStateSequence,
                  M <: ILM{S}}
  lowerbound, upperbound = _bounds(last_event, extents, obs, events, network)
  time = rand(TruncatedNormal(last_event.time,
                              σ,
                              lowerbound,
                              upperbound))
  return Event(time, last_event)
end
