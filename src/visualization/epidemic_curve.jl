function _count_by_state(events::Events{M},
                         state::DiseaseState,
                         time::Float64) where M <: ILM
  if time < 0.0
    @error "Time must be ≥ 0.0"
  elseif state ∉ convert(DiseaseStates, S)
    @error "Invalid state specified"
  end
  n_ids = 0
  if state == convert(DiseaseStates, S)[1] # S
    nextstate = advance(state, T) # Either E or I
    n_ids += sum(events[nextstate] .> Ref(time)) # E/I after `time`
    n_ids += sum(isnan.(events[nextstate])) # Never E/I
  elseif state in convert(DiseaseStates, S)[2:end-1]
    nextstate = advance(state, T) # Either I or R
    n_ids += sum((events[state] .<= Ref(time)) .& (events[nextstate] .> Ref(time))) # E/I at or before time and I/R after time
    n_ids += sum((events[state] .<= Ref(time)) .& isnan.(events[nextstate])) # E/I at or before time and never I/R
  elseif state == convert(DiseaseStates, S)[end] # I or R
    n_ids += sum(events[state] .<= Ref(time)) # I/R at or before time
  end
  @debug "$n_ids individual(s) in state $state at t = $time"
  return n_ids
end

function _epidemic_curve(events::Events{M},
                         state::DiseaseState,
                         min::Float64,
                         max::Float64) where M <: ILM
  if min >= max
    @error "Minimum time must be less than maximum time"
  end
  local times
  if state == convert(DiseaseStates, S)[1]
    nextstate = convert(DiseaseStates, S)[2]
    times = events[nextstate]
  elseif state in convert(DiseaseStates, S)[2:end-1]
    nextstate = advance(state, T)
    times = events[[state; nextstate]][:]
  elseif state == convert(DiseaseStates, S)[end]
    times = events[state]
  else
    @error "Invalid state specified"
  end
  times = times[Ref(min) .< times .< Ref(max)]
  sort!(times)
  insert!(times, 1, min)
  push!(times, max)
  counts = [_count_by_state(events, state, t) for t in times]
  return times, counts
end

@recipe function f(events::Events{M},
                   state::DiseaseState,
                   min::Float64,
                   max::Float64) where T<: EpidemicModel
  xguide --> "Time"
  yguide --> "N"
  xlims --> (min - 1.0, max + 1.0)
  linewidth --> 2.0
  linecolor --> :cornflowerblue
  label --> ""
  seriestype --> :steppost
  _epidemic_curve(events, state, min, max)
end

@recipe function f(events::Events{M},
                   state::DiseaseState) where M <: ILM
  events, state, 0.0, maximum(events)
end

@recipe function f(events::Events{M},
                   min::Float64,
                   max::Float64) where T<: EpidemicModel
  @series begin
    linecolor --> :purple
    label --> "S"
    events, State_S, min, max
  end
  if S in [SEIR; SEI]
    @series begin
      linecolor --> :lightblue4
      label --> "E"
      events, State_E, min, max
    end
  end
  @series begin
    linecolor --> :lightgreen
    label --> "I"
    events, State_I, min, max
  end
  if S in [SEIR; SIR]
    @series begin
      linecolor --> :yellow
      label --> "R"
      events, State_R, min, max
    end
  end
end

@recipe function f(events::Events{M}) where M <: ILM
  events, 0.0, maximum(events)
end
