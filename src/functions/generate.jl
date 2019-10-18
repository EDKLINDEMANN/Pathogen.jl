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

function generate(::Type{Transmission},
                  tr::TransmissionRates,
                  event::Event{M}) where {
                  S <: DiseaseStateSequence,
                  M <: ILM{S}}
  id = event.individual
  if _new_transmission(event)
    external_or_internal = Weights([tr.external[id]; sum(tr.internal[:,id])])
    if sum(external_or_internal) == 0.0
      @error "All transmission rates = 0.0, No transmission can be generated"
      return NoTransmission()
      # return ExogenousTransmission(id)
    elseif sample([true; false], external_or_internal)
      @debug "Exogenous tranmission generated"
      return ExogenousTransmission(id)
    else
      source = sample(1:tr.individuals, Weights(tr.internal[:, id]))
      @debug "Endogenous transmission generated (source id = $source)"
      return EndogenousTransmission(id, source)
    end
  else
    @debug "No transmission generated"
    return NoTransmission()
  end
end

function generate(::Type{Events},
                  obs::EventObservations{M},
                  extents::EventExtents{M}) where {
                  S <: DiseaseStateSequence,
                  M <: ILM{S}}
  events = Events{M}(obs.individuals)
  exposed_state = T in [SEIR; SEI]
  removed_state = T in [SEIR; SIR]
  for i = 1:obs.individuals
    if obs.infection[i] == -Inf
      update!(events, Event{M}(-Inf, i, State_I))
      if exposed_state
        update!(events, Event{M}(-Inf, i, State_E))
      end
      if removed_state && obs.removal[i] == -Inf
        update!(events, Event{M}(-Inf, i, State_R))
      end
    elseif !isnan(obs.infection[i])
      i_lb = obs.infection[i] - extents.infection
      i_ub = obs.infection[i]
      i_time = rand(Uniform(i_lb, i_ub))
      update!(events, Event{M}(i_time, i, State_I))
      if exposed_state
        e_lb = i_time - extents.exposure
        e_ub = i_time
        e_time = rand(Uniform(e_lb, e_ub))
        update!(events, Event{M}(e_time, i, State_E))
      end
      if removed_state && !isnan(obs.removal[i])
        r_lb = maximum([obs.infection[i]; obs.removal[i] - extents.removal])
        r_ub = obs.removal[i]
        r_time = rand(Uniform(r_lb, r_ub))
        update!(events, Event{M}(r_time, i, State_R))
      end
    end
  end
  return events
end

function generate(::Type{Events},
                  mcmc::MCMC{M}) where {
                  S <: DiseaseStateSequence,
                  M <: ILM{S}}
  return generate(Events, mcmc.event_observations, mcmc.event_extents)
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

function generate(::Type{RiskParameters{M}},
                  rpriors::RiskPriors{M}) where {
                  S <: DiseaseStateSequence,
                  M <: ILM{S}}
  sparks = Float64[rand(x) for x in rpriors.sparks]
  susceptibility = Float64[rand(x) for x in rpriors.susceptibility]
  infectivity = Float64[rand(x) for x in rpriors.infectivity]
  transmissibility = Float64[rand(x) for x in rpriors.transmissibility]
  if T in [SEIR; SEI]
    latency = Float64[rand(x) for x in rpriors.latency]
  end
  if T in [SEIR; SIR]
    removal = Float64[rand(x) for x in rpriors.removal]
  end
  if T == SEIR
    return RiskParameters{M}(sparks,
                             susceptibility,
                             infectivity,
                             transmissibility,
                             latency,
                             removal)
  elseif T == SEI
    return RiskParameters{M}(sparks,
                             susceptibility,
                             infectivity,
                             transmissibility,
                             latency)
  elseif T == SIR
    return RiskParameters{M}(sparks,
                             susceptibility,
                             infectivity,
                             transmissibility,
                             removal)
  elseif T == SI
    return RiskParameters{M}(sparks,
                             susceptibility,
                             infectivity,
                             transmissibility)
  end
end

function generate(::Type{RiskParameters},
                  mcmc::MCMC{M}) where {
                  S <: DiseaseStateSequence,
                  M <: ILM{S}}
  return generate(RiskParameters{M}, mcmc.risk_priors)
end

function generate(::Type{RiskParameters{M}},
                  last_rparams::RiskParameters{M},
                  Σ::Array{Float64, 2}) where {
                  S <: DiseaseStateSequence,
                  M <: ILM{S}}
  rparams_vector = rand(MvNormal(convert(Vector{Float64}, last_rparams), Σ))
  return _like(last_rparams, rparams_vector)
end

function generate(::Type{RiskParameters{M}},
                  mc::MarkovChain{M},
                  Σ::Array{Float64, 2}) where {
                  S <: DiseaseStateSequence,
                  M <: ILM{S}}
  return generate(RiskParameters{M}, mc.risk_parameters[end], Σ)
end
