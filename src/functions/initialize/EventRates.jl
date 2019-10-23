function initialize(::Type{EventRates},
                    tr::TransmissionRates,
                    states::Vector{DiseaseState},
                    pop::Population,
                    rf::RiskFunctions{M},
                    rp::RiskParameters{M}) where {
                    S <: DiseaseStateSequence,
                    M <: ILM{S}}
  n_ids = length(states)
  rates = EventRates{M}(n_ids)
  for i = 1:n_ids
    if states[i] == State_S
      if S in [SEIR; SEI]
        rates.exposure[i] = tr.external[i] + sum(tr.internal[:,i])
      elseif S in [SIR; SI]
        rates.infection[i] = tr.external[i] + sum(tr.internal[:,i])
      end
    elseif states[i] == State_E
      rates.infection[i] = rf.latency(rp.latency, pop, i)
    elseif states[i] == State_I
      if S in [SEIR; SIR]
        rates.removal[i] = rf.removal(rp.removal, pop, i)
      end
    end
  end
  @debug "Initialization of $M EventRates complete" rates = rates[_state_progressions[S][2:end]]
  return rates
end
