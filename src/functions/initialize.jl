function initialize(::Type{TransmissionRates},
                    states::Vector{DiseaseState},
                    pop::Population,
                    rf::RiskFunctions{M},
                    rp::RiskParameters{M}) where {
                    S <: DiseaseStateSequence,
                    M <: ILM{S}}
  n_ids = length(states)
  tr = TransmissionRates(n_ids)
  for i in findall(states .== Ref(State_S))
    # External exposure
    #tr.external[i] = rf.susceptibility(rp.susceptibility, pop, i) * rf.sparks(rp.sparks, pop, i)
    tr.external[i] = rf.sparks(rp.sparks, pop, i)
    # Internal exposure
    for k in findall(states .== Ref(State_I))
      tr.internal[k, i] = rf.susceptibility(rp.susceptibility, pop, i) *
                          rf.infectivity(rp.infectivity, pop, i, k) *
                          rf.transmissibility(rp.transmissibility, pop, k)
    end
  end
  @debug "Initialization of $M TransmissionRates complete" external = tr.external ∑external = sum(tr.external) internal = tr.internal ∑internal = sum(tr.internal)
  return tr
end

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

function initialize(::Type{MarkovChain},
                    mcmc::MCMC{M},
                    progress_channel::RemoteChannel;
                    attempts::Int64=1000) where {
                    S <: DiseaseStateSequence,
                    M <: ILM{S}}
  if attempts <= 0
    @error "Must have at least 1 initialization attempt"
  end
  max_lposterior = -Inf
  local markov_chain
  for i in 1:attempts
    events = generate(Events, mcmc)
    rparams = generate(RiskParameters, mcmc)
    lprior = logpriors(rparams, mcmc.risk_priors)
    llikelihood, network = loglikelihood(rparams,
                                         mcmc.risk_functions,
                                         events,
                                         mcmc.population,
                                         mcmc.starting_states,
                                         early_decision_value = max_lposterior - lprior)
    if M <: PhyloILM
      trees, tree_id, obs_leaf_node_id = generate_tree(events, mcmc.event_observations, network)
      sm = generate(mcmc.substitution_model, substitution_model_priors)
      for i = eachindex(trees)
        tree_id
        loglikelihood(trees[i], sm, node_data)



  trees, tree_id, obs_leaf_node_id = generate_tree(sim.events,
                                                   infection,
                                                   sim.transmission_network)
  seq_full = [simulate(RNASeq, t, sim.substitution_model, seq_len) for t in trees]
  seq_obs = [isnan(infection[i])? nothing | seq_full[tree_id[i]][observation_leaf_node_id] for i = eachindex(infection)]


    end

    lposterior = llikelihood + lprior

    if lposterior > max_lposterior
      markov_chain = MarkovChain(events, network, rparams, lposterior)
      max_lposterior = lposterior
    end
    put!(progress_channel, true)
  end
  if max_lposterior == -Inf
    @error "Failed to initialize Markov Chain"
  end
  return markov_chain
end

function initialize(::Type{MarkovChain},
                    mcmc::MCMC{M};
                    attempts::Int64=1000) where {
                    S <: DiseaseStateSequence,
                    M <: ILM{S}}
  if attempts <= 0
    @error "Must have at least 1 initialization attempt"
  end
  max_lposterior = -Inf
  local markov_chain
  pmeter = Progress(attempts, "Initialization progress")
  for i in 1:attempts
    @debug "Beginning MarkovChain initialization attempt $i"
    next!(pmeter)
    events = generate(Events, mcmc)
    rparams = generate(RiskParameters, mcmc)
    lprior = logpriors(rparams, mcmc.risk_priors)
    llikelihood, network = loglikelihood(rparams,
                                         mcmc.risk_functions,
                                         events,
                                         mcmc.population,
                                         mcmc.starting_states,
                                         early_decision_value = max_lposterior - lprior)
    lposterior = llikelihood + lprior
    if lposterior > max_lposterior
      markov_chain = MarkovChain(events, network, rparams, lposterior)
      max_lposterior = lposterior
    end
  end
  finish!(pmeter)
  if max_lposterior == -Inf
    @error "Failed to initialize Markov Chain"
  end
  return markov_chain
end
