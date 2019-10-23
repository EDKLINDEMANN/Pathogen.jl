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
    lprior = logprior(rparams, mcmc.risk_priors)
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
    lprior = logprior(rparams, mcmc.risk_priors)
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
