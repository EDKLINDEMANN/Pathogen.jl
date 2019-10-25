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
      sm = generate(mcmc.substitution_model, substitution_model_priors)
      lprior += logprior(sm, substitution_model_priors)
      tree, obs_nodes = generate(events, mcmc.event_observations, network)
      leaf_data = Dict{Int64, GeneticSeq}()
      for k = eachindex(obs_nodes)
        if !isnothing(obs_nodes[k])
          leaf_data[obs_nodes[k]] = mcmc.event_observations.seq[k]
        end
      end
      llikelihood += loglikelihood(tree, sm, leaf_data)
    end
    lposterior = llikelihood + lprior

    if lposterior > max_lposterior
      if M <: TNILM
        markov_chain = MarkovChain(events, network, rparams, sm, lposterior)
      elseif M <: PhyloILM
        markov_chain = MarkovChain(events, network, rparams, lposterior)
      end
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
    if M <: PhyloILM
      sm = generate(mcmc.substitution_model, substitution_model_priors)
      lprior += logprior(sm, substitution_model_priors)
      tree, obs_nodes = generate(events, mcmc.event_observations, network)
      leaf_data = Dict{Int64, GeneticSeq}()
      for k = eachindex(obs_nodes)
        if !isnothing(obs_nodes[k])
          leaf_data[obs_nodes[k]] = mcmc.event_observations.seq[k]
        end
      end
      llikelihood += loglikelihood(tree, sm, leaf_data)
    end
    lposterior = llikelihood + lprior
    if lposterior > max_lposterior
      if M <: TNILM
        markov_chain = MarkovChain(events, network, rparams, sm, lposterior)
      elseif M <: PhyloILM
        markov_chain = MarkovChain(events, network, rparams, lposterior)
      end
      max_lposterior = lposterior
    end
  end
  finish!(pmeter)
  if max_lposterior == -Inf
    @error "Failed to initialize Markov Chain"
  end
  return markov_chain
end
