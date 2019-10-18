function observe(sim::Simulation{M},
                 delay_infection::UnivariateDistribution,
                 delay_removal::UnivariateDistribution;
                 force::Bool = false) where {S <: Union{SEIR, SIR}, M <: TNILM{S}}
  infection = fill(NaN, sim.events.individuals)
  removal = fill(NaN, sim.events.individuals)
  if force
    for i in findall(.!isnan.(sim.events.infection))
      if sim.events.infection[i] == -Inf
        infection[i] = NaN
        if sim.events.removal[i] == -Inf
          removal[i] = NaN
        end
      elseif isnan(sim.events.removal[i])
        infection[i] = sim.events.infection[i] + rand(delay_infection)
      elseif sim.events.removal[i] > -Inf
        infection_delay_ub = sim.events.removal[i] - sim.events.infection[i]
        infection[i] = sim.events.infection[i] +
                       rand(Truncated(delay_infection, 0.0, infection_delay_ub))
        removal[i] = sim.events.removal[i] + rand(delay_removal)
      end
      @debug "Infection observation of i = $i at t = $(round(infection[i], digits=3)) (actual infection onset at t = $(round(sim.events.infection[i], digits=3)))"
      @debug "Removal observation of i = $i at t = $(round(removal[i], digits=3)) (actual removal at t = $(round(sim.events.removal[i], digits=3)))"
    end
  else
    for i in findall(.!isnan.(sim.events.infection))
      if sim.events.infection[i] == -Inf
        infection[i] = NaN
        if sim.events.removal == -Inf
          removal[i] = NaN
        end
      else
        infection_delay = rand(delay_infection)
        if isnan(sim.events.removal[i])
          infection[i] = sim.events.infection[i] + infection_delay
        else
          if infection_delay + sim.events.infection[i] < sim.events.removal[i]
            infection[i] = sim.events.infection[i] + infection_delay
            removal[i] = sim.events.removal[i] + rand(delay_removal)
          end
        end
      end
      @debug "Infection observation of i = $i t = $(round(infection[i], digits=3)) (actual infection onset at t = $(round(sim.events.infection[i], digits=3)))"
      @debug "Removal observation of i = $i at t = $(round(removal[i], digits=3)) (actual removal at t = $(round(sim.events.removal[i], digits=3)))"
    end
  end
  return EventObservations{M}(infection, removal)
end


function observe(sim::Simulation{M},
                 delay_infection::UnivariateDistribution) where {
                 S <: Union{SEI, SI},
                 M <: TNILM{S}}
  infection = fill(NaN, sim.events.individuals)
  @simd for i in findall(.!isnan.(sim.events.infection))
    if sim.events.infection[i] == -Inf
      infection[i] = -Inf
    else
      infection[i] = sim.events.infection[i] + rand(delay_infection)
    end
    @debug "Infection observation of i = $i at t = $(round(infection[i], digits=3))) (actual infection onset at t = $(round(sim.events.infection[i], digits=3)))"
  end
  return EventObservations{M}(infection)
end


function observe(sim::Simulation{M},
                 delay_infection::UnivariateDistribution,
                 delay_removal::UnivariateDistribution;
                 force::Bool = false) where {S <: Union{SEIR, SIR}, M <: PhyloILM{S}}
  infection = fill(NaN, sim.events.individuals)
  removal = fill(NaN, sim.events.individuals)
  if force
    for i in findall(.!isnan.(sim.events.infection))
      if sim.events.infection[i] == -Inf
        infection[i] = NaN
        if sim.events.removal[i] == -Inf
          removal[i] = NaN
        end
      elseif isnan(sim.events.removal[i])
        infection[i] = sim.events.infection[i] + rand(delay_infection)
      elseif sim.events.removal[i] > -Inf
        infection_delay_ub = sim.events.removal[i] - sim.events.infection[i]
        infection[i] = sim.events.infection[i] +
                       rand(Truncated(delay_infection, 0.0, infection_delay_ub))
        removal[i] = sim.events.removal[i] + rand(delay_removal)
      end
      @debug "Infection observation of i = $i at t = $(round(infection[i], digits=3)) (actual infection onset at t = $(round(sim.events.infection[i], digits=3)))"
      @debug "Removal observation of i = $i at t = $(round(removal[i], digits=3)) (actual removal at t = $(round(sim.events.removal[i], digits=3)))"
    end
  else
    for i in findall(.!isnan.(sim.events.infection))
      if sim.events.infection[i] == -Inf
        infection[i] = NaN
        if sim.events.removal == -Inf
          removal[i] = NaN
        end
      else
        infection_delay = rand(delay_infection)
        if isnan(sim.events.removal[i])
          infection[i] = sim.events.infection[i] + infection_delay
        else
          if infection_delay + sim.events.infection[i] < sim.events.removal[i]
            infection[i] = sim.events.infection[i] + infection_delay
            removal[i] = sim.events.removal[i] + rand(delay_removal)
          end
        end
      end
      @debug "Infection observation of i = $i t = $(round(infection[i], digits=3)) (actual infection onset at t = $(round(sim.events.infection[i], digits=3)))"
      @debug "Removal observation of i = $i at t = $(round(removal[i], digits=3)) (actual removal at t = $(round(sim.events.removal[i], digits=3)))"
    end
  end
  tree, event_nodes = generate(PhyloTree,
                               sim.events,
                               infection,
                               sim.transmission_network)
  seq_full = simulate(RNASeq, tree, sim.substitution_model, seq_len)
  seq_obs = [isnothing(event_nodes[i, 2])? nothing | seq_full[event_nodes[i, 2]] for i = eachindex(infection)]
  return EventObservations{M}(infection, removal, seq_obs)
end


function observe(sim::Simulation{M},
                 delay_infection::UnivariateDistribution,
                 seq_len::Int64) where {
                 S <: Union{SEI, SI},
                 M <: PhyloILM{S}}
  infection = fill(NaN, sim.events.individuals)
  @simd for i in findall(.!isnan.(sim.events.infection))
    if sim.events.infection[i] == -Inf
      infection[i] = -Inf
    else
      infection[i] = sim.events.infection[i] + rand(delay_infection)
    end
    @debug "Infection observation of i = $i at t = $(round(infection[i], digits=3))) (actual infection onset at t = $(round(sim.events.infection[i], digits=3)))"
  end
  trees, tree_id, obs_leaf_node_id = generate_tree(sim.events,
                                                   infection,
                                                   sim.transmission_network)
  seq_full = [simulate(RNASeq, t, sim.substitution_model, seq_len) for t in trees]
  seq_obs = [isnan(infection[i])? nothing | seq_full[tree_id[i]][observation_leaf_node_id] for i = eachindex(infection)]
  return EventObservations{M}(infection, seq_obs)
end
