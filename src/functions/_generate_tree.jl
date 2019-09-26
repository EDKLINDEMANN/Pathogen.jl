function generate_tree(events::Events{T},
                       observations::EventObservations{T},
                       network::TransmissionNetwork) where {T <: SEIR}
  # Initialization
  trees = Tree[]
  eventtimes = [events.exposed observations.infected]
  eventorder = sortperm(eventtimes[:])

  exposure_tree = Dict{Int64, Int64}()
  exposure_internal = Dict{Int64, Int64}()
  observation_internal = Dict{Int64, Int64}()
  observation_leaf = Dict{Int64, Int64}()

  pathways = [pathwayfrom(i, network) for i = 1:events.individuals]

  # Determine significant exposures (results in a observation in the future)
  significant_exposures = [any(observations.infected[pathways[i]] .>= events.exposed[i]) for i = 1:events.individuals]

  # Iterate through all events to build tree
  for i = 1:length(eventorder)

    # Stop when no remaining valid event times
    isnan(eventtimes[eventorder[i]]) && break

    # Determine the individual and event type
    event_individual, event_type = ind2sub(size(eventtimes), eventorder[i])

    # For exposure events...
    if event_type == 1

      # Add exposure event to tree only if it is significant
      if significant_exposures[event_individual]

        # For significant external exposure events...
        if network.external[event_individual]

          # Add a new tree and a node
          tree = length(trees) + 1
          height = eventtimes[eventorder[i]]
          node = 1

          push!(trees, Tree(height))
          addnode!(trees[tree])

          exposure_tree[event_individual] = tree
          exposure_internal[event_individual] = node

        # For significant internal exposure events...
        else

          # Determine exposure source
          source = findfirst(network.internal[:, event_individual])

          # Previous significant exposures from this exposure source
          source_priorexposures = find(network.internal[source, :][:] .& (eventtimes[:, 1] .< eventtimes[eventorder[i]]) .& significant_exposures)

          # For when the exposure source has yet to be detected...
          if isnan(eventtimes[source, 2]) || eventtimes[source, 2] > eventtimes[event_individual, 1]

            # When there have been prior exposures from this undetected source
            if length(source_priorexposures) > 0
              tree = exposure_tree[source]
              parentnode = exposure_internal[source_priorexposures[indmax(eventtimes[source_priorexposures, 1])]]
              branch_length = eventtimes[eventorder[i]] - maximum(eventtimes[source_priorexposures, 1])

            # When there have not been any prior exposures from this undetected source
            else
              tree = exposure_tree[source]
              parentnode = exposure_internal[source]
              branch_length = eventtimes[eventorder[i]] - eventtimes[source, 1]
            end

          # For when the exposure source has been detected...
          else
            # And detection of exposure source is most recent, relevant event...
            if length(source_priorexposures) == 0 || all(eventtimes[source, 2] .> eventtimes[source_priorexposures, 1])
              tree = exposure_tree[source]
              parentnode = observation_internal[source]
              branch_length = eventtimes[eventorder[i]] - eventtimes[source, 2]

            # Or some prior exposure is most recent, relevant event...
            else
              tree = exposure_tree[source]
              parentnode = exposure_internal[source_priorexposures[indmax(eventtimes[source_priorexposures, 1])]]
              branch_length = eventtimes[eventorder[i]] - maximum(eventtimes[source_priorexposures, 1])
            end
          end
          # Add branch and node to tree
          branch!(trees[tree], parentnode, branch_length)

          # Record tree and node ID
          exposure_tree[event_individual] = tree
          exposure_internal[event_individual] = length(trees[tree].nodes)
        end
      end

    # Infection observation event
    elseif event_type == 2

      # Determine if individual has had significant exposures prior to being observed as infected
      priorexposures = find(network.internal[event_individual, :][:] .&
                            (eventtimes[:, 1] .< eventtimes[eventorder[i]]) .&
                            significant_exposures)
      if length(priorexposures) > 0
        tree = exposure_tree[event_individual]

        # Parent node is the internal node associated with the final pre-observation exposure
        parentnode = exposure_internal[priorexposures[indmax(eventtimes[priorexposures, 1])]]
        branch_length = eventtimes[eventorder[i]] - maximum(eventtimes[priorexposures, 1])

      # Individual has not exposed others prior to being observed as infected
      else
        tree = exposure_tree[event_individual]
        parentnode = exposure_internal[event_individual]
        branch_length = eventtimes[eventorder[i]] - eventtimes[event_individual, 1]
      end

      # Individual goes on to have significant exposures
      if any(network.internal[event_individual, :][:] .& (eventtimes[:, 1] .>= eventtimes[eventorder[i]]) .& significant_exposures)
        branch!(trees[tree], parentnode, branch_length)
        observation_internal[event_individual] = length(trees[tree].nodes)

        # Add zero-length branch so observation event is a leaf node
        branch!(trees[tree], observation_internal[event_individual], 0.)

      # Individual does not go one to have any significant exposures
      else
        branch!(trees[tree], parentnode, branch_length)
      end
      # Record node ID
      observation_leaf[event_individual] = length(trees[tree].nodes)
    end
  end
  return trees, exposure_tree, observation_leaf
end


function generate_tree(events::Events{T},
                       observations::EventObservations{T},
                       network::TransmissionNetwork) where {T <: SIR}
  # Initialization
  tree = Tree()
  eventtimes = [events.infected observations.infected]
  eventorder = sortperm(eventtimes[:])
  eventnodes = fill(Nullable{Int64}(), (events.individuals, 3))
  pathways = [pathwayfrom(i, network) for i = 1:events.individuals]

  # Determine significant infections (results in a observation in the future)
  significant_infections = [any(observations.infected[pathways[i]] .>= events.infected[i]) for i = 1:events.individuals]

  # Add a node for each infection observation event
  infection_observations = find(!isnan(observations.infected))
  addnodes!(tree, length(infection_observations))

  for i = 1:length(infection_observations)
    eventnodes[infection_observations[i], 3] = i
  end

  # Add a root node
  addnode!(tree)
  rootnode = length(tree.nodes)
  roottime = minimum(eventtimes)

  # Iterate through all events to build tree
  for i = 1:length(eventorder)

    # Stop when no remaining valid event times
    isnan(eventtimes[eventorder[i]]) && break

    # Determine the individual and event type
    event = ind2sub(size(eventtimes), eventorder[i])

    # For infection events...
    if event[2] == 1

      # Add infection event to tree only if it is significant
      if significant_infections[event[1]]

        # For significant external infection events...
        if network.external[event[1]]

          # Parent node of external infection is the root of the tree...
          parentnode = rootnode
          branch_length = eventtimes[eventorder[i]] - roottime

        # For significant internal infection events...
        else

          # Determine infection source
          source = findfirst(network.internal[:, event[1]])

          # Previous significant infections from this infection source
          source_priorinfections = network.internal[source, :][:] & (eventtimes[:, 1] .< eventtimes[eventorder[i]]) & significant_infections

          # For when the infection source has yet to be detected...
          if isnan(eventtimes[source, 2]) || eventtimes[source, 2] > eventtimes[event[1], 1]

            # When there have been prior infections from this undetected source
            if any(source_priorinfections)
              parentnode = get(eventnodes[source_priorinfections, 1][indmax(eventtimes[source_priorinfections, 1])])
              branch_length = eventtimes[eventorder[i]] - maximum(eventtimes[source_priorinfections, 1])

            # When there have not been any prior infections from this undetected source
            else
              parentnode = get(eventnodes[source, 1])
              branch_length = eventtimes[eventorder[i]] - eventtimes[source, 1]
            end

          # For when the infection source has been detected...
          else
            # And detection of infection source is most recent, relevant event...
            if !any(source_priorinfections) || all(eventtimes[source, 2] .> eventtimes[source_priorinfections, 1])
              parentnode = get(eventnodes[source, 2])
              branch_length = eventtimes[eventorder[i]] - eventtimes[source, 2]

            # Or some prior infection is most recent, relevant event...
            else
              parentnode = get(eventnodes[source_priorinfections, 1][indmax(eventtimes[source_priorinfections, 1])])
              branch_length = eventtimes[eventorder[i]] - maximum(eventtimes[source_priorinfections, 1])
            end
          end
        end

        # Add branch and node to tree
        branch!(tree, parentnode, branch_length)
        eventnodes[eventorder[i]] = length(tree.nodes)
      end

    # Infection observation event
    elseif event[2] == 2

      # Individual has had significant infections prior to being observed as infected
      priorinfections = network.internal[event[1], :][:] & (eventtimes[:, 1] .< eventtimes[eventorder[i]]) & significant_infections
      if any(priorinfections)
        parentnode = get(eventnodes[priorinfections, 1][indmax(eventtimes[priorinfections, 1])])
        branch_length = eventtimes[eventorder[i]] - maximum(eventtimes[priorinfections, 1])

      # Individual has not exposed others prior to being observed as infected
      else
        parentnode = get(eventnodes[event[1], 1])
        branch_length = eventtimes[eventorder[i]] - eventtimes[event[1], 1]
      end

      # Individual goes on to have significant infections
      if any(network.internal[event[1], :][:] & (eventtimes[:, 1] .>= eventtimes[eventorder[i]]) & significant_infections)
        branch!(tree, parentnode, branch_length)
        eventnodes[event[1], 2] = length(tree.nodes)

        # Add zero-length branch so observation event is a leaf node
        addbranch!(tree, get(eventnodes[event[1], 2]), get(eventnodes[event[1], 3]), 0.)

      # Individual does not go one to have any significant infections
      else
        addbranch!(tree, parentnode, get(eventnodes[event[1], 3]), branch_length)
      end
    end
  end
  return tree
end


function generate_tree(events::Events{T},
                       observations::EventObservations{T},
                       network::TransmissionNetwork) where {T <: SEI}
  # Initialization
  tree = Tree()
  eventtimes = [events.exposed observations.infected]
  eventorder = sortperm(eventtimes[:])
  eventnodes = fill(Nullable{Int64}(), (events.individuals, 3))
  pathways = [pathwayfrom(i, network) for i = 1:events.individuals]

  # Determine significant exposures (results in a observation in the future)
  significant_exposures = [any(observations.infected[pathways[i]] .>= events.exposed[i]) for i = 1:events.individuals]

  # Add a node for each infection observation event
  infection_observations = find(!isnan(observations.infected))
  addnodes!(tree, length(infection_observations))

  for i = 1:length(infection_observations)
    eventnodes[infection_observations[i], 3] = i
  end

  # Add a root node
  addnode!(tree)
  rootnode = length(tree.nodes)
  roottime = minimum(eventtimes)

  # Iterate through all events to build tree
  for i = 1:length(eventorder)

    # Stop when no remaining valid event times
    isnan(eventtimes[eventorder[i]]) && break

    # Determine the individual and event type
    event = ind2sub(size(eventtimes), eventorder[i])

    # For exposure events...
    if event[2] == 1

      # Add exposure event to tree only if it is significant
      if significant_exposures[event[1]]

        # For significant external exposure events...
        if network.external[event[1]]

          # Parent node of external exposure is the root of the tree...
          parentnode = rootnode
          branch_length = eventtimes[eventorder[i]] - roottime

        # For significant internal exposure events...
        else

          # Determine exposure source
          source = findfirst(network.internal[:, event[1]])

          # Previous significant exposures from this exposure source
          source_priorexposures = network.internal[source, :][:] & (eventtimes[:, 1] .< eventtimes[eventorder[i]]) & significant_exposures

          # For when the exposure source has yet to be detected...
          if isnan(eventtimes[source, 2]) || eventtimes[source, 2] > eventtimes[event[1], 1]

            # When there have been prior exposures from this undetected source
            if any(source_priorexposures)
              parentnode = get(eventnodes[source_priorexposures, 1][indmax(eventtimes[source_priorexposures, 1])])
              branch_length = eventtimes[eventorder[i]] - maximum(eventtimes[source_priorexposures, 1])

            # When there have not been any prior exposures from this undetected source
            else
              parentnode = get(eventnodes[source, 1])
              branch_length = eventtimes[eventorder[i]] - eventtimes[source, 1]
            end

          # For when the exposure source has been detected...
          else
            # And detection of exposure source is most recent, relevant event...
            if !any(source_priorexposures) || all(eventtimes[source, 2] .> eventtimes[source_priorexposures, 1])
              parentnode = get(eventnodes[source, 2])
              branch_length = eventtimes[eventorder[i]] - eventtimes[source, 2]

            # Or some prior exposure is most recent, relevant event...
            else
              parentnode = get(eventnodes[source_priorexposures, 1][indmax(eventtimes[source_priorexposures, 1])])
              branch_length = eventtimes[eventorder[i]] - maximum(eventtimes[source_priorexposures, 1])
            end
          end
        end

        # Add branch and node to tree
        branch!(tree, parentnode, branch_length)
        eventnodes[eventorder[i]] = length(tree.nodes)
      end

    # Infection observation event
    elseif event[2] == 2

      # Individual has had significant exposures prior to being observed as infected
      priorexposures = network.internal[event[1], :][:] & (eventtimes[:, 1] .< eventtimes[eventorder[i]]) & significant_exposures
      if any(priorexposures)
        parentnode = get(eventnodes[priorexposures, 1][indmax(eventtimes[priorexposures, 1])])
        branch_length = eventtimes[eventorder[i]] - maximum(eventtimes[priorexposures, 1])

      # Individual has not exposed others prior to being observed as infected
      else
        parentnode = get(eventnodes[event[1], 1])
        branch_length = eventtimes[eventorder[i]] - eventtimes[event[1], 1]
      end

      # Individual goes on to have significant exposures
      if any(network.internal[event[1], :][:] & (eventtimes[:, 1] .>= eventtimes[eventorder[i]]) & significant_exposures)
        branch!(tree, parentnode, branch_length)
        eventnodes[event[1], 2] = length(tree.nodes)

        # Add zero-length branch so observation event is a leaf node
        addbranch!(tree, get(eventnodes[event[1], 2]), get(eventnodes[event[1], 3]), 0.)

      # Individual does not go one to have any significant exposures
      else
        addbranch!(tree, parentnode, get(eventnodes[event[1], 3]), branch_length)
      end
    end
  end
  return tree
end


function generate_tree(events::Events{T},
                       observations::EventObservations{T},
                       network::TransmissionNetwork) where {T <: SI}
  # Initialization
  tree = Tree()
  eventtimes = [events.infected observations.infected]
  eventorder = sortperm(eventtimes[:])
  eventnodes = fill(Nullable{Int64}(), (events.individuals, 3))
  pathways = [pathwayfrom(i, network) for i = 1:events.individuals]

  # Determine significant infections (results in a observation in the future)
  significant_infections = [any(observations.infected[pathways[i]] .>= events.infected[i]) for i = 1:events.individuals]

  # Add a node for each infection observation event
  infection_observations = find(!isnan(observations.infected))
  addnodes!(tree, length(infection_observations))

  for i = 1:length(infection_observations)
    eventnodes[infection_observations[i], 3] = i
  end

  # Add a root node
  addnode!(tree)
  rootnode = length(tree.nodes)
  roottime = minimum(eventtimes)

  # Iterate through all events to build tree
  for i = 1:length(eventorder)

    # Stop when no remaining valid event times
    isnan(eventtimes[eventorder[i]]) && break

    # Determine the individual and event type
    event = ind2sub(size(eventtimes), eventorder[i])

    # For infection events...
    if event[2] == 1

      # Add infection event to tree only if it is significant
      if significant_infections[event[1]]

        # For significant external infection events...
        if network.external[event[1]]

          # Parent node of external infection is the root of the tree...
          parentnode = rootnode
          branch_length = eventtimes[eventorder[i]] - roottime

        # For significant internal infection events...
        else

          # Determine infection source
          source = findfirst(network.internal[:, event[1]])

          # Previous significant infections from this infection source
          source_priorinfections = network.internal[source, :][:] & (eventtimes[:, 1] .< eventtimes[eventorder[i]]) & significant_infections

          # For when the infection source has yet to be detected...
          if isnan(eventtimes[source, 2]) || eventtimes[source, 2] > eventtimes[event[1], 1]

            # When there have been prior infections from this undetected source
            if any(source_priorinfections)
              parentnode = get(eventnodes[source_priorinfections, 1][indmax(eventtimes[source_priorinfections, 1])])
              branch_length = eventtimes[eventorder[i]] - maximum(eventtimes[source_priorinfections, 1])

            # When there have not been any prior infections from this undetected source
            else
              parentnode = get(eventnodes[source, 1])
              branch_length = eventtimes[eventorder[i]] - eventtimes[source, 1]
            end

          # For when the infection source has been detected...
          else
            # And detection of infection source is most recent, relevant event...
            if !any(source_priorinfections) || all(eventtimes[source, 2] .> eventtimes[source_priorinfections, 1])
              parentnode = get(eventnodes[source, 2])
              branch_length = eventtimes[eventorder[i]] - eventtimes[source, 2]

            # Or some prior infection is most recent, relevant event...
            else
              parentnode = get(eventnodes[source_priorinfections, 1][indmax(eventtimes[source_priorinfections, 1])])
              branch_length = eventtimes[eventorder[i]] - maximum(eventtimes[source_priorinfections, 1])
            end
          end
        end

        # Add branch and node to tree
        branch!(tree, parentnode, branch_length)
        eventnodes[eventorder[i]] = length(tree.nodes)
      end

    # Infection observation event
    elseif event[2] == 2

      # Individual has had significant infections prior to being observed as infected
      priorinfections = network.internal[event[1], :][:] & (eventtimes[:, 1] .< eventtimes[eventorder[i]]) & significant_infections
      if any(priorinfections)
        parentnode = get(eventnodes[priorinfections, 1][indmax(eventtimes[priorinfections, 1])])
        branch_length = eventtimes[eventorder[i]] - maximum(eventtimes[priorinfections, 1])

      # Individual has not exposed others prior to being observed as infected
      else
        parentnode = get(eventnodes[event[1], 1])
        branch_length = eventtimes[eventorder[i]] - eventtimes[event[1], 1]
      end

      # Individual goes on to have significant infections
      if any(network.internal[event[1], :][:] & (eventtimes[:, 1] .>= eventtimes[eventorder[i]]) & significant_infections)
        branch!(tree, parentnode, branch_length)
        eventnodes[event[1], 2] = length(tree.nodes)

        # Add zero-length branch so observation event is a leaf node
        addbranch!(tree, get(eventnodes[event[1], 2]), get(eventnodes[event[1], 3]), 0.)

      # Individual does not go one to have any significant infections
      else
        addbranch!(tree, parentnode, get(eventnodes[event[1], 3]), branch_length)
      end
    end
  end
  return tree
end
