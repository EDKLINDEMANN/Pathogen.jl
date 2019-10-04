function generate_tree(events::Events{T},
                       observations::EventObservations{T},
                       network::TransmissionNetwork) where {T <: EpidemicModel}
  # Initialization
  trees = Tree[]
  local event_times

  if T in [SEIR; SEI]
    event_times = [events.exposure observations.infection]
  elseif T in [SIR; SI]
    event_times = [events.infection observations.infection]
  end

  event_order = sortperm(event_times[:])
  event_lookup = CartesianIndices(size(event_times))

  tree_id = Dict{Int64, Int64}()
  transmission_node_id = Dict{Int64, Int64}()
  observation_internal_node_id = Dict{Int64, Int64}()
  observation_leaf_node_id = Dict{Int64, Int64}()

  pathways = [_pathway_from(i, network) for i = 1:events.individuals]

  # Determine significant transmissions (results in a observation in the future)
  significant_transmissions = [event_times[i, 1] == -Inf || any(event_times[pathways[i], 2] .>= Ref(event_times[i, 1])) for i = 1:events.individuals]

  # Iterate through all events to build tree
  for i = 1:length(event_order)

    # Stop when no remaining valid event times
    isnan(event_times[event_order[i]]) && break

    # Determine the individual and event type
    event_individual, event_type = Tuple(event_lookup[event_order[i]])

    @debug "Building phylogenetic tree..." Event = i Individual = event_individual Type = event_type == 1 ? "Transmission" : "Observation"

    # For transmission events...
    if event_type == 1

      # Add transmission event to tree only if it is significant
      if significant_transmissions[event_individual]

        # For significant external transmission events...
        if network.external[event_individual]
          @debug "Event $i is a significant external transmission"

          # Add a new tree and a node
          tree = length(trees) + 1
          height = event_times[event_order[i]]
          node = 1

          push!(trees, Tree(height))
          addnode!(trees[tree])

          tree_id[event_individual] = tree
          transmission_node_id[event_individual] = node

        # For other significant transmission events...
        else

          # Determine source of transmission
          source = findfirst(network.internal[:, event_individual])

          # If infected initial condition, no source
          if source == nothing
            @debug "Event $i is an initial condition"

            # Add a new tree and a node
            tree = length(trees) + 1
            height = events.start_time
            node = 1

            push!(trees, Tree(height))
            addnode!(trees[tree])

            tree_id[event_individual] = tree
            transmission_node_id[event_individual] = node

          # Previous significant transmissions from this source of transmission
          else
            @debug "Event $i is a significant internal transmission"

            source_prior_transmissions = findall(network.internal[source, :][:] .& (event_times[:, 1] .< event_times[event_order[i]]) .& significant_transmissions)

            # For when the source of transmission has yet to be detected...
            if isnan(event_times[source, 2]) || event_times[source, 2] > event_times[event_individual, 1]

              # When there have been prior transmissions from this undetected source
              if length(source_prior_transmissions) > 0
                tree = tree_id[source]
                parentnode = transmission_node_id[source_prior_transmissions[argmax(event_times[source_prior_transmissions, 1])]]
                branch_length = event_times[event_order[i]] - maximum(event_times[source_prior_transmissions, 1])

              # When there have not been any prior transmissions from this undetected source
              else
                tree = tree_id[source]
                parentnode = transmission_node_id[source]
                branch_length = event_times[event_order[i]] - event_times[source, 1]
              end

            # For when the source of transmission has been detected...
            else
              # And detection of source of transmission is most recent, relevant event...
              if length(source_prior_transmissions) == 0 || all(event_times[source, 2] .> event_times[source_prior_transmissions, 1])
                tree = tree_id[source]
                parentnode = observation_internal_node_id[source]
                branch_length = event_times[event_order[i]] - event_times[source, 2]

              # Or some prior transmission is most recent, relevant event...
              else
                tree = tree_id[source]
                parentnode = transmission_node_id[source_prior_transmissions[argmax(event_times[source_prior_transmissions, 1])]]
                branch_length = event_times[event_order[i]] - maximum(event_times[source_prior_transmissions, 1])
              end
            end
            # Add branch and node to tree
            branch!(trees[tree], parentnode, branch_length)

            # Record tree and node ID
            tree_id[event_individual] = tree
            transmission_node_id[event_individual] = length(trees[tree].nodes)
          end
        end
      else
        @debug "Event $i is an insignificant transmission"
      end

    # Infection observation event
    elseif event_type == 2

      # Determine if individual has had significant transmissions prior to being observed as infected
      prior_transmissions = findall(network.internal[event_individual, :][:] .&
                            (event_times[:, 1] .< event_times[event_order[i]]) .&
                            significant_transmissions)
      if length(prior_transmissions) > 0
        tree = tree_id[event_individual]

        # Parent node is the internal node associated with the final pre-observation transmission
        parentnode = transmission_node_id[prior_transmissions[argmax(event_times[prior_transmissions, 1])]]
        branch_length = event_times[event_order[i]] - maximum(event_times[prior_transmissions, 1])

      # Individual has not transmitted to others prior to being observed as infected
      else
        tree = tree_id[event_individual]
        parentnode = transmission_node_id[event_individual]
        branch_length = event_times[event_order[i]] - event_times[event_individual, 1]
      end

      # Individual goes on to have significant transmissions
      if any(network.internal[event_individual, :][:] .& (event_times[:, 1] .>= event_times[event_order[i]]) .& significant_transmissions)
        branch!(trees[tree], parentnode, branch_length)
        observation_internal_node_id[event_individual] = length(trees[tree].nodes)

        # Add zero-length branch so observation event is a leaf node
        branch!(trees[tree], observation_internal_node_id[event_individual], 0.)

      # Individual does not go one to have any significant transmissions
      else
        branch!(trees[tree], parentnode, branch_length)
      end
      # Record node ID
      observation_leaf_node_id[event_individual] = length(trees[tree].nodes)
    end
  end
  return trees, tree_id, observation_leaf_node_id
end
