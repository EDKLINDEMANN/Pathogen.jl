# TODO: left off implementing adaptive COV for phyloilms


function iterate!(mc::MarkovChain{M},
                  mcmc::MCMC{M},
                  n::Int64,
                  Σ::Array{Float64, 2},
                  σ::Float64;
                  condition_on_network::Bool=false,
                  event_batches::Int64=1,
                  adapt_cov::Int64=100) where {
                  S <: DiseaseStateSequence,
                  M <: TNILM{S}}
  if adapt_cov < 0
    @warn "`adapt_cov` argument indicates the increment in iterations in which the covariance matrix is updated and must be ≧ 0. Setting to 0 for non-adaptive Metropolis-Hastings."
    adapt_cov = 0
  end
  pmeter = Progress(n, "MCMC progress")
  use_adapted_cov = false
  local adapted_cov
  if (adapt_cov != 0 & mc.iterations > adapt_cov) && LinearAlgebra.isposdef(value(mc.Σrp))
    use_adapted_cov = true
    adapted_cov = OnlineStats.value(mc.Σrp) * 2.38^2 / length(mc.risk_parameters[1])
  end
  for i = 1:n
    if use_adapted_cov
      next!(mc, mcmc, adapted_cov, σ, condition_on_network = condition_on_network, event_batches = event_batches)
    else
      next!(mc, mcmc, Σ, σ, condition_on_network = condition_on_network, event_batches = event_batches)
    end
    next!(pmeter)
    @logmsg LogLevel(-5000) "MCMC progress" progress = i/n
    if adapt_cov != 0 && mod(i, adapt_cov) == 0
      OnlineStats.fit!(mc.Σrp, convert(Array{Float64, 2}, mc.risk_parameters[end-adapt_cov:end]))
      if LinearAlgebra.isposdef(value(mc.Σrp))
        if !use_adapted_cov
          use_adapted_cov = true
          @debug "Now using adapted covariance matrix for adaptive MCMC on core $(Distributed.myid()) (updated every $adapt_cov iterations)"
        end
        adapted_cov = OnlineStats.value(mc.Σrp) * 2.38^2 / length(mc.risk_parameters[1])
        @debug "Covariance matrix updated for Adaptive Metropolis-Hastings MCMC" Covariance = adapted_cov
      end
    end
  end
  return mc
end

function iterate!(mc::MarkovChain{M},
                  mcmc::MCMC{M},
                  n::Int64,
                  Σ::Array{Float64, 2},
                  σ::Float64,
                  progress_channel::RemoteChannel;
                  condition_on_network::Bool=false,
                  event_batches::Int64=1,
                  adapt_cov::Int64=100) where {
                  S <: DiseaseStateSequence,
                  M <: TNILM{S}}
  use_adapted_cov = false
  local adapted_cov
  if (adapt_cov != 0 & mc.Σrp.n >= adapt_cov) && LinearAlgebra.isposdef(value(mc.Σrp))
    use_adapted_cov = true
    adapted_cov = OnlineStats.value(mc.Σrp) * 2.38^2 / length(mc.risk_parameters[1])
  end
  for i = 1:n
    if use_adapted_cov
      next!(mc, mcmc, adapted_cov, σ, condition_on_network = condition_on_network, event_batches = event_batches)
    else
      next!(mc, mcmc, Σ, σ, condition_on_network = condition_on_network, event_batches = event_batches)
    end
    put!(progress_channel, true)
    if adapt_cov != 0 && mod(i, adapt_cov) == 0
      OnlineStats.fit!(mc.Σrp, eachrow(convert(Array{Float64, 2}, mc.risk_parameters[end-adapt_cov:end])))
      if LinearAlgebra.isposdef(value(mc.Σrp))
        if !use_adapted_cov
          use_adapted_cov = true
          @debug "Now using adapted covariance matrix for adaptive MCMC on core $(Distributed.myid()) (updated every $adapt_cov iterations)"
        end
        adapted_cov = OnlineStats.value(mc.Σrp) * 2.38^2 / length(mc.risk_parameters[1])
        @debug "Covariance matrix updated for Adaptive Metropolis-Hastings MCMC" Covariance = adapted_cov
      end
    end
  end
  return mc
end

function iterate!(mc::MarkovChain{M},
                  mcmc::MCMC{M},
                  n::Int64,
                  σ::Float64;
                  condition_on_network::Bool=false,
                  event_batches::Int64=1,
                  adapt_cov::Int64=100) where {
                  S <: DiseaseStateSequence,
                  M <: TNILM{S}}
  return iterate!(mc,
                  mcmc,
                  n,
                  Matrix(LinearAlgebra.Diagonal([var(mcmc.risk_priors[i]) for i in 1:length(mcmc.risk_priors)])),
                  σ,
                  condition_on_network = condition_on_network,
                  event_batches = event_batches,
                  adapt_cov = adapt_cov)
end

function iterate!(mc::MarkovChain{M},
                  mcmc::MCMC{M},
                  n::Int64,
                  Σ::Array{Float64, 2},
                  σ::Float64;
                  condition_on_network::Bool=false,
                  event_batches::Int64=1,
                  adapt_cov::Int64=100) where {
                  S <: DiseaseStateSequence,
                  M <: TNILM{S}}
  if adapt_cov < 0
    @warn "`adapt_cov` argument indicates the increment in iterations in which the covariance matrix is updated and must be ≧ 0. Setting to 0 for non-adaptive Metropolis-Hastings."
    adapt_cov = 0
  end
  pmeter = Progress(n, "MCMC progress")
  use_adapted_cov = false
  local adapted_cov
  if (adapt_cov != 0 & mc.iterations > adapt_cov) && LinearAlgebra.isposdef(value(mc.Σrp))
    use_adapted_cov = true
    adapted_cov = OnlineStats.value(mc.Σrp) * 2.38^2 / length(mc.risk_parameters[1])
  end
  for i = 1:n
    if use_adapted_cov
      next!(mc, mcmc, adapted_cov, σ, condition_on_network = condition_on_network, event_batches = event_batches)
    else
      next!(mc, mcmc, Σ, σ, condition_on_network = condition_on_network, event_batches = event_batches)
    end
    next!(pmeter)
    @logmsg LogLevel(-5000) "MCMC progress" progress = i/n
    if adapt_cov != 0 && mod(i, adapt_cov) == 0
      OnlineStats.fit!(mc.Σrp, convert(Array{Float64, 2}, mc.risk_parameters[end-adapt_cov:end]))
      if LinearAlgebra.isposdef(value(mc.Σrp))
        if !use_adapted_cov
          use_adapted_cov = true
          @debug "Now using adapted covariance matrix for adaptive MCMC on core $(Distributed.myid()) (updated every $adapt_cov iterations)"
        end
        adapted_cov = OnlineStats.value(mc.Σrp) * 2.38^2 / length(mc.risk_parameters[1])
        @debug "Covariance matrix updated for Adaptive Metropolis-Hastings MCMC" Covariance = adapted_cov
      end
    end
  end
  return mc
end

function iterate!(mc::MarkovChain{M},
                  mcmc::MCMC{M},
                  n::Int64,
                  Σrp::Array{Float64, 2},
                  Σsm::Array{Float64, 2},
                  σ::Float64,
                  progress_channel::RemoteChannel;
                  condition_on_network::Bool=false,
                  event_batches::Int64=1,
                  adapt_cov::Int64=100) where {
                  S <: DiseaseStateSequence,
                  M <: PhyloILM{S}}
  use_adapted_cov = false
  local adapted_cov
  if (adapt_cov != 0 & mc.Σrp.n >= adapt_cov) && LinearAlgebra.isposdef(value(mc.Σrp)) &&
    use_adapted_cov = true
    adapted_cov_rp = OnlineStats.value(mc.Σrp) * 2.38^2 / length(mc.risk_parameters[1])
    adapted_cov_sm = OnlineStats.value(mc.Σsm) * 2.38^2 / length(propertynames(mc.substitution_model[1]))
  end
  for i = 1:n
    if use_adapted_cov
      next!(mc, mcmc, adapted_cov, σ, condition_on_network = condition_on_network, event_batches = event_batches)
    else
      next!(mc, mcmc, Σrp, Σsm, σ, condition_on_network = condition_on_network, event_batches = event_batches)
    end
    put!(progress_channel, true)
    if adapt_cov != 0 && mod(i, adapt_cov) == 0
      OnlineStats.fit!(mc.Σrp, eachrow(convert(Array{Float64, 2}, mc.risk_parameters[end-adapt_cov:end])))
      if LinearAlgebra.isposdef(value(mc.Σrp))
        if !use_adapted_cov
          use_adapted_cov = true
          @debug "Now using adapted covariance matrix for adaptive MCMC on core $(Distributed.myid()) (updated every $adapt_cov iterations)"
        end
        adapted_cov = OnlineStats.value(mc.Σrp) * 2.38^2 / length(mc.risk_parameters[1])
        @debug "Covariance matrix updated for Adaptive Metropolis-Hastings MCMC" Covariance = adapted_cov
      end
    end
  end
  return mc
end

function iterate!(mc::MarkovChain{M},
                  mcmc::MCMC{M},
                  n::Int64,
                  σ::Float64;
                  condition_on_network::Bool=false,
                  event_batches::Int64=1,
                  adapt_cov::Int64=100) where {
                  S <: DiseaseStateSequence,
                  M <: TNILM{S}}
  return iterate!(mc,
                  mcmc,
                  n,
                  Matrix(LinearAlgebra.Diagonal([var(mcmc.risk_priors[i]) for i in 1:length(mcmc.risk_priors)])),
                  σ,
                  condition_on_network = condition_on_network,
                  event_batches = event_batches,
                  adapt_cov = adapt_cov)
end
