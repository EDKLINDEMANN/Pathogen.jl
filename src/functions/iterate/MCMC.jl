function iterate!(mcmc::MCMC{T},
                  n::Int64,
                  Σ::Array{Float64, 2},
                  σ::Float64;
                  condition_on_network::Bool=false,
                  event_batches::Int64=1,
                  adapt_cov::Int64=100) where T <: EpidemicModel
  if adapt_cov < 0
    @warn "`adapt_cov` argument indicates the increment in iterations in which the covariance matrix is updated and must be ≧ 0. Setting to 0 for non-adaptive Metropolis-Hastings."
    adapt_cov = 0
  end
  pmeter = Progress(n*length(mcmc.markov_chains), "MCMC progress")
  pchannel = RemoteChannel(()->Channel{Bool}(10), 1)
  mc_Futures = Future[]
  for mc in mcmc.markov_chains
    @debug "Starting MCMC..."
    push!(mc_Futures, @spawn iterate!(mc, mcmc, n, Σ, σ, pchannel, condition_on_network = condition_on_network, event_batches = event_batches, adapt_cov = adapt_cov))
  end
  @debug "MCMC in progress..."
  while !all(isready.(mc_Futures))
    take!(pchannel) && next!(pmeter)
  end
  close(pchannel)
  finish!(pmeter)
  mcmc.markov_chains = [fetch(i) for i in mc_Futures]
  @debug "MCMC complete..."
  return mcmc
end

function iterate!(mcmc::MCMC{T},
                  n::Int64,
                  σ::Float64;
                  condition_on_network::Bool=false,
                  event_batches::Int64=1,
                  adapt_cov::Int64=100) where T <: EpidemicModel
  return iterate!(mcmc,
                  n,
                  Matrix(LinearAlgebra.Diagonal([var(mcmc.risk_priors[i]) for i in 1:length(mcmc.risk_priors)]))/(10*length(mcmc.risk_priors)),
                  σ,
                  condition_on_network = condition_on_network,
                  event_batches = event_batches,
                  adapt_cov = adapt_cov)
end
