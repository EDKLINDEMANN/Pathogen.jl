function next!(mcmc::MCMC{M},
               Σ::Array{Float64, 2},
               σ::Float64) where M <: TNILM
  @simd for mc in mcmc.markov_chains
    next!(mc, mcmc, Σ, σ)
  end
  return mcmc
end

function next!(mcmc::MCMC{M},
               Σrp::Array{Float64, 2},
               Σsm::Array{Float64, 2},
               σ::Float64) where M <: PhyloILM
  @simd for mc in mcmc.markov_chains
    next!(mc, mcmc, Σrp, Σsm, σ)
  end
  return mcmc
end
