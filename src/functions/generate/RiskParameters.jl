

function generate(::Type{RiskParameters{M}},
                  rpriors::RiskPriors{M}) where {
                  S <: DiseaseStateSequence,
                  M <: ILM{S}}
  sparks = Float64[rand(x) for x in rpriors.sparks]
  susceptibility = Float64[rand(x) for x in rpriors.susceptibility]
  infectivity = Float64[rand(x) for x in rpriors.infectivity]
  transmissibility = Float64[rand(x) for x in rpriors.transmissibility]
  if T in [SEIR; SEI]
    latency = Float64[rand(x) for x in rpriors.latency]
  end
  if T in [SEIR; SIR]
    removal = Float64[rand(x) for x in rpriors.removal]
  end
  if T == SEIR
    return RiskParameters{M}(sparks,
                             susceptibility,
                             infectivity,
                             transmissibility,
                             latency,
                             removal)
  elseif T == SEI
    return RiskParameters{M}(sparks,
                             susceptibility,
                             infectivity,
                             transmissibility,
                             latency)
  elseif T == SIR
    return RiskParameters{M}(sparks,
                             susceptibility,
                             infectivity,
                             transmissibility,
                             removal)
  elseif T == SI
    return RiskParameters{M}(sparks,
                             susceptibility,
                             infectivity,
                             transmissibility)
  end
end

function generate(::Type{RiskParameters},
                  mcmc::MCMC{M}) where {
                  S <: DiseaseStateSequence,
                  M <: ILM{S}}
  return generate(RiskParameters{M}, mcmc.risk_priors)
end

function generate(::Type{RiskParameters{M}},
                  last_rparams::RiskParameters{M},
                  Σ::Array{Float64, 2}) where {
                  S <: DiseaseStateSequence,
                  M <: ILM{S}}
  rparams_vector = rand(MvNormal(convert(Vector{Float64}, last_rparams), Σ))
  return _like(last_rparams, rparams_vector)
end

function generate(::Type{RiskParameters{M}},
                  mc::MarkovChain{M},
                  Σ::Array{Float64, 2}) where {
                  S <: DiseaseStateSequence,
                  M <: ILM{S}}
  return generate(RiskParameters{M}, mc.risk_parameters[end], Σ)
end
