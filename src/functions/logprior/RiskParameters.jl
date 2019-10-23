function logprior(rparams::RiskParameters{M},
                  rpriors::RiskPriors{M}) where {
                  S <: DiseaseStateSequence,
                  M <: ILM{S}}
  lprior = 0.0
  for i in 1:length(rpriors.sparks)
    lprior == -Inf && break
    lprior += loglikelihood(rpriors.sparks[i], rparams.sparks[i])
  end
  for i in 1:length(rpriors.susceptibility)
    lprior == -Inf && break
    lprior += loglikelihood(rpriors.susceptibility[i], rparams.susceptibility[i])
  end
  for i in 1:length(rpriors.infectivity)
    lprior == -Inf && break
    lprior += loglikelihood(rpriors.infectivity[i], rparams.infectivity[i])
  end
  for i in 1:length(rpriors.transmissibility)
    lprior == -Inf && break
    lprior += loglikelihood(rpriors.transmissibility[i], rparams.transmissibility[i])
  end
  if S in [SEIR; SEI]
    for i in 1:length(rpriors.latency)
      lprior == -Inf && break
      lprior += loglikelihood(rpriors.latency[i], rparams.latency[i])
    end
  end
  if S in [SEIR; SIR]
    for i in 1:length(rpriors.removal)
      lprior == -Inf && break
      lprior += loglikelihood(rpriors.removal[i], rparams.removal[i])
    end
  end
  return lprior
end
