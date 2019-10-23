function logprior(sm::NucleicAcidSubstitutionModel,
                  sm_priors::Vector{UnivariateDistributions})
  lprior = 0.0
  for param in [getproperty(x1, θ) for θ in [propertynames(x1)...]], dist in sm_priors
    lprior += loglikelihood(dist, param)
    lprior == -Inf && break
  end
  return lprior
end
