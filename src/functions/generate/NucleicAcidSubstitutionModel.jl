function generate(::Type{SM},
                  priors::Vector{UnivariateDistribution) where {
                  SM <: NucleicAcidSubstitutionModel}
  return SM([rand(x) for x in priors])
end
