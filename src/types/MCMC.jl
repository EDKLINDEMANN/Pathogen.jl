mutable struct MCMC{M <: ILM}
  event_observations::EventObservations{M}
  event_extents::EventExtents{M}
  population::Population
  starting_states::Vector{DiseaseState}
  risk_functions::RiskFunctions{M}
  risk_priors::RiskPriors{M}
  substitution_model::Union{Nothing, Type{<:NucleicAcidSubstitutionModel}}
  substitution_model_priors::Union{Nothing, Vector{UnivariateDistribution}}
  markov_chains::Vector{MarkovChain{M}}

  function MCMC(obs::EventObservations{M},
                ee::EventExtents{M},
                pop::Population,
                states::Vector{DiseaseState},
                rf::RiskFunctions{M},
                rp::RiskPriors{M}) where M <: TNILM
    return new{M}(obs, ee, pop, states, rf, rp, nothing, nothing, MarkovChain{M}[])
  end

  function MCMC(obs::EventObservations{M},
                ee::EventExtents{M},
                pop::Population,
                states::Vector{DiseaseState},
                rf::RiskFunctions{M},
                rp::RiskPriors{M},
                sm::Type{N},
                smp::Vector{UnivariateDistribution}) where {
                M <: PhyloILM,
                N <: NucleicAcidSubstitutionModel}
    return new{M}(obs, ee, pop, states, rf, rp, sm, smp, MarkovChain{M}[])
  end
end

function MCMC(obs::EventObservations{M},
              ee::EventExtents{M},
              pop::Population,
              rf::RiskFunctions{M},
              rp::RiskPriors{M}) where M <: TNILM
  return MCMC(obs, ee, pop, fill(State_S, pop.individuals), rf, rp)
end

function MCMC(obs::EventObservations{M},
              ee::EventExtents{M},
              pop::Population,
              rf::RiskFunctions{M},
              rp::RiskPriors{M},
              sm::Type{N},
              smp::Vector{UnivariateDistribution}) where {
              M <: PhyloILM,
              N <: NucleicAcidSubstitutionModel}
  return MCMC(obs, ee, pop, fill(State_S, pop.individuals), rf, rp, sm, smp)
end

function Base.show(io::IO, x::MCMC{M}) where M <: ILM
  return print(io, "$M MCMC with $(length(x.markov_chains)) chains")
end
