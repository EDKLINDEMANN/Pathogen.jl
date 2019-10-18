mutable struct MarkovChain{M <: ILM}
  iterations::Int64
  events::Vector{Events{M}}
  transmission_network::Vector{TransmissionNetwork}
  risk_parameters::Vector{RiskParameters{M}}
  substitutionmodel::Union{Nothing, Vector{NucleicAcidSubstitutionModel}}
  log_posterior::Vector{Float64}
  cov::OnlineStats.CovMatrix

  function MarkovChain(e::Events{M},
                       n::TransmissionNetwork,
                       rp::RiskParameters{M},
                       lp::Float64) where {
                       M <: PhyloILM}
    return new{M}(0, [e], [n], [rp], nothing, [lp], OnlineStats.CovMatrix())
  end

  function MarkovChain(e::Events{M},
                       n::TransmissionNetwork,
                       rp::RiskParameters{M},
                       sm::NucleicAcidSubstitutionModel,
                       lp::Float64) where {
                       M <: PhyloILM}
    return new{M}(0, [e], [n], [rp], [sm], [lp], OnlineStats.CovMatrix())
  end
end

function Base.length(x::MarkovChain{M}) where M <: ILM
  return x.iterations
end

function Base.show(io::IO, x::MarkovChain{M}) where M <: ILM
  return print(io, "$M Markov chain (iterations = $(length(x)))")
end
