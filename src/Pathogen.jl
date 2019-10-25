module Pathogen

  # Dependencies
  using Distributed,
      DataFrames,
      Distributions,
      Logging,
      StatsBase,
      Statistics,
      ProgressMeter,
      LinearAlgebra,
      OnlineStats,
      PhyloModels

  import ProgressMeter.next!

  # Types
  include("types/ILM.jl")
  include("types/DiseaseState.jl")
  include("types/Population.jl")
  include("types/Events/Event.jl")
  include("types/Events/Events.jl")
  include("types/Events/EventExtents.jl")
  include("types/Events/EventObservations.jl")
  include("types/Events/EventRates.jl")
  include("types/Risks/RiskFunctions.jl")
  include("types/Risks/RiskParameters.jl")
  include("types/Risks/RiskPriors.jl")
  include("types/Transmissions/Transmission.jl")
  include("types/Transmissions/TransmissionNetwork.jl")
  include("types/Transmissions/TransmissionRates.jl")
  include("types/Simulation.jl")
  include("types/MarkovChain.jl")
  include("types/MCMC.jl")

  # Functions
  include("functions/_pathway_from.jl")
  include("functions/_pathway_to.jl")
  include("functions/_count_by_state.jl")
  include("functions/_accept.jl")
  include("functions/_bounds.jl")

  # Generate functions
  include("functions/generate/Event.jl")
  include("functions/generate/Events.jl")
  include("functions/generate/Transmission.jl")
  include("functions/generate/PhyloTree.jl")
  include("functions/generate/RiskParameters.jl")
  include("functions/generate/NucleicAcidSubstitutionModel.jl")

  # Logprior functions
  include("functions/logprior/RiskParameters.jl")
  include("functions/logprior/NucleicAcidSubstitutionModel.jl")

  # Initialize functions
  include("functions/initialize/TransmissionRates.jl")
  include("functions/initialize/EventRates.jl")
  include("functions/initialize/MarkovChain.jl")

  # Next! functions
  include("functions/next!/Simulation.jl")
  include("functions/next!/MarkovChain.jl")
  include("functions/next!/MCMC.jl")

  # Iterate functions
  include("functions/iterate/MarkovChain.jl")
  include("functions/iterate/MCMC.jl")

  include("functions/observe.jl")
  include("functions/update!.jl")
  include("functions/next!.jl")
  include("functions/simulate!.jl")
  include("functions/loglikelihood.jl")
  include("functions/start!.jl")
  include("functions/iterate!.jl")

  for name in names(PhyloModels)
    @eval export $(name)
  end

  export
    SEIR, SEI, SIR, SI,
    DiseaseState,
    State_S, State_E, State_I, State_R,
    RiskFunctions, RiskParameters, RiskPriors,
    Population,
    TransmissionNetwork,
    Simulation,
    next!, simulate!,
    Events, EventObservations, EventExtents,
    observe,
    MCMC, start!, iterate!,
    generate_tree
end
