mutable struct Simulation{M <: ILM}
  start_time::Float64
  current_time::Float64
  iterations::Int64
  population::Population
  risk_functions::RiskFunctions{M}
  risk_parameters::RiskParameters{M}
  substitution_model::NucleicAcidSubstitutionModel
  disease_states::Vector{DiseaseState}
  transmission_rates::TransmissionRates
  event_rates::EventRates{M}
  events::Events{M}
  transmission_network::TransmissionNetwork

  function Simulation(pop::Population,
                      rf::RiskFunctions{M},
                      rp::RiskParameters{M}) where M <: ILM
    states = fill(State_S, pop.individuals)
    tr = initialize(TransmissionRates, states, pop, rf, rp)
    rates = initialize(EventRates, tr, states, pop, rf, rp)
    events = Events{M}(pop.individuals)
    net = TransmissionNetwork(pop.individuals)
    return new{M}(0.0, 0.0, 0, pop, rf, rp, states, tr, rates, events, net)
  end

  function Simulation(pop::Population,
                      states::Vector{DiseaseState},
                      time::Float64,
                      rf::RiskFunctions{M},
                      rp::RiskParameters{M}) where M <: ILM
    @debug "Initializing $M Simulation with the following starting states:" states
    if length(states) != pop.individuals
      @error "Length of initial disease state vector must match number of individuals"
    elseif !all(in.(states, Ref(_state_progressions[T])))
      @error "All states in initial disease state vector must be valid within specified epidemic model"
    end
    tr = initialize(TransmissionRates, states, pop, rf, rp)
    rates = initialize(EventRates, tr, states, pop, rf, rp)
    events = Events{M}(states, start_time = time)
    net = TransmissionNetwork(pop.individuals)
    return new{M}(time, time, 0, pop, rf, rp, copy(states), tr, rates, events, net)
  end

  function Simulation(pop::Population,
                      states::Vector{DiseaseState},
                      rf::RiskFunctions{M},
                      rp::RiskParameters{M}) where M <: ILM
    @debug "Initializing $M Simulation with the following starting states:" states
    if length(states) != pop.individuals
      @error "Length of initial disease state vector must match number of individuals"
    elseif !all(in.(states, Ref(_state_progressions[T])))
      @error "All states in initial disease state vector must be valid within specified epidemic model"
    end
    tr = initialize(TransmissionRates, states, pop, rf, rp)
    rates = initialize(EventRates, tr, states, pop, rf, rp)
    events = Events{M}(states)
    net = TransmissionNetwork(pop.individuals)
    return new{M}(0.0, 0.0, 0, pop, rf, rp, copy(states), tr, rates, events, net)
  end
end

function Base.show(io::IO, x::Simulation{M}) where M <: ILM
  return print(io, "$M epidemic simulation")
end
