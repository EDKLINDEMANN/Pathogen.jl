struct EventObservations{M <: ILM}
  infection::Union{Nothing, Vector{Float64}}
  removal::Union{Nothing, Vector{Float64}}
  seq::Union{Nothing, Vector{Union{Nothing, GeneticSeq}}}
  individuals::Int64

  function EventObservations{M}(in::VF,
                                re::VF,
                                seq::VG) where {
                                VF <: Vector{Float64},
                                G  <: GeneticSeq,
                                VG <: Vector{G},
                                S  <: Union{SEIR, SIR},
                                M  <: PhyloILM{S}}
  if length(unique(length.([in, re, seq])) != 1)
    throw(DimensionMismatch("Observation vector lengths do not match"))
  end
    return new{M}(in, re, seq, length(in))
  end

  function EventObservations{M}(in::VF,
                                re::VF) where {
                                VF <: Vector{Float64},
                                S  <: Union{SEIR, SIR},
                                M  <: TNILM{S}}
  if length(unique(length.([in, re])) != 1)
    throw(DimensionMismatch("Observation vector lengths do not match"))
  end
    return new{M}(in, re, nothing, length(in))
  end

  function EventObservations{M}(in::Vector{Float64},
                                seq::VG) where {
                                G  <: GeneticSeq,
                                VG <: Vector{G},
                                S  <: Union{SEI, SI},
                                M  <: PhyloILM{S}}
  if length(unique(length.([in, seq])) != 1)
    throw(DimensionMismatch("Observation vector lengths do not match"))
  end
    return new{M}(in, nothing, seq, length(in))
  end

  function EventObservations{M}(in::Vector{Float64}) where {
                                S  <: Union{SEI, SI},
                                M  <: TNILM{S}}
    return new{M}(in, nothing, nothing, length(in))
  end
end


function Base.show(io::IO,
                   x::EventObservations{M}) where {
                   S <: DiseaseStateSequence,
                   M <: ILM{S}}
  print(io, "$M event observations (n=$(x.individuals))")
end
