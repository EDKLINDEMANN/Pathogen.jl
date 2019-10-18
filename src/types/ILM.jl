abstract type DiseaseStateSequence end
abstract type SEIR <: DiseaseStateSequence end
abstract type SIR <: DiseaseStateSequence end
abstract type SEI <: DiseaseStateSequence end
abstract type SI <: DiseaseStateSequence end

abstract type TransmissionNetworkILM{T <: DiseaseStateSequence} end
const TNILM{T} = TransmissionNetworkILM{T}

abstract type PhylodynamicILM{T <: DiseaseStateSequence} end
const PhyloILM{T} = PhylodynamicILM{T}

const IndividualLevelModel{T} = Union{TNILM{T}, PhyloILM{T}}
const ILM{T} = IndividualLevelModel{T}
