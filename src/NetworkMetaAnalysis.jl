module NetworkMetaAnalysis

using SIMDPirates, ReverseDiffExpressionsBase, PaddedMatrices

export NetworkMetaAnalysis, RaggedNetwork, OffsetNetwork

abstract type AbstractNetwork end

"""
Represents a network as a ragged matrix.
"""
struct RaggedNetwork <: AbstractNetwork

end
"""
Represents a network as a series of offsets.
"""
struct OffsetNetwork <: AbstractNetwork

end

abstract type AbstractNetworkEffect{P,T,L} <: PaddedMatrices.AbstractMutableFixedSizeVector{P,T,L}
struct FixedEffect <: AbstractNetworkEffect{P,T,L}

end

struct RandomEffect <: AbstractNetworkEffect{P,T,L}

end

"""

"""
function network_meta_analysis_quote()


end

"""
α = FixedEffect(  baselines, δα )
θ = RandomEffect( effects, δθ, stdev )
(α,θ) ~ NetworkMetaAnalysis( network, transforms )
"""
@generated function NetworkMetaAnalysis()

end

end # module
