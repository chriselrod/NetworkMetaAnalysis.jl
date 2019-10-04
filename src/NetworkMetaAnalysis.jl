module NetworkMetaAnalysis

using SIMDPirates, ReverseDiffExpressions

export RaggedNetworkMetaAnalysis, OffsetNetworkMetaAnalysis, Zero

abstract type AbstractNetworkMetaAnalysis end

"""
Represents a network as a ragged matrix.
"""
struct RaggedNetworkMetaAnalysis <: AbstractNetworkMetaAnalysis

end
"""
Represents a network as a series of offsets.
"""
struct OffsetNetworkMetaAnalysis <: AbstractNetworkMetaAnalysis

end

"""

"""
function network_meta_analysis_quote()


end

end # module
