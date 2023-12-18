using Astrometry
using Test

@testset "Astrometry.jl" begin
    include("astrotests.jl")
    include("sofatests.jl")
end
