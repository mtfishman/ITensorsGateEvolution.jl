using ITensorsGateEvolution
using ITensors
using Test

@testset "ITensorsGateEvolution.jl" begin
  N = 3

  pos = ProductOps()
  pos *= ("Z", 3)
  pos *= ("Y", 2)
  pos *= ("X", 1)

  s = siteinds("qubit", N)
  gates = ops(s, pos)
  ψ0 = productMPS(s, ["0" for n in 1:N])

  # Apply the gates
  ψ = apply(gates, ψ0)

  # Move site 1 to position 3
  ψ′ = movesite(ψ, 1, 3)
  @test siteind(ψ′, 1) == s[2]
  @test siteind(ψ′, 2) == s[3]
  @test siteind(ψ′, 3) == s[1]
  @test prod(ψ) ≈ prod(ψ′)

  # Move the site back
  ψ′′ = movesite(ψ′, 3, 1)
  @test siteind(ψ′′, 1) == s[1]
  @test siteind(ψ′′, 2) == s[2]
  @test siteind(ψ′′, 3) == s[3]
  @test prod(ψ) ≈ prod(ψ′′)
end
