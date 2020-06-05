include("productops.jl")

N = 3

pos = ProductOps()
pos *= ("Z", 3)
pos *= ("Y", 2)
pos *= ("X", 1)

s = siteinds("qubit", N)
gates = ops(s, pos)
ψ0 = productMPS(s, ["0" for n in 1:N])

# Apply the gates
ψ = product(gates, ψ0)

# Move site 1 to position 3
ψ′ = movesite(ψ, 1, 3)
@show siteind(ψ′, 1) == s[2]
@show siteind(ψ′, 2) == s[3]
@show siteind(ψ′, 3) == s[1]
@show norm(prod(ψ) - prod(ψ′))

# Move the site back
ψ′′ = movesite(ψ′, 3, 1)
@show siteind(ψ′′, 1) == s[1]
@show siteind(ψ′′, 2) == s[2]
@show siteind(ψ′′, 3) == s[3]
@show norm(prod(ψ) - prod(ψ′′))

