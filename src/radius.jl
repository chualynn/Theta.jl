"""
    radius_ellipsoid(ϵ, T, Y, ρ, nderivs, tol)

Compute radius of the ellipsoid used for computing the Riemann theta function or its derivatives.

# Arguments
- `ϵ::Real`: absolute error in value of theta function or its derivatives
- `T::Array{<:Real}`: upper triangular matrix in Cholesky factorization of Y
- `Y::Array{<:Real}`: imaginary part of Riemann matrix
- `ρ::Real`: norm of shortest vector in lattice generated by T
- `nderivs::Integer`: order of derivative to be computed
- `tol::Real=1e-15`: error tolerance in optimization routine
"""
function radius_ellipsoid(ϵ::Real, T::Array{<:Real}, Y::Array{<:Real}, ρ::Real, nderivs::Integer, tol::Real=1e-15)
    m = size(T,2); # number of columns in matrix 
    norm_Tinv = norm(m > 1 ? inv(T) : inv.(T)); # norm of T⁻¹
    r = 0.5*(sqrt(m + 2*nderivs + sqrt(m^2 + 8*nderivs))+ρ); 
    opt = Opt(:LN_COBYLA, 1); # use the COBYLA algorithm (Constrained Optimization BY Linear Approximations) in NLopt
    min_objective!(opt, (x, grad) -> radius_ellipsoid_helper(x[1], ϵ, ρ, m, norm_Tinv, nderivs)); # set objective function
    lower_bounds!(opt, [0.]); # add lower bound so that optimization parameter is positive
    inequality_constraint!(opt, (x, grad) -> - radius_ellipsoid_helper(x[1], ϵ, ρ, m, norm_Tinv, nderivs)); # set objective to be positive
    ftol_abs!(opt, tol); # set tolerance on value of objective function
    optx = optimize(opt, [r])[2];
    if optx[1] > r
        r = optx[1];
    end
    return r;
end

"""
    radius_ellipsoid_helper(x, grad, ϵ, ρ, g, T, L, nderivs)

Helper function for computing radius of ellipsoid for computing the theta function or its derivatives, in [`radius_ellipsoid`](@ref).
"""
function radius_ellipsoid_helper(r::Real, ϵ::Real, ρ::Real, g::Integer, norm_Tinv::Real, nderivs::Integer)
    if nderivs == 0
        return abs(g*2^(g-1)*sf_gamma_inc(g/2, (r-ρ/2)^2)/ρ^g) - ϵ;
    else
        return 2^(nderivs-1) * g * (2/ρ)^g * sum([binomial(nderivs, k) * π^(nderivs-k/2) * norm_Tinv^k * g^((nderivs-k)/2) * sf_gamma_inc((g+k)/2, (r-ρ/2)^2) for k=0:nderivs]) - ϵ;
    end
end


