"""
    min_even_theta_char(R)

Compute the even theta characteristic for which the theta constant at the input matrix in R has the smallest absolute value.

# Arguments
- `R::RiemannMatrix`: RiemannMatrix type containing matrix.

# Examples
```julia
julia> R = RiemannMatrix(random_siegel(5));
julia> min_even_theta_char(R)
```
"""
function min_even_theta_char(R::RiemannMatrix)
    z = zeros(R.g);
    even_chars = even_theta_char(R.g);
    min_theta, char_index = findmin([abs(theta(z, R, char=c)) for c in even_chars])
    char = even_chars[char_index];
    return char
end

"""
    hessian(z, R, char=[])

Compute the hessian matrix of the theta function at input z and input matrix in R, with optional characteristics. 

# Arguments
- `z::Array{<:Number}`: Array of size g, each entry is a complex number.
- `R::RiemannMatrix`: RiemannMatrix type containing matrix.
- `char::Array{}=[]`: 2-element array containing 2 arrays of size g, where entries are either 0 or 1.

# Examples
```julia
julia> hessian([0,0], RiemannMatrix([1+im -1; -1 1+im]), char=[[1,0],[1,1]])
```
"""
function hessian(z::Array{<:Number}, R::RiemannMatrix, char::Array{}=[])
    M = zeros(Complex{Float64}, R.g, R.g);
    for i=1:R.g
        for j=i:R.g
            derivs = [zeros(R.g), zeros(R.g)];
            derivs[1][i] = 1;
            derivs[2][j] = 1;
            M[i,j] = theta(z, R, char=char, derivs=derivs);
            M[j,i] = M[i,j];
        end
    end
    return M
end


"""
    hessian(z, τ, char=[], siegel=true, ϵ=1.0e-12)

Compute the hessian matrix of the theta function at inputs z and τ, with optional characteristics, and optional inputs specifying whether to compute the Siegel transformation and the error in the value of the theta function. 

# Arguments
- `z::Array{<:Number}`: Array of size g, each entry is a complex number.
- `τ::Array{<:Number}`: 2-dimensional array of complex numbers
- `char::Array{}=[]`: 2-element array containing 2 arrays of size g, where entries are either 0 or 1.
- `siegel::Bool=false`: do a siegel transformation on τ if true, and use the input matrix otherwise
- `ϵ::Real=1.0e-12`: absolute error in value of theta function or its derivatives


# Examples
```julia
julia> hessian([0,0], [1+im -1; -1 1+im], char=[[1,0],[1,1]])
```
"""
function hessian(z::Array{<:Number}, τ::Array{<:Number}, char::Array{}=[]; siegel::Bool=false, ϵ::Real=1.0e-12)
    R = RiemannMatrix(τ, siegel=siegel, ϵ=ϵ, nderivs=2);
    return hessian(z, R, char)
end


"""
    schottky_null(R, tol=1.0e-8)

Compute the even theta characteristic where the theta constant vanishes and the hessian matrix at the characteristic, as well as its nonzero eigenvalues, up to the input tolerance. Returns an array containing the characteristic, the hessian and the nonzero eigenvalues. If there is no such characteristic, return the absolute value of the theta constant at the smallest even characteristic.

# Arguments
- `R::RiemannMatrix`: RiemannMatrix type containing matrix.
- `tol::Real=1.0e-8`: tolerance for deciding whether eigenvalues are nonzero.

# Examples
```julia
julia> schottky_null(RiemannMatrix([1+im -1; -1 1+im]))
```
"""
function schottky_null(R::RiemannMatrix, tol::Real=1.0e-8)
    char = min_even_theta_char(R);
    z = zeros(R.g);
    theta_value = abs(theta(z, R, char=char));
    if theta_value < tol # check that there is a vanishing theta null
        H = hessian(z, R, char);
        hessian_eigvals = filter(x -> abs(x) > tol, eigvals(H));
        return [char, H, hessian_eigvals]
    end
    return theta_value
end


"""
    schottky_null(τ, tol=1.0e-8)

Compute the even theta characteristic where the theta constant vanishes and the hessian matrix at the characteristic, as well as its nonzero eigenvalues, up to the input tolerance. Returns an array containing the characteristic, the hessian and the nonzero eigenvalues. If there is no such characteristic, return the absolute value of the theta constant at the smallest even characteristic.

# Arguments
- `τ::Array{<:Number}`: 2-dimensional array of complex numbers
- `tol::Real=1.0e-8`: tolerance for deciding whether eigenvalues are nonzero.

# Examples
```julia
julia> schottky_null([1+im -1; -1 1+im])
```
"""
function schottky_null(τ::Array{<:Number}, tol::Real=1.0e-8)
    R = RiemannMatrix(τ, ϵ=tol, nderivs=2);
    return schottky_null(R, tol)
end
