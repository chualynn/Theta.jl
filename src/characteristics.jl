"""
    parity(char)

Compute whether a characteristic is even or odd.

# Arguments
- `char::Array{}`: array consisting of two arrays of the same length with 0 or 1 entries.
"""
function parity_char(char::Array{})
    return rem.((transpose(char[1])*char[2]),2);
end


"""
    remainder_char(char)

Compute remainder modulo 2 of all the entries of a characteristic.

# Arguments
- `char::Array{}`: array consisting of two arrays of the same length with 0 or 1 entries.
"""
function remainder_char(char::Array{})
    return [rem.(char[1], 2), rem.(char[2], 2)];
end

"""
    theta_char(g)

Compute all theta characteristics of genus g.

# Arguments
- `g::Integer`: genus

# Examples
```julia
julia> theta_char(4)
```
"""
function theta_char(g::Integer)
    chars = [digits(i, base=2, pad=2*g) for i = 0:2^(2*g)-1];
    return [[c[1:g], c[g+1:2*g]] for c in chars];
end

"""
    even_theta_char(g)

Compute all even theta characteristics of genus g.

# Arguments
- `g::Integer`: genus

# Examples
```julia
julia> even_theta_char(4)
```
"""
function even_theta_char(g::Integer)
    char = theta_char(g);
    even_char = filter(x -> parity_char(x) == 0, char);
    return even_char;
end

"""
    odd_theta_char(g)

Compute all odd theta characteristics of genus g.

# Arguments
- `g::Integer`: genus

# Examples
```julia
julia> odd_theta_char(4)
```
"""
function odd_theta_char(g::Integer)
    char = theta_char(g);
    odd_char = filter(x -> parity_char(x) == 1, char);
    return odd_char;
end

"""
    check_azygetic(chars)

Check if a list of characteristics is azygetic.

# Arguments
- `chars::Array{}`: array where each entry is an array consisting of two arrays of the same length with 0 or 1 entries.

# Examples
```julia
julia> check_azygetic([[[1,0,1,0], [1,0,1,0]], [[0,0,0,1], [1,0,0,0]], [[0,0,1,1], [1,0,1,1]]])
"""
function check_azygetic(chars::Array{})
    for i = 1:length(chars), j = i+1:length(chars), k = j+1:length(chars)
        if parity_char(remainder_char(chars[i] + chars[j] + chars[k])) == 0
            return false
        end
    end
    return true
end


