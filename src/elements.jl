using Colors

dir = "elements/"
filename_colors = joinpath(dir, "element_colors.txt")
filename_radii = joinpath(dir, "element_radii.txt")
filename_numbers = joinpath(dir, "element_numbers.txt")

function read_element_properties{T}(f::Function, result_type::Type{T}, filename::AbstractString)
    data, header = readdlm(filename, '\t', AbstractString, header = true, comments = false)

    props = [prop => data[:, i] for (i, prop) in enumerate(header)]
    elements = props["element"]
    delete!(props, "element")

    result = Dict{Symbol, Dict{Symbol, T}}()
    for (i, el) in enumerate(elements)
        element = Symbol(el)
        result[element] = Dict{Symbol, T}()
        for (key, val) in props
            length(val[i]) > 0 || continue
            result[element][Symbol(key)] = f(val[i])::T
        end
    end

    return result
end

const element_colors = read_element_properties(NTuple{3, Float64}, filename_colors) do color
    h = "#" * color
    c = parse(Colorant, h)
    (Float64(c.r), Float64(c.g), Float64(c.b))
end

const element_radii = read_element_properties(Float64, filename_radii) do radius
    parse(Float64, radius)
end

const element_numbers = atomic_numbers = read_element_properties(Int, filename_numbers) do number
    parse(Int, number)
end
