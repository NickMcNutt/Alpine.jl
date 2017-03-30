using SQLite

dir_elements = joinpath(Pkg.dir("Alpine"), "src/elements")

const db = SQLite.DB(joinpath(dir_elements, "elements.db"))

function element_properties()
    df = SQLite.columns(db, "elements")
    keys = SQLite.columns(Alpine.db, "elements")[:name]
    vals = SQLite.columns(Alpine.db, "elements")[:type]
    Dict(zip(keys, vals))
end

function element_property(symbol::Symbol, property::Symbol)
    df = SQLite.query(db, "SELECT * FROM elements WHERE symbol == ?;", values = [string(symbol)])
    df[property][1]
end
