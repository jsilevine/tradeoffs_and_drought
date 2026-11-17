##-------------------------------------------------------------------
## Utility functions 
##-------------------------------------------------------------------

module Utils

export projectdir, datadir, figdir, srcdir, tabledir

"""
    projectdir(parts...)

Return an absolute path inside the project root.  
`projectdir()` alone gives the root;  
`projectdir("data", "file.csv")` gives `<root>/data/file.csv`.
"""
function projectdir(parts...)
    # climb up from current file until Project.toml is found
    dir = @__DIR__
    while !isfile(joinpath(dir, "Project.toml"))
        newdir = dirname(dir)
        newdir == dir && error("No Project.toml found above $(@__DIR__)")
        dir = newdir
    end
    return joinpath(dir, parts...)
end

# convenience wrappers
datadir(parts...)  = projectdir("data", parts...)
figdir(parts...)   = projectdir("figures", parts...)
srcdir(parts...)   = projectdir("src", parts...)
tabledir(parts...) = projectdir("tables", parts...)

end # module

