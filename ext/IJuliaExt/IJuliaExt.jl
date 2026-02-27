module IJuliaExt

using AbstractAlgebra
using IJulia

AbstractAlgebra._set_ijulia_inited_hook!(() -> IJulia.inited)

end # module
