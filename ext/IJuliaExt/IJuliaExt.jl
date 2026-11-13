module IJuliaExt

using AbstractAlgebra
using IJulia

function __init__()
  # this needs to be called in __init__ as it modifies a global variable in another package
  # (AbstractAlgebra is considered a different package than this extension for precompilation purposes)
  AbstractAlgebra._set_ijulia_inited_hook!(() -> IJulia.inited)
end

end # module
