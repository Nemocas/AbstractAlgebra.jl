module ModA

import AbstractAlgebra: should_show_banner

function __init__()
    @debug "__init__ of ModA"
    if should_show_banner()
        println("Banner of ModA")
    end
end

end # module ModA
