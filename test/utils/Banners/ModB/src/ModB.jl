module ModB

import AbstractAlgebra: should_show_banner

using ModA

function __init__()
    @debug "__init__ of ModB"
    if should_show_banner()
        println("Banner of ModB")
    end
end

end # module ModB
