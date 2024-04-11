module ModB

import AbstractAlgebra: should_show_banner

using ModA

function __init__()
    if should_show_banner()
        println("Banner of ModB")
    end
end

end # module ModB
