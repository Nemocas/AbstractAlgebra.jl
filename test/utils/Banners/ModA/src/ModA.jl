module ModA

import AbstractAlgebra: should_show_banner

function __init__()
    if should_show_banner()
        println("Banner of ModA")
    end
end

end # module ModA
