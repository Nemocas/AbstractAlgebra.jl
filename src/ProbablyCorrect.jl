# [[prototype]]
# Curt doc for "dodgy mode" in Oscar

# DESIGN/MAINTAINER DOC

# We use the following global variables:
#   GLOBAL_VARIABLE_ListOfDodgyFuncs::Vector{Symbol}
#   -  a list of symbols: each symbol corresponds to a
#      family of functions which offer a dodgy mode

#   GLOBAL_VARIABLE_ActiveDodgyFuncs::Vector{Symbol}
#   -  a list of symbols: each symbol corresponds to a
#      family of functions for which dodgy mode is enabled

#   GLOBAL_VARIABLE_DodgySteps::Vector{DodgyStepInfo}
#   -  contains a record of dodgy steps executed
#   -  must be cleared explicitly

#   GLOBAL_VARIABLE_DodgySteps_MaxSize::Int
#   - upper limit for length of GLOBAL_VARIABLE_DodgySteps


# The values of the global variables should be accessed or
# changed only via the relevant setter/getter functions:

#  setters/getters for GLOBAL_VARIABLE_ActiveDodgyFuncs
#  -  get_dodgy_mode(func::Symbol)
#     Returns `true` iff dodgy mode is enabled for family  `func`
#  -  set_dodgy_mode(func::Symbol, b::Bool)
#     Enable dodgy mode for family `func` if `b` is `true`,
#     otherwise disable dodgy mode for that family.
#     Returns new value of GLOBAL_VARIABLE_ListOfDodgyFuncs.
#  -  set_dodgy_mode(b::Bool)
#     enable/disable dodgy mode for all families.

# setter/getter for GLOBAL_VARIABLE_DodgySteps
#  -  dodgy_steps_clear()    [[clears the vector/log]]
#  -  dodgy_steps_get()      [[value of GLOBAL_VARIABLE_DodgySteps]
#  -  @RecordDodgyStep(fn_name, argv, result)
#     append a new dodgy step to the log
#     [[does nothing if the length limit has been reached]]



# Global stuff for the "probably correct" prototype.

# Initialization (fairly KISS):
GLOBAL_VARIABLE_ListOfDodgyFuncs::Vector{Symbol} = Symbol[]   # see add_dodgy_func below
GLOBAL_VARIABLE_ActiveDodgyFuncs::Vector{Symbol} = Symbol[]


@doc raw"""
    add_dodgy_func(FuncFamily::Symbol)

This registers the name `FuncFamily` as a name for a function family having dodgy variants.
It would normally be called in the initialization of the package.
"""
function add_dodgy_func(func::Symbol)
  !(func in GLOBAL_VARIABLE_ListOfDodgyFuncs) && push!(GLOBAL_VARIABLE_ListOfDodgyFuncs, func)
  return nothing
end


@doc raw"""
    get_dodgy_funcs()

Returns a list of symbols for dodgy function families.
"""
function get_active_dodgy_funcs()
  global GLOBAL_VARIABLE_ActiveDodgyFuncs
  return GLOBAL_VARIABLE_ActiveDodgyFuncs
end


@doc raw"""
    set_dodgy_mode(b::Bool = true)
    set_dodgy_mode(func:Symbol, b::Bool = true)

Enable dodgy mode: globally (1st version), for the given family (2nd version).
"""
function set_dodgy_mode end

# Activate (or deactivate) **all** dodgy functions
function set_dodgy_mode(b::Bool = true)
  global GLOBAL_VARIABLE_ListOfDodgyFuncs
  global GLOBAL_VARIABLE_ActiveDodgyFuncs
  if !b
    GLOBAL_VARIABLE_ActiveDodgyFuncs = Symbol[]
  else
    GLOBAL_VARIABLE_ActiveDodgyFuncs = GLOBAL_VARIABLE_ListOfDodgyFuncs
  end
  return GLOBAL_VARIABLE_ActiveDodgyFuncs
end

function set_dodgy_mode(FuncFamily::Symbol, b::Bool = true)
  global GLOBAL_VARIABLE_ListOfDodgyFuncs
  global GLOBAL_VARIABLE_ActiveDodgyFuncs
  !(func in GLOBAL_VARIABLE_ListOfDodgyFuncs) && error("Unregistered function family ID")
  if b
    !(func in GLOBAL_VARIABLE_ActiveDodgyFuncs) && push!(GLOBAL_VARIABLE_ActiveDodgyFuncs, func)
  else
    GLOBAL_VARIABLE_ActiveDodgyFuncs = filter(!=(func), GLOBAL_VARIABLE_ActiveDodgyFuncs)
  end
  return GLOBAL_VARIABLE_ActiveDodgyFuncs  ## or just nothing?
end

@doc raw"""
    get_dodgy_mode(func:Symbol)

Return `true` iff dodgy mode is enable for family `func`.
"""
function get_dodgy_mode(FuncFamily::Symbol)
  global GLOBAL_VARIABLE_ListOfDodgyFuncs
  global GLOBAL_VARIABLE_ActiveDodgyFuncs
  !(func in GLOBAL_VARIABLE_ListOfDodgyFuncs) && error("Unregistered function family ID")
  return (func in GLOBAL_VARIABLE_ActiveDodgyFuncs)
end



# Global variable DodgySteps is just a list of DodgyStepInfo structs
# (KISS goodbye) The struct now contains several fields, but the protoype
# currently uses only some of them.
struct DodgyStepInfo
# fn_family::Symbol    the "marker" used to identify the function family
  fn_name::Symbol
  #  methodID::Union{Nothing,Method}
  SourceLocation::Tuple{String,Int}
  argv::Tuple
  result::Any  # either simple value or a tuple
  stack::Vector{Union{Ptr{Nothing}, Base.InterpreterIP}}

  #  CALLED ONLY BY macro RecordDodgyStep (see below)
  function DodgyStepInfo(fn_name::Symbol,
                         SourceLocation::Tuple{String,Int},
                         argv::Tuple,
                         result::Any,
                         stack)
    return new(fn_name, SourceLocation, argv, result, stack)
  end
end

# KISS for now:
function Base.show(io::IO, info::DodgyStepInfo)
  print(io, "DodgyStepInfo(fn_name=$(info.fn_name), argv=$(info.argv))")
end


GLOBAL_VARIABLE_DodgySteps::Vector{DodgyStepInfo} = DodgyStepInfo[]
GLOBAL_VARIABLE_DodgySteps_MaxSize::Int = 20  #set this to 1 if you just want to know whether any (potentially) dodgy result was generated [[BUG: not currently user-settable]]
const GLOBAL_VARIABLE_EmptyStackTrace = Union{Ptr{Nothing}, Base.InterpreterIP}[]


function dodgy_steps_clear()
  global GLOBAL_VARIABLE_DodgySteps
  empty!(GLOBAL_VARIABLE_DodgySteps)
  return nothing
end


function dodgy_steps_get()
  global GLOBAL_VARIABLE_DodgySteps
  return GLOBAL_VARIABLE_DodgySteps
end


# Use macro RecordDodgyStep in code which wants to record a potentially dodgy step
# The macro uses the following global variables:
#   AbstractAlgebra.GLOBAL_VARIABLE_DodgySteps
#   AbstractAlgebra.GLOBAL_VARIABLE_DodgySteps_MaxSize

macro RecordDodgyStep(fn_name, argv, result)  # fn_name::Symbol;  argv::Vector{Any};  result::Any
  SourceLine = __source__.line
  SourceFile = String(__source__.file)  # In Julia 1.13 there will be @__FUNCTION__ (I have julia 1.10)
  return :(
    if length(AbstractAlgebra.GLOBAL_VARIABLE_DodgySteps) < GLOBAL_VARIABLE_DodgySteps_MaxSize

      push!(AbstractAlgebra.GLOBAL_VARIABLE_DodgySteps,
            AbstractAlgebra.DodgyStepInfo($(esc(fn_name)),
                                          ($(esc(SourceFile)), $(esc(SourceLine))),
                                          $(esc(argv)),
                                          $(esc(result)),
                                          AbstractAlgebra.GLOBAL_VARIABLE_EmptyStackTrace))
    end
  )  # end of quoted part
end
