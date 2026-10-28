# [[prototype]]
# Curt doc for "dodgy mode" in Oscar

# **USER DOC**
# _How to tell Oscar to use dodgy mode_

# Dodgy mode permits some computations to be faster at the cost of
# potentially produing an incorrect result (with low probability).
# By default "dodgy mode" is inactive -- it must be explicitly enabled:
# this can be done either globally or individually for each applicable
# function.  Here is a typical workflow:

#    AbstractAlgebra.set_dodgy_mode(true);
#    AbstractAlgebra.dodgy_steps_clear()
#    ...RUN YOUR COMPUTATION...
#    AbstractAlgebra.dodgy_steps_get()     -- list of possibly dodgy steps
#    AbstractAlgebra.set_dodgy_mode(false);

# If `dodgy_steps_get()` returns an empty vector then there were no dodgy
# steps, so all results computed are reliable.
# If the vector is not empty then there were some dodgy steps, so the results
# produced by those steps are potentially incorrect (though this is typically
# extremely unlikely).

# Dodgy mode for a particular function may be enabled/disabled individually.
# For instance, we can set dodgy more for the ":is_prime" family by running:

#    AbstractAlgebra.set_dodgy_mode(:is_prime, true);   # ENABLE
#    AbstractAlgebra.set_dodgy_mode(:is_prime, false);  # DISABLE


# **USER/DEVELOPER DOC**
# _How to write/modify a function to exploit dodgy mode_

# A simple example can be found in the implementation of is_prime(::ZZRingElem)
# in the file  Nemo/src/flint/fmpz.jl  near line 1960.

# Check whether dodgy mode is active for the corresponding family, by calling
#    AbstractAlgebra.get_dodgy_mode(family_symbol)
# This returns `true` if dodgy mode is active for that family; otherwise `false`.

# If dodgy mode is active then a possibly dodgy step must be registered by calling
#    AbstractAlgebra.@RegisterDodgyStep(fn_name, argv)
# where
#    fn_name::Symbol   indicates the name of function which was called
#    argv::Vector{Any} are the argument values supplied.
# If the argument values are "simple", and it is known that for these values
# the result will surely be correct, then @RegisterDodgyStep should not be called.
# For instance, is_prime is always correct for values below 1000000, so even if
# dodgy_mode is enabled, no dodgy step should be registered when the given value
# is below 1000000 -- actually, the true safe limit for is_prime is rather higher.


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
#  -  dodgy_steps_clear()    [[clears the vector/register]]
#  -  dodgy_steps_get()      [[value of GLOBAL_VARIABLE_DodgySteps]
#  -  @RegisterDodgyStep(fn_name, argv)
#     append a new dodgy step to the register
#     [[does nothing if the length limit has been reached]]



# Global stuff for the "probably correct" prototype.

# Initialization (fairly KISS):
GLOBAL_VARIABLE_ListOfDodgyFuncs::Vector{Symbol} = Symbol[]   # see add_dodgy_func below
GLOBAL_VARIABLE_ActiveDodgyFuncs::Vector{Symbol} = Symbol[]


# This is for "registering" a dodgy function:
# it would normally be called in the initialization of the package.
function add_dodgy_func(func::Symbol)
  !(func in GLOBAL_VARIABLE_ListOfDodgyFuncs) && push!(GLOBAL_VARIABLE_ListOfDodgyFuncs, func)
  return nothing
end


function get_active_dodgy_funcs()
  global GLOBAL_VARIABLE_ActiveDodgyFuncs
  return GLOBAL_VARIABLE_ActiveDodgyFuncs
end


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

function set_dodgy_mode(func::Symbol, b::Bool = true)
  global GLOBAL_VARIABLE_ListOfDodgyFuncs
  global GLOBAL_VARIABLE_ActiveDodgyFuncs
  !(func in GLOBAL_VARIABLE_ListOfDodgyFuncs) && error("Unregistered function ID")
  if b
    !(func in GLOBAL_VARIABLE_ActiveDodgyFuncs) && push!(GLOBAL_VARIABLE_ActiveDodgyFuncs, func)
  else
    GLOBAL_VARIABLE_ActiveDodgyFuncs = filter(!=(func), GLOBAL_VARIABLE_ActiveDodgyFuncs)
  end
  return GLOBAL_VARIABLE_ActiveDodgyFuncs  ## or just nothing?
end

function get_dodgy_mode(func::Symbol)
  global GLOBAL_VARIABLE_ListOfDodgyFuncs
  global GLOBAL_VARIABLE_ActiveDodgyFuncs
  !(func in GLOBAL_VARIABLE_ListOfDodgyFuncs) && error("Unregistered function ID")
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
  resv::Any  # either simple value or a tuple
  stack::Vector{Union{Ptr{Nothing}, Base.InterpreterIP}}

  #  CALLED ONLY BY macro RegisterDodgyStep (see below)
  function DodgyStepInfo(fn_name::Symbol,
                         SourceLocation::Tuple{String,Int},
                         argv::Tuple,
                         resv::Any,
                         stack)
    return new(fn_name, SourceLocation, argv, resv, stack)
  end
end

# KISS for now:
function Base.show(io::IO, info::DodgyStepInfo)
  print(io, "DodgyStepInfo(fn_name=$(info.fn_name), argv=$(info.argv))")
end


GLOBAL_VARIABLE_DodgySteps::Vector{DodgyStepInfo} = DodgyStepInfo[]
GLOBAL_VARIABLE_DodgySteps_MaxSize::Int = 10  #set this to 1 if you just want to know whether any (potentially) dodgy result was generated [[BUG: not currently user-settable]]
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


# Use macro RegisterDodgyStep in code which wants to register a potentially dodgy step
# The macro uses the following global variables:
#   AbstractAlgebra.GLOBAL_VARIABLE_DodgySteps
#   AbstractAlgebra.GLOBAL_VARIABLE_DodgySteps_MaxSize

macro RegisterDodgyStep(fn_name, argv)  # fn_name::Symbol;  argv::Vector{Any}
  SourceLine = __source__.line
  SourceFile = String(__source__.file)  # In Julia 1.13 there will be @__FUNCTION__ (I have julia 1.10)
  return :(
    if length(AbstractAlgebra.GLOBAL_VARIABLE_DodgySteps) < GLOBAL_VARIABLE_DodgySteps_MaxSize

      push!(AbstractAlgebra.GLOBAL_VARIABLE_DodgySteps,
            AbstractAlgebra.DodgyStepInfo($(esc(fn_name)),
                                          ($(esc(SourceFile)), $(esc(SourceLine))),
                                          $(esc(argv)),
                                          nothing,
                                          AbstractAlgebra.GLOBAL_VARIABLE_EmptyStackTrace))
    end
  )  # end of quoted part
end
