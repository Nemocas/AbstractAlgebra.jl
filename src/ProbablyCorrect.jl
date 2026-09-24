# [[prototype]]
# Curt doc for "dodgy mode" in Oscar

# USER DOC: how to tell Oscar to use "dodgy mode"

# Dodgy mode permits some computations to be faster at the cost of
# potentially produing an incorrect result (with low probability).

# AbstractAlgebra.set_dodgy_mode(true)  -- returns previous setting
# AbstractAlgebra.dodgy_steps_clear()   -- clear log of dodgy steps
# ...RUN YOUR COMPUTATION...
# AbstractAlgebra.dodgy_steps_get()     -- list of potentially dodgy steps

# If `dodgy_steps_get` returns an empty vector then there were no dodgy steps!
# If the vector is not empty then there were some dodgy steps, so the results
# produced by those steps are potentially incorrect (though this is typically
# extremely unlikely)


# USER DOC: how to modify a function to exploit dodgy mode

# Look at the implementation of is_prime in Nemo/src/flint/fmpz.jl near line 1970

# Check whether dodgy mode is active: AbstractAlgebra.get_dodgy_mode();
# this returns `true` if dodgy mode is active; otherwise `false`.

# AbstractAlgebra.@RegisterDodgyStep(fn_name, argv)  where fn_name::Symbol and argv::Vector{Any}
# records that there was a call to `fn_name` with arguments `argv` where the result may
# potentially be incorrect.


# DESIGN/MAINTAINER DOC

# We use 3 global variables:
#   GLOBAL_VARIABLE_DodgyMode::Bool
#   -  `true` iff dodgy mode is active
#   GLOBAL_VARIABLE_DodgySteps::Vector{DodgyStepInfo}
#   -  contains a record of dodgy steps executed
#   - must be cleared explicitly
#   GLOBAL_VARIABLE_DodgySteps_MaxSize::Int
#   - upper limit for length of GLOBAL_VARIABLE_DodgySteps

#  setter/getter for GLOBAL_VARIABLE_DodgyMode
#  - get_dodgy_mode()
#  - set_dodgy_mode(b::Bool)   [[returns previous setting]]

# setter/getter for GLOBAL_VARIABLE_DodgySteps
#  - dodgy_steps_clear()    [[clears the vector]]
#  - dodgy_steps_get()
#  - @RegisterDodgyStep(fn_name, argv)
#    [[does nothing if the length limit has been reached]]



# Global stuff for the "probably correct" prototype.

# Fairly KISS: 
GLOBAL_VARIABLE_ListOfDodgyFuncs::Vector{Symbol} = Symbol[];
GLOBAL_VARIABLE_ActiveDodgyFuncs::Vector{Symbol} = Symbol[];


# This is for "registering" a dodgy function:
# it would normally be called in the initialization of the package.
function add_dodgy_func(func::Symbol)
  !(func in GLOBAL_VARIABLE_ListOfDodgyFuncs) && push!(GLOBAL_VARIABLE_ListOfDodgyFuncs, func);
  return nothing
end


function get_active_dodgy_funcs()
  global GLOBAL_VARIABLE_ActiveDodgyFuncs;
  return GLOBAL_VARIABLE_ActiveDodgyFuncs;
end;


# Activate (or deactivate) **all** dodgy functions
function set_dodgy_mode(b::Bool = true)
  global GLOBAL_VARIABLE_ListOfDodgyFuncs;
  global GLOBAL_VARIABLE_ActiveDodgyFuncs;
  if !b
    GLOBAL_VARIABLE_ActiveDodgyFuncs = Symbol[];
  else
    GLOBAL_VARIABLE_ActiveDodgyFuncs = GLOBAL_VARIABLE_ListOfDodgyFuncs;
  end
  return GLOBAL_VARIABLE_ActiveDodgyFuncs;
end;

function set_dodgy_mode(func::Symbol, b::Bool = true)
  global GLOBAL_VARIABLE_ListOfDodgyFuncs;
  global GLOBAL_VARIABLE_ActiveDodgyFuncs;
  !(func in GLOBAL_VARIABLE_ListOfDodgyFuncs) && error("Unregistered function ID");
  if b
    !(func in GLOBAL_VARIABLE_ActiveDodgyFuncs) && push!(GLOBAL_VARIABLE_ActiveDodgyFuncs, func);
  else
    GLOBAL_VARIABLE_ActiveDodgyFuncs = filter(!=(func), GLOBAL_VARIABLE_ActiveDodgyFuncs);
  end;
  return GLOBAL_VARIABLE_ActiveDodgyFuncs;  ## or just nothing?
end;

function get_dodgy_mode(func::Symbol)
  global GLOBAL_VARIABLE_ListOfDodgyFuncs;
  global GLOBAL_VARIABLE_ActiveDodgyFuncs;
  !(func in GLOBAL_VARIABLE_ListOfDodgyFuncs) && error("Unregistered function ID");
  return (func in GLOBAL_VARIABLE_ActiveDodgyFuncs);
end;



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
  print(io, "DodgyStepInfo(fn_name=$(info.fn_name), argv=$(info.argv))");
end


GLOBAL_VARIABLE_DodgySteps::Vector{DodgyStepInfo} = DodgyStepInfo[];
GLOBAL_VARIABLE_DodgySteps_MaxSize::Int = 10;  #set this to 1 if you just want to know whether any (potentially) dodgy result was generated
const GLOBAL_VARIABLE_EmptyStackTrace = Union{Ptr{Nothing}, Base.InterpreterIP}[];


function dodgy_steps_clear()
  global GLOBAL_VARIABLE_DodgySteps;
  empty!(GLOBAL_VARIABLE_DodgySteps);
  return nothing;
end;


function dodgy_steps_get()
  global GLOBAL_VARIABLE_DodgySteps;
  return GLOBAL_VARIABLE_DodgySteps;
end;


# Use macro RegisterDodgyStep in code which wants to register a potentially dodgy step
# The macro uses the following global variables:
#   AbstractAlgebra.GLOBAL_VARIABLE_DodgySteps
#   AbstractAlgebra.GLOBAL_VARIABLE_DodgySteps_MaxSize

macro RegisterDodgyStep(fn_name, argv)  # fn_name::Symbol;  argv::Vector{Any}
  SourceLine = __source__.line;
  SourceFile = String(__source__.file);  # In Julia 1.13 there will be @__FUNCTION__ (I have julia 1.10)
  return :(
    if length(AbstractAlgebra.GLOBAL_VARIABLE_DodgySteps) < GLOBAL_VARIABLE_DodgySteps_MaxSize

      push!(AbstractAlgebra.GLOBAL_VARIABLE_DodgySteps,
            AbstractAlgebra.DodgyStepInfo($(esc(fn_name)),
                                          ($(esc(SourceFile)), $(esc(SourceLine))), 
                                          $(esc(argv)),
                                          nothing,
                                          AbstractAlgebra.GLOBAL_VARIABLE_EmptyStackTrace));
    end
  )  # end of quoted part
end

