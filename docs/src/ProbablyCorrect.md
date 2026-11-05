```@meta
CurrentModule = AbstractAlgebra
CollapsedDocStrings = true
DocTestSetup = AbstractAlgebra.doctestsetup()
```

# Prototype

**This foundation is still prototypical: future UI/API changes are likely!**

_**We shall use the interim name "Dodgy Mode"**_


# Probably correct computations

AbstractAlgebra.jl provides the foundations upon which _fast and probably correct_ implementations
can be built (in a coherent manner).  The reason behind this is that there are some operations
where we can compute a _probably correct_ result much faster than a guaranteed correct result:
a classical example is primality testing, where the iterated Miller-Rabin method gives an almost
surely correct result must faster than the currently best-known, guaranteed algorithms.

This documentation comprises three parts: guidance for users, for developers, and for maintainers.


## USER DOC: _How to tell Oscar to use dodgy mode_

Unless you say otherwise, all user-facing OSCAR functions will use guaranteed-correct algorithms.
Some OSCAR functions offer a dodgy variant, which computes the same result with high probability
and taking less time.  The documentation for the function states whether it has a dodgy variant.

Reiteratng: dodgy mode permits some computations to be faster at the cost of
potentially produing an incorrect result (with low probability).  And by default
"dodgy mode" is disabled -- it must be explicitly enabled.  Dodgy mode can be enabled
either globally or individually for each applicable function family.  Here is a typical workflow:
```
    AbstractAlgebra.set_dodgy_mode(true)   # enable dodgy mode globally
    AbstractAlgebra.dodgy_steps_clear()
    ###  RUN YOUR COMPUTATION  ###
    AbstractAlgebra.dodgy_steps_get()      # list of possibly dodgy steps
    AbstractAlgebra.set_dodgy_mode(false)  # disable dodgy mode globally
```

### Enabling/Disabling dodgy mode

`set_dodgy_mode` is the function for enabling or disabling dodgy mode, it can
do this globally or for a specific function family.

For instance, we can enable dodgy mode for the `:is_prime` family by running:
```
    AbstractAlgebra.set_dodgy_mode(:is_prime, true);   # ENABLE
```
And to disable it for the same function family, we do
```
    AbstractAlgebra.set_dodgy_mode(:is_prime, false);  # DISABLE
```

If the function family is not specified then dodgy mode is enabled/disabled globally,
_i.e._ all function families are enabled/disabled.


### Log of dodgy steps

A log of potentially dodgy results is kept.  Each entry in the log indicates where a
potentially incorrect result was produced: the entry records  which function was called,
and the arguments there were passed.  If the log is empty then there were no dodgy steps,
so the result of your computation is guaranteed correct (as much as can be).

There are two functions associated to this log:
```
dodgy_steps_clear()
dodgy_steps_get()
```
The first deletes all entries in the log; the second returns the log as a Julia
`Vector{DodgyStepInfo}`.  Each `DodgyStepInfo` data-structure records the function name,
the values of the supplied arguments, and the returned value -- **the design of this
data-structure is likely to be revised!**

#### Interpreting and using the log

If the log is non-empty then one could, in principle, verify each dodgy result by some
independent means -- this would then certify that the final result is indeed correct.

Some functions with a dodgy variant make log entries in an "intelligent" manner.
For instance, if the function `is_prime` always gives a guaranteed result for small
arguments (_e.g._ below 1000000).  So if dodgy mode for the `:is_prime` is enabled
and `is_prime` is called with a small argument then there no log entry is recorded
since the result is surely correct.  Indeed, regardless of argument size, if `is_prime`
returns `false` then the result is surely correct; so no log entry will be made.


###  Example of use

**This example requires use of OSCAR** since the code here in `AbstractAlgebra` lays only the foundations.
We shall use the function `is_prime` in this example.
```
# First we disable dodgy mode globally:
AbstractAlgebra.set_dodgy_mode(false)

# Clear any pre-existing log entries:
AbstractAlgebra.dodgy_steps_clear()

AbstractAlgebra.dodgy_steps_get()  # empty list

# Set-up the values we shall use
p1 = 97;                    # "small" prime
p2 = 73742412689492826067;  # prime, below FLINT 80-bit bound
N = P1^20;                  # bigger than FLINT 80-bit bound, clearly not prime
p3 = next_prime(N);         # prime above FLINT 80-bit bound

# Call is_prime several times:
is_prime(p1)  # true
is_prime(p2)  # true
is_prime(N)   # false
is_prime(p3)  # true

# Since dodgy mode was disabled the log is empty
isempty(AbstractAlgebra.dodgy_steps_get())  || println("OOPS!");

# Now repeat the above calls but with dodgy mode enabled for function
# family `:is_prime` (includes both `is_prime` and `next_prime` functions)
AbstractAlgebra.set_dodgy_mode(:is_prime)  # also covers next_prime

is_prime(p1)  # true
isempty(AbstractAlgebra.dodgy_steps_get())  || println("OOPS!");

is_prime(p2)  # true
isempty(AbstractAlgebra.dodgy_steps_get())  || println("OOPS!");

is_prime(N)   # false
isempty(AbstractAlgebra.dodgy_steps_get())  || println("OOPS!");

is_prime(p3)  # true
(length(AbstractAlgebra.dodgy_steps_get()) == 1)  || println("OOPS!");

next_prime(p3);
(length(AbstractAlgebra.dodgy_steps_get()) == 2)  || println("OOPS!");

# Now we look at the dodgy steps recorded in the log.
# We expect one entry from the call `is_prime(p3)` and another one from `next_prime(p3)`:
AbstractAlgebra.dodgy_steps_get()
```



## DEVELOPER DOC: _How to modify a function to exploit dodgy mode_

We shall need the following function and macro:
```
   AbstractAlgebra.get_dodgy_mode(func_family_symbol)
   AbstractAlgebra.@RecordDodgyStep(fn_name, argv, result)
```
The first returns `true` iff dodgy mode for the specified family is enabled.
The second records a dodgy step, where
*   `fn_name::Symbol`   indicates the name of function which was called
*   `argv::Vector{Any}` are the argument values supplied
*   `result::Any`       is the esult produced

It is worthwhile looking at a concrete example: a simple example can be found in
the implementation of `is_prime(::ZZRingElem)` in the file  `Nemo/src/flint/fmpz.jl`.

There are two steps to complete:
0. decide a name for the dodgy function family
1. declare a dodgy function family by calling `add_dodgy_func` in the initialization code for the package
2. inside the function(s) you are modifying, insert a check whether dodgy mode is enabled for the function family, and if so, call the macro `RegsiterDodgyStep`

Here is a copy of the `is_prime` function as illustrative example:
```
   function is_prime(x::ZZRingElem)
     !is_probable_prime(x) && return false
     (x <= (Int128(1) << 80)) && return true  # guaranteed correct by FLINT doc
     if AbstractAlgebra.get_dodgy_mode(:is_prime)
       AbstractAlgebra.@RecordDodgyStep(:is_prime, (x,));
       return true
     end
     return Bool(@ccall libflint.fmpz_is_prime(x::Ref{ZZRingElem})::Cint)
   end
```

The function `is_probable_prime` is fast and probably correct.	If `is_probable_prime(x)`
returns `false` then we know the result is correct, so can return immediately (without
recording the call in the log of dodgy steps).  The FLINT documentation for `is_probable_prime`
states that the result is guaranteed when `x <= typemax(Int64)`; we check that condition
and return the guaranteed result.  Now we handle the case where `is_probable_prime` could
rarely return a wrong result.  If dodgy mode for the family `is_prime` is enabled, we
record the potentially dodgy step, and then return `true` (which is very probably correct).
If dodgy mode is not enabled, we use FLINT's potentially costly `fmpz_is_prime` function
to provide a guaranteed result.

Here is a simplistic but _less desirable_ dodgy variant of `is_prime`:
```
   function is_prime(x::ZZRingElem)
     AbstractAlgebra.@RecordDodgyStep(:is_prime, (x,));
     !is_probable_prime(x) && return false
     return Bool(@ccall libflint.fmpz_is_prime(x::Ref{ZZRingElem})::Cint)
   end
```
The disadvantage with this simplistic approach is that we record dodgy steps even
when we know (from the FLINT documentation) that the result is guaranteed correct.
This could fill the log with "false alarms".


## MAINTAINER DOC

The maintainer doc is in the source file: namely, `src/ProbablyCorrect.jl`
