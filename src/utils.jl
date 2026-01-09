function is_loaded_directly()
  try
    @static if VERSION < v"1.11.0-"
      @debug "is_loaded_directly: VERSION < 1.11.0-"
      # Check if were loaded from another package
      # if VERSION < 1.7.*, only the "other" package will have the
      # _tryrequire_from_serialized in the backtrace.
      # if VERSION >= 1.8, also doing 'using Package' will have
      # _tryrequire_from_serialized the backtrace.
      #
      # To still distinguish both scenarios, notice that
      # 'using OtherPackage' will either have _tryrequire_from_serialized at least twice,
      # or one with four arguments (hence five as the function name is the first argument)
      # 'using Package' serialized will have a version with less arguments
      bt = Base.process_backtrace(Base.backtrace())
      @debug "is_loaded_directly: full backtrace:\n$(sprint(show, "text/plain", bt))"
      Base.filter!(sf -> sf[1].func === :_tryrequire_from_serialized, bt)
      length_bt = length(bt)
      @debug "is_loaded_directly: `_tryrequire_from_serialized` appears $(length_bt) times in backtrace"
      length_bt == 0 && return true
      length_bt != 1 && return false
      params = only(bt)[1].linfo.specTypes.parameters
      @debug "is_loaded_directly: `_tryrequire_from_serialized` gets called with parameters of type $(params)"
      return length(params) < 4
    else
      @debug "is_loaded_directly: VERSION >= 1.11.0-"
      # Starting with julia 1.11, the package loading was completely revamped.
      # The only difference in the callstack is the line number of the call to _include_from_serialized
      # inside of the _require_search_from_serialized function.
      # To make it a bit more robust, we check the difference between the line number of the beginning
      # of _require_search_from_serialized and the call to _include_from_serialized.
      st = Base.stacktrace(Base.backtrace())
      @debug "is_loaded_directly: full backtrace:\n$(sprint(show, "text/plain", st))"
      Base.filter!(sf -> contains(string(sf.func), "_require_search_from_serialized") && !startswith(string(sf), "kwcall(") && !contains(string(sf), "[inlined]"), st)
      length_st = length(st)
      @debug "is_loaded_directly: `_require_search_from_serialized` appears $(length_st) times in backtrace; expected 1"
      sf = only(st)
      line_call = sf.line
      line_funcbegin = (sf.linfo::Core.MethodInstance).def.line
      line_difference = line_call - line_funcbegin
      @debug "is_loaded_directly: `_require_search_from_serialized` called at line $line_call, function begins at line $line_funcbegin, difference $(line_difference)"
      # difference for `using Package` / `using OtherPackage`
      # 1.11.0-alpha1:    75 /  61
      # 1.11.0-alpha2:    75 /  61
      # 1.11.0-beta1:     75 /  61
      # 1.11.0-beta2:     75 /  61
      # 1.11.0-rc1:       75 /  61
      # 1.11.0-rc2:       78 /  61
      # 1.11.0-rc3:       78 /  61
      # 1.11.0-rc4:       77 /  61
      # 1.11.0:           77 /  61
      # 1.11.1:           77 /  61
      # 1.12.0-DEV.896:   78 /  61 # ignored
      # 1.12.0-DEV.1322:  88 /  72 # ignored
      # 1.12.0-DEV.1506: 106 /  93
      @static if v"1.11.0-" < VERSION < v"1.12.0-DEV.1506"
        return line_difference >= 73
      else # v"1.12.0-DEV.1506" <= VERSION
        return line_difference >= 100
      end
    end
  catch e
    @debug "Error while checking if loaded directly" exception=(e, Base.catch_backtrace())
    return true
  end
end

function should_show_banner()
  # Respect the -q flag
  isquiet = Bool(Base.JLOptions().quiet)
  show_banner = !isquiet && is_loaded_directly() && isinteractive() && Base.JLOptions().banner != 0
  return show_banner
end
