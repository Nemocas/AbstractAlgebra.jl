function is_loaded_directly()
  try
    @static if VERSION < v"1.11.0-"
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
      Base.filter!(sf -> sf[1].func === :_tryrequire_from_serialized, bt)
      return length(bt) == 0 ||
            (length(bt) == 1 && length(only(bt)[1].linfo.specTypes.parameters) < 4)
    else
      # Starting with julia 1.11, the package loading was completely revamped.
      # The only difference in the callstack is the line number of the call to _include_from_serialized
      # inside of the _require_search_from_serialized function.
      # To make it a bit more robust, we check the difference between the line number of the beginning
      # of _require_search_from_serialized and the call to _include_from_serialized.
      # For `using OtherPackage`, the difference is 61, while for `using Package`, the difference is 75 or 78
      # (on all 1.11 pre-releases up to 1.11.0-rc1 and 1.12.0-DEV.896, which are the newest at the time of writing this).
      bt = Base.process_backtrace(Base.backtrace())
      Base.filter!(sf -> contains(string(sf[1].func), "_require_search_from_serialized"), bt)
      bt_entry = only(bt)[1]
      return bt_entry.line - bt_entry.linfo.def.line >= 70
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
