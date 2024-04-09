function should_show_banner()
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
    isinteractive_manual =
      length(bt) == 0 || (length(bt) == 1 && length(only(bt)[1].linfo.specTypes.parameters) < 4)

    # Respect the -q flag
    isquiet = Bool(Base.JLOptions().quiet)
    show_banner = !isquiet && isinteractive_manual && isinteractive() &&
        Base.JLOptions().banner != 0
    return show_banner
end
