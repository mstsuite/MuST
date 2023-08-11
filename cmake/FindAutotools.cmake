
find_program(AUTORECONF_EXECUTABLE
        NAMES autoreconf
        DOC "Autoreconf" REQUIRED)

find_program(AUTOCONF_EXECUTABLE
        NAMES autoconf
        DOC "Autoconf" REQUIRED)

find_program(AUTOMAKE_EXECUTABLE
        NAMES automake
        DOC "Automake" REQUIRED)

find_program(MAKE_EXECUTABLE
        NAMES gmake make
        NAMES_PER_DIR
        DOC "GNU Make")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Autotools
        REQUIRED_VARS AUTOCONF_EXECUTABLE AUTOMAKE_EXECUTABLE MAKE_EXECUTABLE)

mark_as_advanced(AUTOCONF_EXECUTABLE AUTOMAKE_EXECUTABLE MAKE_EXECUTABLE)
