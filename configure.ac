#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT(bamkit, 1.0, mingfu.shao@cs.cmu.edu)
AM_INIT_AUTOMAKE()
AC_CONFIG_SRCDIR([src/main.cc])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

### require environment variables BOOST_HOME 
AC_ARG_WITH(htslib, AS_HELP_STRING([--with-htslib], [home directory for htslib]), HTSLIB_HOME=$withval, HTSLIB_HOME=)
AS_IF([test "x$HTSLIB_HOME" != "x"], [AC_SUBST([CXXFLAGS], ["-I$HTSLIB_HOME/include $CXXFLAGS"])], [])
AS_IF([test "x$HTSLIB_HOME" != "x"], [AC_SUBST([LDFLAGS], ["-L$HTSLIB_HOME/lib $LDFLAGS"])])
AC_CHECK_LIB([hts], [hts_open], [], [AC_MSG_ERROR([could not find htslib])])

# Checks for header files.
AC_CHECK_HEADERS([stdint.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_UINT32_T

# Checks for library functions.

AC_CONFIG_FILES([Makefile
                 src/Makefile])
AC_OUTPUT
