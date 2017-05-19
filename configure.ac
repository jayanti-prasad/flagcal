#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
#AC_PROG_LIBTOOL
#AC_DISABLE_SHARED
AC_INIT([flagcal], [1.0], [jayanti@iucaa.ernet.in])
AC_PREREQ([2.68])
AM_INIT_AUTOMAKE([1.10 no-define])
AC_CONFIG_SRCDIR([src/io.c])
AC_CONFIG_HEADERS([config.h])

AC_OPENMP

# Checks for programs.
AC_PROG_CC

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lX11':
AC_CHECK_LIB([X11], [main])
# FIXME: Replace `main' with a function in `-lcfitsio':
AC_CHECK_LIB([cfitsio], [ffvers])
# FIXME: Replace `main' with a function in `-lgfortran':
AC_CHECK_LIB([gfortran], [main])
# FIXME: Replace `main' with a function in `-lpng':
AC_CHECK_LIB([png], [main])
# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m], [main])
# FIXME: Replace `main' with a function in `-lpgplot':
AC_CHECK_LIB([pgplot], [main])
# FIXME: Replace `main' with a function in `-lcpgplot':
AC_CHECK_LIB([cpgplot], [cpgbeg])

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h sys/time.h])

AC_CHECK_HEADER([fitsio.h],
               [AC_DEFINE([HAVE_FITSIO_H], [1],
               [Define to 1 if you have <fitsio.h>.])],
               [AC_MSG_ERROR([sorry, cant do anything for you])])

AC_CHECK_HEADER([omp.h],
               [AC_DEFINE([HAVE_OMP_H], [1],
               [Define to 1 if you have <omp.h>.])],
               [AC_MSG_ERROR([sorry no OpenMP :(, cant do anything for you])])

AC_CHECK_HEADER([cpgplot.h],
               [AC_DEFINE([HAVE_CPGPLOT_H], [1],
               [Define to 1 if you have <cpgplot.h>.])],
               [AC_MSG_ERROR([sorry no cpgplot.h :(, cant do anything for you])])


# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([floor ftime gettimeofday pow sqrt strchr strstr])

AC_CONFIG_FILES([Makefile
                 src/Makefile])
AC_OUTPUT
