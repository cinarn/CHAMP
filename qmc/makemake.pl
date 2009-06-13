#! /usr/bin/perl
#
# Usage: makemake {<program name> {<F95 compiler or fc or f77 or cc or c>}}
#
# Generate a Makefile from the sources in the current directory.  The source
# files may be in either C, FORTRAN 77, Fortran 90 or some combination of
# these languages.  If the F95 compiler specified is cray or parasoft, then
# the Makefile generated will conform to the conventions of these compilers.
# To run makemake, it will be necessary to modify the first line of this script
# to point to the actual location of Perl on your system.
#
# Written by Michael Wester <wester@math.unm.edu> February 16, 1995
# Cotopaxi (Consulting), Albuquerque, New Mexico
#
open(MAKEFILE, "> Makefile");
#
print MAKEFILE "# Warning: this 'Makefile' is automatically generated by typing 'make make'.\n";
print MAKEFILE "# Warning: DO NOT MODIFY THIS 'Makefile' BY HAND!\n";
print MAKEFILE "# Warning: However, you may need to modify the file 'makefile.inc' for machine-specific options.\n\n";

print MAKEFILE "include ../makefile.inc\n\n";
print MAKEFILE "PROG =\t\$(CHAMP_EXE)\n\n";
#
# Source listing
#
print MAKEFILE "SRCS =\t";
@srcs = <commons/*.f90 *.f90 *.f95 *.f *.F *.c fit/*.f fit/MPI/*f vmc/*.f vmc/*.f90 vmc/MPI/*.f dmc/*.f dmc/dmc_elec/*.f dmc/dmc_elec/MPI/*.f dmc/dmc_elec/MPI_global_pop/*.f dmc/dmc_elec/MPI_global_pop_big/*.f>;
&PrintWords(8, 0, @srcs);
print MAKEFILE "\n\n";
#
# Object listing
#
print MAKEFILE "OBJS =\t";
@objs = @srcs;
foreach (@objs) { s/\.[^.]+$/.o/ };
&PrintWords(8, 0, @objs);
print MAKEFILE "\n\n";
#
# Define common macros
#
# make
#

print MAKEFILE "# PHONY is a Make keyword that tells it to execute the rule regardless of\n";
print MAKEFILE "# the dependencies.  We use it here for 2 reasons:\n";
print MAKEFILE "# 1) prevents it from getting confused if there happens to be a file named\n";
print MAKEFILE "#    'clean' or 'clean_all' in the directory.\n";
print MAKEFILE "# 2) Also, for the libraries, linpack etc. it does the make even though there\n";
print MAKEFILE "#    are directories named linpack etc.\n";
print MAKEFILE ".PHONY: clean_local clean clean_all clean_all_lib make\n\n";

#print MAKEFILE "ifdef EINSPLINE\n";
#print MAKEFILE "LIBS = ../lib/lib/libcyrus.a ../lib/lib2/blas/libblas.a ../lib/lib2/lapack/liblapack.a ../lib/lib2/linpack/liblinpack.a ../lib/lib2/einspline/lib/libeinspline.a ../lib/lib2/pspline/pspline/libpspline.a ../lib/SimulatedAnnealing/quench_anneal/lib/libquench.a ../lib/SimulatedAnnealing/quench_anneal/lib/libquench_seq.a\n";
#print MAKEFILE "else\n";
#print MAKEFILE "LIBS = ../lib/lib/libcyrus.a ../lib/lib2/blas/libblas.a ../lib/lib2/lapack/liblapack.a ../lib/lib2/linpack/liblinpack.a ../lib/lib2/pspline/pspline/libpspline.a ../lib/SimulatedAnnealing/quench_anneal/lib/libquench.a ../lib/SimulatedAnnealing/quench_anneal/lib/libquench_seq.a\n";
#print MAKEFILE "endif\n\n";

print MAKEFILE "\$(PROG): \$(CHAMP_LIBS_MAKE) \$(OBJS)\n";
print MAKEFILE "\t\$(CHAMP_LD) \$(CHAMP_LD_FLAGS) -o \$@ \$(OBJS) \$(CHAMP_LIBS) \$(CHAMP_LD_END)\n\n";
print MAKEFILE "";
print MAKEFILE "../lib/lib/libcyrus.a:\n";
print MAKEFILE "\tcd ../lib ; make\n";
print MAKEFILE "\n";
print MAKEFILE "../lib/lib2/blas/libblas.a:\n";
print MAKEFILE "\tcd ../lib ; make\n";
print MAKEFILE "\n";
print MAKEFILE "../lib/lib2/lapack/liblapack.a:\n";
print MAKEFILE "\tcd ../lib ; make\n";
print MAKEFILE "\n";
print MAKEFILE "../lib/lib2/linpack/liblinpack.a:\n";
print MAKEFILE "\tcd ../lib ; make\n";
print MAKEFILE "\n";
print MAKEFILE "../lib/lib2/einspline/lib/libeinspline.a:\n";
print MAKEFILE "\tcd ../lib ; make\n";
print MAKEFILE "\n";
print MAKEFILE "../lib/SimulatedAnnealing/quench_anneal/lib/libquench.a: ../lib/SimulatedAnnealing/quench_anneal/lib/libquench.a(main/anneal.o)\n";
print MAKEFILE "\tcd ../lib ; make\n";
print MAKEFILE "\n";
#
# make various cleans
#
print MAKEFILE "clean_local:\n";
print MAKEFILE "\trm -f *.o *.dif *.lst *.sav *.mod\n";
print MAKEFILE "\n";
print MAKEFILE "clean:\n";
print MAKEFILE "\trm -f \$(OBJS) *.dif *.lst *.sav *.mod\n";
print MAKEFILE "\n";
print MAKEFILE "clean_all:\n";
print MAKEFILE "\tmake clean\n";
print MAKEFILE "\trm -f *exe\n";
print MAKEFILE "\n";
print MAKEFILE "clean_all_lib:\n";
print MAKEFILE "\tmake clean_all\n";
print MAKEFILE "\tcd ../lib ; make clean_all\n";

#
# make make
#
print MAKEFILE "make:\n";
print MAKEFILE "\tperl makemake.pl $ARGV[0]\n";
print MAKEFILE "\t/usr/bin/ctags *.f *.f90\n";
#for etags want all .f .f90 files cataloged
print MAKEFILE "\tfind . -iregex '.*\\.f\\(90\\)?' | etags -\n";
#
# libraries
#print MAKEFILE "../lib/lib/libcyrus.a:\n";
#print MAKEFILE "../lib/lib/libcyrus.a:\n";
#print MAKEFILE "\tcd ../lib ; make\n";
#print MAKEFILE "../lib/lib2/blas/libblas.a:\n";
#print MAKEFILE "\tcd ../lib ; make\n";
#print MAKEFILE "../lib/lib2/lapack/liblapack.a:\n";
#print MAKEFILE "\tcd ../lib ; make\n";
#print MAKEFILE "../lib/lib2/linpack/liblinpack.a:\n";
#print MAKEFILE "\tcd ../lib ; make\n";
#print MAKEFILE "../lib/lib2/einspline/lib/libeinspline.a:\n";
#print MAKEFILE "\tcd ../lib ; make\n";
#print MAKEFILE "../lib/lib2/pspline/pspline/libpspline.a:\n";
#print MAKEFILE "\tcd ../lib ; make\n";
#print MAKEFILE "../lib/SimulatedAnnealing/quench_anneal/lib/libquench.a:\n";
#print MAKEFILE "\tcd ../lib ; make\n";
#print MAKEFILE "../lib/SimulatedAnnealing/quench_anneal/lib/libquench_seq.a:\n";
#print MAKEFILE "\tcd ../lib ; make\n\n";

# Make .f95 a valid suffix
#
print MAKEFILE ".SUFFIXES: \n\n";
print MAKEFILE ".SUFFIXES: .f90 .f95 .o .f .c\n\n";
#
# .c -> .o
#
print MAKEFILE ".c.o:\n";
print MAKEFILE "\t\$(CHAMP_C) \$(CHAMP_C_FLAGS) -o \$@ -c \$<\n\n";
#
# .f95 -> .o
#
print MAKEFILE ".f95.o:\n";
print MAKEFILE "\t\$(CHAMP_F95) \$(CHAMP_F95_FLAGS) -o \$@ -c \$<\n\n";
#
# .f90 -> .o
#
print MAKEFILE ".f90.o:\n";
print MAKEFILE "\t\$(CHAMP_F95) \$(CHAMP_F95_FLAGS) -o \$@ -c \$<\n\n";

#
# .f   -> .o
#
print MAKEFILE ".f.o:\n";
print MAKEFILE "\t\$(CHAMP_F77) \$(CHAMP_F77_FLAGS) -o \$@ -c \$<\n\n";
#
# debug:     
#
print MAKEFILE "debug:\n";
print MAKEFILE "\tmake \"CHAMP_EXE=\$\{CHAMP_DEBUG_EXE\}\" \"CHAMP_F95_FLAGS=\$\{CHAMP_F95_DEBUG\}\" \"CHAMP_F77_FLAGS=\$\{CHAMP_F77_DEBUG\}\" \"CHAMP_C_FLAGS=\$\{CHAMP_C_DEBUG\}\" \"CHAMP_LD_FLAGS=\$\{CHAMP_LD_DEBUG\}\"\n\n";
#
# prof:     
#
print MAKEFILE "prof:\n";
print MAKEFILE "\tmake \"CHAMP_F95_FLAGS=\$\{CHAMP_F95_PROF\}\" \"CHAMP_F77_FLAGS=\$\{CHAMP_F77_PROF\}\" \"CHAMP_C_FLAGS=\$\{CHAMP_C_PROF\}\"  \"CHAMP_LD_FLAGS=\$\{CHAMP_LD_PROF\}\"\n\n";
#
# mpi:     
#
print MAKEFILE "mpi:\n";
print MAKEFILE "\tmake \"CHAMP_EXE=\$\{CHAMP_MPI_EXE\}\" \"CHAMP_F95=\$\{CHAMP_MPIF95\}\" \"CHAMP_F77=\$\{CHAMP_MPIF77\}\" \"CHAMP_C=\$\{CHAMP_MPIC\}\" \"CHAMP_F77_FLAGS=\$\{CHAMP_F77_MPI_FLAGS\}\" \"CHAMP_F95_FLAGS=\$\{CHAMP_F95_MPI_FLAGS\}\" \"CHAMP_C_FLAGS=\$\{CHAMP_C_MPI_FLAGS\}\" \"CHAMP_LD=\$\{CHAMP_LD_MPI\}\" \"CHAMP_LD_END=\$\{CHAMP_LD_END_MPI\}\"\n\n";
#
#
# mpi_debug:     
#
print MAKEFILE "mpi_debug:\n";
print MAKEFILE "\tmake \"CHAMP_EXE=\$\{CHAMP_MPI_EXE\}\" \"CHAMP_F95=\$\{CHAMP_MPIF95\}\" \"CHAMP_F77=\$\{CHAMP_MPIF77\}\" \"CHAMP_F77_FLAGS=\$\{CHAMP_F77_MPI_FLAGSDEBUG\}\" \"CHAMP_F95_FLAGS=\$\{CHAMP_F95_MPI_FLAGSDEBUG\}\" \"CHAMP_LD=\$\{CHAMP_LD_MPI\}\" \"CHAMP_LD_END=\$\{CHAMP_LD_END_MPI\}\"\n\n";
#
#
#
# Dependency listings
#
&MakeDependsf95($ARGV[1]);
&MakeDepends("*.f *.F", '^\s*include\s+["\']([^"\']+)["\']');
&MakeDepends("*.c",     '^\s*#\s*include\s+["\']([^"\']+)["\']');
#
# &PrintWords(current output column, extra tab?, word list); --- print words
#    nicely
#
sub PrintWords {
   local($columns) = 78 - shift(@_);
   local($extratab) = shift(@_);
   local($wordlength);
   #
   print MAKEFILE @_[0];
   $columns -= length(shift(@_));
   foreach $word (@_) {
      $wordlength = length($word);
      if ($wordlength + 1 < $columns) {
         print MAKEFILE " $word";
         $columns -= $wordlength + 1;
         }
      else {
         #
         # Continue onto a new line
         #
         if ($extratab) {
            print MAKEFILE " \\\n\t\t$word";
            $columns = 62 - $wordlength;
            }
         else {
            print MAKEFILE " \\\n\t$word";
            $columns = 70 - $wordlength;
            }
         }
      }
   }

#
# &LanguageCompiler(compiler, sources); --- determine the correct language
#    compiler
#
sub LanguageCompiler {
   local($compiler) = &toLower(shift(@_));
   local(@srcs) = @_;
   #
   if (length($compiler) > 0) {
      CASE: {
         grep(/^$compiler$/, ("fc", "f77")) &&
            do { $compiler = "F77"; last CASE; };
         grep(/^$compiler$/, ("cc", "c"))   &&
            do { $compiler = "CHAMP_C"; last CASE; };
         $compiler = "CHAMP_F95";
         }
      }
   else {
      CASE: {
         grep(/\.(f90|f95)$/, @srcs)   && do { $compiler = "CHAMP_F95"; last CASE; };
         grep(/\.(f|F)$/, @srcs) && do { $compiler = "F77";  last CASE; };
         grep(/\.c$/, @srcs)     && do { $compiler = "CHAMP_C";  last CASE; };
         $compiler = "???";
         }
      }
   $compiler;
   }

#
# &toLower(string); --- convert string into lower case
#
sub toLower {
   local($string) = @_[0];
   $string =~ tr/A-Z/a-z/;
   $string;
   }

#
# &uniq(sorted word list); --- remove adjacent duplicate words
#
sub uniq {
   local(@words);
   foreach $word (@_) {
      if ($word ne $words[$#words]) {
         push(@words, $word);
         }
      }
   @words;
   }

#
# &MakeDepends(language pattern, include file sed pattern); --- dependency
#    maker
#
sub MakeDepends {
   local(@incs);
   local($lang) = @_[0];
   local($pattern) = @_[1];
   #
   foreach $file (<${lang}>) {
      open(FILE, $file) || warn "Cannot open $file: $!\n";
      while (<FILE>) {
         /$pattern/i && push(@incs, $1);
         }
      if (defined @incs) {
         $file =~ s/\.[^.]+$/.o/;
         print MAKEFILE "$file: ";
         &PrintWords(length($file) + 2, 0, @incs);
         print MAKEFILE "\n";
         undef @incs;
         }
      }
   }

#
# &MakeDependsf95(f95 compiler); --- FORTRAN 90 dependency maker
#
sub MakeDependsf95 {
   local($compiler) = &toLower(@_[0]);
   local(@dependencies);
   local(%filename);
   local(@incs);
   local(@modules);
   local($objfile);
   #
   # Associate each module with the name of the file that contains it
   #
   foreach $file (<*.f90 *.f95 commons/*.f90>) {
      open(FILE, $file) || warn "Cannot open $file: $!\n";
      while (<FILE>) {
         /^\s*module\s+([^\s!]+)/i &&
            ($filename{&toLower($1)} = $file) =~ s/\.(f90|f95)$/.o/;
         }
      }
   #
   # Print the dependencies of each file that has one or more include's or
   # references one or more modules
   #
   foreach $file (<*.f90 *.f95 commons/*.f90>) {
      open(FILE, $file);
      while (<FILE>) {
         /^\s*include\s+["\']([^"\']+)["\']/i && push(@incs, $1);
         /^\s*include \'mpif\.h\'/i && pop(@incs);  # JT: remove mpif.h from the list of dependancies
         /^\s*use\s+([^\s,!]+)/i && push(@modules, &toLower($1));
         }
      if (defined @incs || defined @modules) {
         ($objfile = $file) =~ s/\.(f90|f95)$/.o/;
         print MAKEFILE "$objfile: ";
         undef @dependencies;
         foreach $module (@modules) {
            push(@dependencies, $filename{$module});
            }
         @dependencies = &uniq(sort(@dependencies));
         &PrintWords(length($objfile) + 2, 0,
                     @dependencies, &uniq(sort(@incs)));
         print MAKEFILE "\n";
         undef @incs;
         undef @modules;
         #
         # Cray CHAMP_F95 compiler
         #
         if ($compiler eq "cray") {
            print MAKEFILE "\t\$(CHAMP_F95) \$(CHAMP_F95_FLAGS) -c ";
            foreach $depend (@dependencies) {
               push(@modules, "-p", $depend);
               }
            push(@modules, $file);
            &PrintWords(30, 1, @modules);
            print MAKEFILE "\n";
            undef @modules;
            }
         #
         # ParaSoft CHAMP_F95 compiler
         #
         if ($compiler eq "parasoft") {
            print MAKEFILE "\t\$(CHAMP_F95) \$(CHAMP_F95_FLAGS) -c ";
            foreach $depend (@dependencies) {
               $depend =~ s/\.o$/.(f90|f95)/;
               push(@modules, "-module", $depend);
               }
            push(@modules, $file);
            &PrintWords(30, 1, @modules);
            print MAKEFILE "\n";
            undef @modules;
            }
         }
      }
   }
