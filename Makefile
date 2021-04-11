#!/usr/bin/env make -f
# Makefile to compile SPM MB GMM lib C-MEX files (hacked together by MB)
#
# Copyright (C) 1991-2021 Wellcome Centre for Human Neuroimaging
#
# $Id: Makefile 8058 2021-02-10 10:38:31Z guillaume $
#
###############################################################################
#
# This Makefile has been tested under Linux, Windows and MacOS.
#
# If you have to tweak this Makefile or Makefile.var to compile the SPM
# mex-files for your platform, please send the details to <fil.spm@ucl.ac.uk>
# so they can be included here.
#
# To compile and install SPM, type the following in a Unix console:
# >  make distclean && make && make install
#
# You can specify a particular platform with the following syntax:
# >  make PLATFORM=Your_Platform
# The standard targets are 'all', 'clean', 'distclean', 'doc' and 'install'.
#
# For a list of compatible compilers, see
#    https://www.mathworks.com/support/compilers.html
#
###############################################################################

include Makefile.var

###############################################################################
# Objects to go in the archive and mexfiles
###############################################################################

SPMMEX  =\
	spm_gmmlib.$(MEXEXT)

###############################################################################
# Public make targets
###############################################################################

all: verb.$(MEXEXT) main-all verb.all.end

clean: verb.distclean main-distclean verb.distclean.end

###############################################################################
# Private make targets
###############################################################################

main-all: $(SPMMEX)

main-distclean:
	$(DEL) $(SPMMEX)

###############################################################################
# Compile the mex files themselves
###############################################################################

spm_gmmlib.$(MEXEXT): spm_gmmlib.c gmmlib.c gmmlib.h 
	$(MEX) spm_gmmlib.c gmmlib.c $(MEXEND)

###############################################################################
# Display Messages
###############################################################################

verb.clean:
	$(call verb, "Deleting object (.o) files")

verb.distclean:
	$(call verb, "Deleting MEX (.$(MEXEXT)) and archive (.a) files")

verb.install:
	$(call verb, "Installing MEX files")

verb.tarball:
	$(call verb, "Creating archive spm_mex.tar.gz")

verb.mexw32:
	$(call verb, "Windows compilation (32 bit)")

verb.mexw64:
	$(call verb, "Windows compilation (64 bit)")

verb.mexglx:
	$(call verb, "Linux compilation (x86-32)")

verb.mexa64:
	$(call verb, "Linux compilation (x86-64)")

verb.mexmaci:
	$(call verb, "MacOS compilation (Intel 32 bit)")

verb.mexmaci64:
	$(call verb, "MacOS compilation (Intel 64 bit)")

verb.mex:
	$(call verb, "${PLATFORM} compilation (`${MEXBIN} -v | head -n 1`)")

verb.all.end:
	$(call verb, "Compilation: done")

verb.distclean.end:
	$(call verb, "Distclean: done")

verb.install.end:
	$(call verb, "Installation: done")

verb.external:
	$(call verb, "In external")
