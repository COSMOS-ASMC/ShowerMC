# Copyright (c) 1991 Peter J. Nicklin.
# Copyright (c) 1991 Version Technology.
# All Rights Reserved.
#
# $License: VT.1.1 $
# Redistribution and use in source and binary forms,  with or without
# modification,  are permitted provided that the following conditions
# are met:  (1) Redistributions of source code must retain the  above
# copyright  notice,  this  list  of  conditions  and  the  following
# disclaimer.  (2) Redistributions in binary form must reproduce  the
# above  copyright notice,  this list of conditions and the following
# disclaimer in the  documentation  and/or other  materials  provided
# with  the  distribution.  (3) All advertising materials  mentioning
# features or  use  of  this  software  must  display  the  following
# acknowledgement:  ``This  product  includes  software  developed by
# Version Technology.''  Neither the name of Version  Technology  nor
# the  name  of  Peter J. Nicklin  may  be used to endorse or promote
# products derived from this software without specific prior  written
# permission.
#
# THIS SOFTWARE IS PROVIDED BY VERSION TECHNOLOGY ``AS IS''  AND  ANY
# EXPRESS OR IMPLIED WARRANTIES,  INCLUDING,  BUT NOT LIMITED TO, THE
# IMPLIED  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL  VERSION  TECHNOLOGY  BE
# LIABLE  FOR ANY DIRECT,  INDIRECT,  INCIDENTAL, SPECIAL, EXEMPLARY,
# OR  CONSEQUENTIAL DAMAGES   (INCLUDING,   BUT   NOT   LIMITED   TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;  LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
# OF  LIABILITY,  WHETHER  IN  CONTRACT,  STRICT LIABILITY,  OR  TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING  IN ANY WAY OUT OF THE
# USE OF THIS SOFTWARE,  EVEN  IF  ADVISED OF THE POSSIBILITY OF SUCH
# DAMAGE.
#
# Report problems and direct questions to nicklin@netcom.com
#
# $Header: Makefile,v 4.2 91/11/26 07:24:37 nicklin Exp $
#
# Author: Peter J. Nicklin
#
#INSTALLDIR    = /usr/contrib
INSTALLDIR    = /usr/local

MAKEFILE      = Makefile

SHELL	      = /bin/sh

SUBDIR	      = src templates man/man1

all:;		@for i in $(SUBDIR);\
		do (\
			echo Making $$i ...;\
			cd $$i;\
			$(MAKE) DESTDIR=$(DESTDIR)\
			 	INSTALLDIR=$(INSTALLDIR)\
		); done

clean:;		@for i in $(SUBDIR);\
		do (\
			echo Cleaning $$i ...;\
			cd $$i;\
			$(MAKE) clean\
		); done

clobber:;	@for i in $(SUBDIR);\
		do (\
			echo Clobbering $$i ...;\
			cd $$i;\
			$(MAKE) clobber\
		); done

depend:;	@for i in $(SUBDIR);\
		do (\
			echo Creating dependencies for $$i ...;\
			cd $$i;\
			$(MAKE) DESTDIR=$(DESTDIR)\
				INSTALLDIR=$(INSTALLDIR) depend\
		); done

install:;	@for i in $(SUBDIR);\
		do (\
			echo Installing $$i ...;\
			cd $$i;\
			$(MAKE) DESTDIR=$(DESTDIR)\
				INSTALLDIR=$(INSTALLDIR) install\
		); done
