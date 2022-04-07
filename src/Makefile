PKG_NAME=qap-code


CC=gcc
FC=gfortran
AR=ar

#CFLAGS=-Wall -g
CFLAGS=-Wall -fomit-frame-pointer -O3 -W -Wno-unused-parameter

# for profiling

#CFLAGS=-O2 -pg -Wall


#FCFLAGS = -g -fbounds-check
FCFLAGS = -O3

OBJS=main.o qap-utils.o tools.o 
EXECS=eo-qap rots-qap sa-qap fant-qap brute-force qap-new-format check-sol mk-grey


#EXECS+=qapeli  qapglb


all: $(EXECS)

eo-qap: eo-qap.c $(OBJS) eo-pdf.o

%: %.c $(OBJS)
	$(CC) -o $@ $(CFLAGS) $^ -lm

%: %.f
	$(FC) -o $@ $(FCFLAGS) $^



qap-new-format: qap-new-format.c qap-utils.o
	$(CC) -o $@ $(CFLAGS) $^

check-sol: check-sol.c qap-utils.o tools.o
	$(CC) -o $@ $(CFLAGS) $^

mk-grey: mk-grey.c qap-utils.o
	$(CC) -o $@ $(CFLAGS) $^


qap-utils.o: qap-utils.h

main.o: main.h qap-utils.h

tools.o: tools.h

eo-pdf.o: eo-pdf.h tools.h

# utils

find_tau: find_tau.c
	$(CC)  -o $@ $(CFLAGS) $^


# distribution

ROOT_DIR=$(shell pwd)

distdir=/tmp
distpath=$(distdir)/$(PKG_NAME)

#TAR_NAME=$(PKG_NAME).tgz
TAR_NAME=$(PKG_NAME)-$(shell date +'%Y%m%d').tgz

dist: $(TAR_NAME).tgz

$(TAR_NAME).tgz: copy_dist_tree
	@(cd $(distdir); tar cfz $(ROOT_DIR)/$(TAR_NAME) $(PKG_NAME))
	@banner="$(TAR_NAME) is ready for distribution"; \
	dashes=`echo "$$banner" | sed s/./=/g`; \
	echo "$$dashes"; \
	echo "$$banner"; \
	echo "$$dashes"


copy_dist_tree:
	@-rm -rf $(distpath); mkdir -p $(distpath)
	@(tar cf - `cat DISTRIB_FILES` | (cd $(distpath); tar xf -))


# cleaning

clean:
	rm -f *.o *.a *.d *~ .*~ \#*# .#*# $(EXECS)
