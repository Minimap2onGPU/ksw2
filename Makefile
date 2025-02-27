CC=			gcc
CXX=		g++
CFLAGS=		-g -Wall -Wextra -Wc++-compat -O2
CXXFLAGS=	-g -Wall -Wextra
CPPFLAGS=	-g -Wall -Wextra -O2 #-DHAVE_KALLOC 
INCLUDES=	-I.
OBJS=		ksw2_gg.o ksw2_gg2.o ksw2_gg2_sse.o ksw2_extz.o ksw2_extz2_sse.o \
			ksw2_extd.o ksw2_extd2_sse.o ksw2_extf2_sse.o ksw2_exts2_sse.o \
			ksw2_extd2.o ksw2_extd2_cpp.o
PROG=		ksw2-test
LIBS=		-lz
coverage = n

ifeq ($(coverage),y)
	CFLAGS += -fprofile-arcs -ftest-coverage
	CXXFLAGS += -fprofile-arcs -ftest-coverage
else
	CXXFLAGS += -O2
endif

ifneq ($(gaba),) # gaba source code directory
	CPPFLAGS += -DHAVE_GABA
	INCLUDES += -I$(gaba)
	LIBS_MORE += -L$(gaba)/build -lgaba
	CFLAGS += -msse4
endif

ifneq ($(parasail),) # parasail install prefix
	CPPFLAGS += -DHAVE_PARASAIL
	INCLUDES += -I$(parasail)/include
	LIBS_MORE += $(parasail)/lib/libparasail.a # don't link against the dynamic library
endif

ifeq ($(sse2),)
	CFLAGS += -march=native
endif

ifneq ($(avx2),)
	CFLAGS += -mavx2
endif


.SUFFIXES:.c .o _cpp.o .cpp

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

debug: CPPFLAGS += -DDEBUG
debug: $(PROG)

.cpp_cpp.o:
		$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

#ksw2-test:cli.o kalloc.o $(OBJS)
#		$(CC) $(CFLAGS) $^ -o $@ $(LIBS_MORE) $(LIBS)

report:
	mkdir -p coverage
	lcov --capture --directory . --output-file coverage.info
	genhtml coverage.info --output-directory coverage

ksw2-test:cli.o kalloc.o $(OBJS)
		$(CXX) $(CFLAGS) $^ -o $@ $(LIBS_MORE) $(LIBS)
		
clean:
		rm -fr gmon.out *.o a.out $(PROG) $(PROG_EXTRA) *~ *.a *.dSYM session*
		rm -rf *.gcda *.gcno

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c)

# DO NOT DELETE

cli.o: ksw2.h kseq.h
kalloc.o: kalloc.h
ksw2_extd.o: ksw2.h
ksw2_extd2_sse.o: ksw2.h
ksw2_extf2_sse.o: ksw2.h
ksw2_extz.o: ksw2.h
ksw2_extz2_sse.o: ksw2.h
ksw2_gg.o: ksw2.h
ksw2_gg2.o: ksw2.h
ksw2_gg2_sse.o: ksw2.h
ksw2_extd2.o: ksw2.h
ksw2_extd2_cpp.o: ksw2.h
