# New Makefile that automatically depends itself
#
# $Id: Makefile,v 1.3 1996/12/17 19:52:37 chase Exp $
#

IFLAGS = -I. -IGCC_4.5
DFLAGS = 
CXX = g++
CC  = $(GCC)
GCC = gcc
LD  = $(CXX)

LIBS = 

WFLAGS = -Wall 
SYMFLAGS = -g
SPECIALFLAGS = -std=c++0x

PROFILE = #-pg 
OPTFLAGS =#-O
CFLAGS = $(OPTFLAGS) $(SPECIALFLAGS) $(PROFILE) $(WFLAGS) $(IFLAGS) $(SYMFLAGS)
CXXFLAGS = $(CFLAGS)
CPPFLAGS = $(IFLAGS) $(DFLAGS)
LDFLAGS = $(PROFILE) -g 

PROGRAM = vtest
CXXSRCS = vtest.cpp
          
YACCSRCS = 
LEXSRCS = 
CSRCS =
F77SRCS =
ASMSRCS =

SRCS = $(CXXSRCS) $(CSRCS) $(F77SRCS) $(ASMSRCS) 

OBJS = $(CXXSRCS:.cpp=.o) $(CSRCS:.c=.o) $(F77SRCS:.F=.o) $(ASMSRCS:.S=.o) \
       $(YACCSRCS:.y=.o) $(LEXSRCS:.l=.o)

all: $(PROGRAM)


$(PROGRAM): $(OBJS)
	$(LD) -o $@ $(LDFLAGS) $(OBJS) $(LIBS)

clean:
	rm -f $(OBJS) $(PROGRAM)

tidy:
	rm -f *.BAK *.bak *.CKP

undepend:
	rm -f $(OBJS:%.o=.%.d) 

spotless: tidy clean undepend

.y.cpp:
	$(BISON) $(BISONFLAGS) -o $@ $<
	mv $@.h $*.h
	mv $@.output $*.output
.l.cpp:
	$(FLEX) ${FLEXFLAGS} -t $< > $@

# auto depend stuff for GNU make only
depend: undepend
	@echo ""
	@echo "Dependences are handled automatically, just \"make\""

ifneq ($(strip $(CSRCS)),)
.%.d: %.c 
	$(SHELL) -ec '$(GCC) -MM $(CPPFLAGS) $< > $@'


include $(CSRCS:%.c=.%.d)
endif 

ifneq ($(strip $(CXXSRCS)),)
.%.d: %.cpp
	$(SHELL) -ec '$(GCC) -MM $(CPPFLAGS) $< > $@'

include $(CXXSRCS:%.cpp=.%.d) 
endif 

# Not known to work for F77 (.F) or ASM (.S)
#
# ifneq ($(F77SRCS),)
# include $(F77SRCS:%.F=.%.d) 
# endif # F77SRCS

# ifneq ($(ASMSRCS),)
# include $(ASMSRCS:%.S=.%.d)
# endif # ASMSRCS
