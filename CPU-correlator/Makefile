CXX		?= g++
CFLAGS	        ?= -O3 -msse3 -Wall -march=nocona -m64 -mfpmath=sse -g
LDFLAGS         ?= -m64 -lpthread

CC_SOURCES	= main.cc timer.cc cpu_correlator_1x1.cc cpu_correlator_1x1_sse3.cc cpu_correlator_1x1_sse3.cc cpu_correlator_1x1_time_sse3.cc reference_correlator.cc cpu_correlator_2x2_sse3.cc cpu_correlator_2x2_time_sse3.cc cpu_correlator_3x2_time_sse3.cc

CC_OBJECTS	= $(CC_SOURCES:%.cc=%.o)

CC_AS		= $(CC_SOURCES:%.cc=%.s)

%.o:		%.cc
		$(CXX) $(CFLAGS) -c $< -o $@

%.s:		%.cc
		$(CXX) -S $(CFLAGS) -c $< -o $@

correlator:	$(CC_OBJECTS)
		$(CXX) $(LDFLAGS) $^ -o $@

s:		$(CC_AS)

clean:
		rm -f *.o correlator *~
