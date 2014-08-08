all: simulator

OPTFLAGS = -O3

#LDFLAGS = -static

CFLAGS += $(OPTFLAGS)

#CFLAGS += $(LDFLAGS)


Astroid:
	g++ -std=c++0x $(CFLAGS) -o simulator src/general_functions.cpp src/def_sample.cpp src/def_transcript.cpp src/sim_gap.cpp

clean: