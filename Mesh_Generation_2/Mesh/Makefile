LIBS=Flotation2 USFs USF_GroundedAreaCoeff

all: $(LIBS)

%: %.F90
	elmerf90 $^ -o $@


clean:
	rm -f $(LIBS)
