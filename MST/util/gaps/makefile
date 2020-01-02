VPATH = $(ODIR)/../util/gaps

OBJECTS=$(ODIR)/fill.o $(ODIR)/gapslib.o $(ODIR)/radii.o $(ODIR)/symbols.o

all: $(ODIR)/gaps

$(ODIR)/gaps: gaps.c $(OBJECTS) gaps.h 
	$(MPICC) $(CFLAGS) -lm gaps.c $(OBJECTS) -o $@

$(ODIR)/fill.o: fill.c fill.h
	$(MPICC) $(CFLAGS) -c fill.c

$(ODIR)/gapslib.o: gapslib.c gapslib.h 
	$(MPICC) $(CFLAGS) -c gapslib.c

$(ODIR)/radii.o: radii.cpp
	$(MPICC) $(CFLAGS) -c $<

$(ODIR)/symbols.o: symbols.c
	$(MPICC) $(CFLAGS) -c $<

# Preview the man page
man:
	man ./gaps.1

clean:
	rm -f $(ODIR)/gaps $(OBJECTS)
