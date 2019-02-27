LIBS=-lm
INCLUDES=-Iinclude
# DEBUG=-g

# For some reason, using the default -O2 optimization level generate errors in the CSV decoding.
# TODO: Fix that with appropriate strncompare in cffeps.c
CFLAGS=-O0

BIN=bin/cffeps
OBJS=obj/bdlffmc.o \
	 obj/emissions.o \
	 obj/energy.o \
	 obj/eqffmc.o \
	 obj/hrffmc.o \
	 obj/fbp_2009.o \
	 obj/fbp_2009_pro.o \
	 obj/fwi84.o \
	 obj/thermo.o \
	 obj/julian.o \
	 obj/plume.o \
	 obj/reademissions.o \
	 obj/readfeps.o \
	 obj/readprofile.o \
	 obj/writefeps.o 

CC=icc

.PHONY: clean distclean

cffeps: src/CFFEPS/cffeps.c $(OBJS)
	$(CC) $(CFLAGS) $(DEBUG) src/CFFEPS/cffeps.c -o $(BIN) $(OBJS) $(LIBS) $(INCLUDES)

obj/bdlffmc.o: src/diurnal/bdlffmc.c 
	$(CC) -c $(CFLAGS) $(DEBUG) src/diurnal/bdlffmc.c -o obj/bdlffmc.o $(INCLUDES)

obj/cffeps.o: src/CFFEPS/cffeps.c 
	$(CC) -c $(CFLAGS) $(DEBUG) cffeps.c -o obj/cffeps.o  $(INCLUDES)

obj/emissions.o: src/CFFEPS/emissions.c 
	$(CC) -c $(CFLAGS) $(DEBUG) src/CFFEPS/emissions.c -o obj/emissions.o  $(INCLUDES)

obj/energy.o: src/CFFEPS/energy.c 
	$(CC) -c $(CFLAGS) $(DEBUG) src/CFFEPS/energy.c -o obj/energy.o  $(INCLUDES)

obj/eqffmc.o: src/diurnal/eqffmc.c 
	$(CC) -c $(CFLAGS) $(DEBUG) src/diurnal/eqffmc.c -o obj/eqffmc.o $(INCLUDES)
	
obj/hrffmc.o: src/diurnal/hrffmc.c 
	$(CC) -c $(CFLAGS) $(DEBUG) src/diurnal/hrffmc.c -o obj/hrffmc.o $(INCLUDES)

obj/fbp_2009.o: src/fbp2009/fbp_2009.c 
	$(CC) -c $(CFLAGS) $(DEBUG) src/fbp2009/fbp_2009.c -o obj/fbp_2009.o  $(INCLUDES)

obj/fbp_2009_pro.o: src/fbp2009/fbp_2009_pro.c 
	$(CC) -c $(CFLAGS) $(DEBUG) src/fbp2009/fbp_2009_pro.c -o obj/fbp_2009_pro.o  $(INCLUDES)

obj/fwi84.o: src/fwi84/fwi84.c 
	$(CC) -c $(CFLAGS) $(DEBUG) src/fwi84/fwi84.c -o obj/fwi84.o $(INCLUDES)

obj/julian.o: src/utils/julian.c 
	$(CC) -c $(CFLAGS) $(DEBUG) src/utils/julian.c -o obj/julian.o $(INCLUDES)

obj/plume.o: src/CFFEPS/plume.c 
	$(CC) -c $(CFLAGS) $(DEBUG) src/CFFEPS/plume.c -o obj/plume.o  $(INCLUDES)

obj/reademissions.o: src/utils/reademissions.c 
	$(CC) -c $(CFLAGS) $(DEBUG) src/utils/reademissions.c -o obj/reademissions.o  $(INCLUDES)

obj/readfeps.o: src/utils/readfeps.c 
	$(CC) -c $(CFLAGS) $(DEBUG) src/utils/readfeps.c -o obj/readfeps.o  $(INCLUDES)

obj/readprofile.o: src/utils/readprofile.c 
	$(CC) -c $(CFLAGS) $(DEBUG) src/utils/readprofile.c -o obj/readprofile.o  $(INCLUDES)

obj/thermo.o: src/thermo/thermo.c 
	$(CC) -c $(CFLAGS) $(DEBUG) src/thermo/thermo.c -o obj/thermo.o $(INCLUDES)

obj/writefeps.o: src/utils/writefeps.c 
	$(CC) -c $(CFLAGS) $(DEBUG) src/utils/writefeps.c -o obj/writefeps.o  $(INCLUDES)

clean:
	rm -f $(OBJS)

distclean: clean
	rm -f $(BIN)
