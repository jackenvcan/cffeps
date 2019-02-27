LIBS=-lm
INCLUDES=-I../../include
OBJS=../../obj/bdlffmc.o \
	 ../../obj/emissions.o \
	 ../../obj/energy.o \
	 ../../obj/eqffmc.o \
	 ../../obj/hrffmc.o \
	 ../../obj/fbp_2009.o \
	 ../../obj/fbp_2009_pro.o \
	 ../../obj/fwi84.o \
	 ../../obj/thermo.o \
	 ../../obj/julian.o \
	 ../../obj/plume.o \
	 ../../obj/reademissions.o \
	 ../../obj/readfeps.o \
	 ../../obj/readprofile.o \
	 ../../obj/writefeps.o \

CFFEPS.exe: CFFEPS.c $(OBJS)
	gcc CFFEPS.c -o CFFEPS.exe $(OBJS) $(LIBS) $(INCLUDES)

../../obj/bdlffmc.o: ../diurnal/bdlffmc.c 
	gcc -c ../diurnal/bdlffmc.c -o ../../obj/bdlffmc.o $(INCLUDES)

../../obj/CFFEPS.o: CFFEPS.c 
	gcc -c CFFEPS.c -o ../../obj/CFFEPS.o  $(INCLUDES)
    
../../obj/emissions.o: emissions.c 
	gcc -c emissions.c -o ../../obj/emissions.o  $(INCLUDES)

../../obj/energy.o: energy.c 
	gcc -c energy.c -o ../../obj/energy.o  $(INCLUDES)

../../obj/eqffmc.o: ../diurnal/eqffmc.c 
	gcc -c ../diurnal/eqffmc.c -o ../../obj/eqffmc.o $(INCLUDES)
	
../../obj/hrffmc.o: ../diurnal/hrffmc.c 
	gcc -c ../diurnal/hrffmc.c -o ../../obj/hrffmc.o $(INCLUDES)

../../obj/fbp_2009.o: ../fbp2009/fbp_2009.c 
	gcc -c ../fbp2009/fbp_2009.c -o ../../obj/fbp_2009.o  $(INCLUDES)

../../obj/fbp_2009_pro.o: ../fbp2009/fbp_2009_pro.c 
	gcc -c ../fbp2009/fbp_2009_pro.c -o ../../obj/fbp_2009_pro.o  $(INCLUDES)

../../obj/fwi84.o: ../fwi84/fwi84.c 
	gcc -c ../fwi84/fwi84.c -o ../../obj/fwi84.o $(INCLUDES)

../../obj/julian.o: ../utils/julian.c 
	gcc -c ../utils/julian.c -o ../../obj/julian.o $(INCLUDES)

../../obj/plume.o: plume.c 
	gcc -c plume.c -o ../../obj/plume.o  $(INCLUDES)

../../obj/reademissions.o: ../utils/reademissions.c 
	gcc -c ../utils/reademissions.c -o ../../obj/reademissions.o  $(INCLUDES)

../../obj/readfeps.o: ../utils/readfeps.c 
	gcc -c ../utils/readfeps.c -o ../../obj/readfeps.o  $(INCLUDES)

../../obj/readprofile.o: ../utils/readprofile.c 
	gcc -c ../utils/readprofile.c -o ../../obj/readprofile.o  $(INCLUDES)

../../obj/thermo.o: ../thermo/thermo.c 
	gcc -c ../thermo/thermo.c -o ../../obj/thermo.o $(INCLUDES)

../../obj/writefeps.o: ../utils/writefeps.c 
	gcc -c ../utils/writefeps.c -o ../../obj/writefeps.o  $(INCLUDES)


				
