#CC=mpicc 
CC=h5pcc #-DH5_USE_16_API
CFLAGS= -Wall #-g #-pg 
LFLAGS= -lm -lhdf5
SOURCES= parallel.c grid.c parameters.c lb.c boundary_conditions.c initial_conditions.c fluid.c output.c output_h5.c melting.c
OBJECTS=$(SOURCES:.c=.o)
HEADERS= common_object.h common_main.h typedef.h functions.h define.h
EXECUTABLE= ch4

%.o: %.c $(HEADERS) 
	$(CC) -c $(CFLAGS) $< -o $@ 

$(EXECUTABLE):  main.c $(OBJECTS) $(HEADERS)
	$(CC) $(CFLAGS) -o $(EXECUTABLE) main.c $(OBJECTS) $(LFLAGS)

clean:
	rm -rf ./*.o ./*~ $(EXECUTABLE) ./*.dat ./*.out ./RUN

clean-dat:
	rm -rf ./*.dat

tar:
	tar cvzf ch4-project_backup.tgz ./*.c ./*.h Makefile 