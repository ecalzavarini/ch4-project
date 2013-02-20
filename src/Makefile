CC=mpicc 
CFLAGS= -Wall  
LFLAGS= -lm 
SOURCES= parallel.c grid.c parameters.c lb.c initial_conditions.c #fluid.c
OBJECTS=$(SOURCES:.c=.o)
HEADERS= common_object.h common_main.h typedef.h functions.h define.h
EXECUTABLE= ch4

%.o: %.c $(HEADERS) 
	$(CC) -c $(CFLAGS) $< -o $@ 

$(EXECUTABLE):  main.c $(OBJECTS) $(HEADERS)
	$(CC) $(CFLAGS) -o $(EXECUTABLE) main.c $(OBJECTS) $(LFLAGS)

clean:
	rm -rf ./*.o ./*~ $(EXECUTABLE) ./*.dat ./RUN

tar:
	tar cvzf ch4-project_backup.tgz ./*.c ./*.h ./*.sh Makefile 