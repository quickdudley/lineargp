OBJECTS = vmgenome.o vm.o levenshtein.o pareto.o -lm

all: vmgenome runvm

clean:
	rm vmgenome $(OBJECTS)

vmgenome: main.o $(OBJECTS)
	gcc -g main.o $(OBJECTS) -o vmgenome

runvm: runvm.o $(OBJECTS)
	gcc -g runvm.o $(OBJECTS) -o runvm

%.o : %.c
	gcc -g -std=gnu99 -c $< -o $@
