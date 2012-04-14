OBJECTS = main.o vmgenome.o vm.o

all: vmgenome

vmgenome: $(OBJECTS)
	gcc $(OBJECTS) -o vmgenome

%.o : %.cc
	gcc -c $< -o $@
