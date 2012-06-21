OBJECTS = main.o vmgenome.o vm.o

all: vmgenome

vmgenome: $(OBJECTS)
	gcc -g $(OBJECTS) -o vmgenome

%.o : %.c
	gcc -g -c $< -o $@
