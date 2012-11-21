OBJECTS = main.o vmgenome.o vm.o levenshtein.o

all: vmgenome

vmgenome: $(OBJECTS)
	gcc -g $(OBJECTS) -o vmgenome

%.o : %.c
	gcc -g -std=gnu99 -c $< -o $@
