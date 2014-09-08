CC = clang
CC_FLAGS = -Wall -O0 -fPIC -g -std=c99 -pipe -I/usr/share/R/include 
#CC_FLAGS = -w -std=c11 -fPIC -S -save-temps -fverbose-asm -masm=intel -I/usr/share/R/include


LIB = src/cec.so
EXEC = src/cec_test
SOURCES = $(wildcard src/*.c)
OBJECTS = $(SOURCES:.c=.o)

# Main target
$(LIB): $(OBJECTS)
	$(CC) -shared -O0 -fPIC -g $(OBJECTS) -o $(LIB) -L/usr/lib/R/lib -lR -lm -llapack

$(EXEC): $(OBJECTS)
	$(CC) -O0 -g $(OBJECTS) -o $(EXEC) -L/usr/lib/R/lib -lR -lm -llapack



%.o: %.c
	$(CC) -c $(CC_FLAGS) $< -o $@


clean:
	rm -f $(EXEC) $(LIB) $(OBJECTS) *.s *.so src/*.s src/*.so