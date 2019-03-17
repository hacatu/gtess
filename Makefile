CC := gcc
LD := gcc
AR := ar
#CFLAGS_RELEASE := -Wall -std=c99 -Iinclude -Ofast -ffast-math -mtune=native -march=native
CFLAGS_RELEASE := -Wall -std=c99 -Iinclude -Og -g3
#LDFLAGS_RELEASE := -Ofast -ffast-math -mtune=native -march=native -flto -fuse-ld=gold -lm -lgsl -lgslcblas -lSDL2 -lSDL2_image -lSDL2_gfx
LDFLAGS_RELEASE := -Og -g3 -lm -lgsl -lgslcblas -lSDL2 -lSDL2_image -lSDL2_gfx
CFLAGS_DEBUG := -Wall -std=c99 -Iinclude -g3 --coverage
LDFLAGS_DEBUG := --coverage -lm -lgsl  -lgslcblas -lSDL2 -lSDL2_image -lSDL2_gfx
TEST_CFLAGS := -Wall -std=c99 -Iinclude -g3
TEST_LDFLAGS := --coverage -lm -lgsl -lgslcblas -lSDL2 -lSDL2_image -lSDL2_gfx
BUILD_ROOT := $(shell pwd)
SOURCES := $(shell find src -maxdepth 1 -name '*.c')
TEST_SOURCES := $(shell find src/test -maxdepth 1 -name '*.c')
BINARY_SOURCES := $(shell find src/binary -maxdepth 1 -name '*.c')
OBJS_DEBUG := $(patsubst src/%.c,obj/debug/%.o,$(SOURCES))
OBJS_RELEASE := $(patsubst src/%.c,obj/release/%.o,$(SOURCES))
TEST_BINARIES := $(patsubst src/test/%.c,bin/test/%,$(TEST_SOURCES))
BINARIES := $(patsubst src/binary/%.c,bin/%,$(BINARY_SOURCES))

all: coverage binaries | Makefile

.SECONDARY:
obj/bin obj/test obj/debug obj/release bin bin/test notes cov:
	mkdir -p $@

bin/%: obj/bin/%.o $(OBJS_RELEASE) | bin
	$(CC) $(LDFLAGS_RELEASE) -o $@ $^

bin/test/%: obj/test/%.o $(OBJS_DEBUG) | bin/test 
	$(CC) $(TEST_LDFLAGS) -o $@ $^

obj/release/%.o: src/%.c | obj/release
	$(CC) $(CFLAGS_RELEASE) -c -o $@ $^

obj/debug/%.o: src/%.c | obj/debug notes
	$(CC) $(CFLAGS_DEBUG) -c $^ -o notes/$*.o
	mv notes/$*.o $@

obj/test/%.o: src/test/%.c | obj/test
	$(CC) $(TEST_CFLAGS) -c $^ -o $@

obj/bin/%.o: src/binary/%.c | obj/bin
	$(CC) $(CFLAGS_RELEASE) -c -o $@ $^

.PHONY:
binaries: $(BINARIES)

.PHONY:
test: $(TEST_BINARIES)
	for test in $^ ; do ./$$test; done

.PHONY:
coverage: test cov
	gcovr --html-details -e src/test/.\*\.c -o cov/index.html notes

clean:
	rm -rf bin obj notes cov

