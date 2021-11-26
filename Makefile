TARGETS = interp_lib.so
OBJS=$(TARGETS:=.o)

all: $(TARGETS)

interp_lib.o: interp_lib.c
	gcc -I. -Wall -fPIC $(shell python3-config --cflags) -c -o $@ $<

interp_lib.so: interp_lib.o
	gcc -Wall -fPIC -I. $(shell python3-config --cflags) $< -o $@ -shared $(shell python3-config --cflags)

$(OBJS): Makefile

test: ${TARGETS}
	python3 interp_test.py