cat := /bin/cat
all-tests := $(addsuffix .test, $(basename $(wildcard test/*.in)))

.PHONY : clean test %.test

all : segA

segA : segA.c lnfac.h genwin.h genwin.o
	cc -O -o segA segA.c genwin.o -lm

genwin.o : genwin.c genwin.h alphabet.h
	cc -O -c genwin.c
 
alphabet.h : alphabet.txt
	xxd -i alphabet.txt > alphabet.h



clean:
	rm -f segA segA.o genwin.o alphabet.h

test : $(all-tests)
	@ >&2 echo "SUCCESS" 

test/%.test : test/%.in test/%.cmp test/seq.fna segA
	$(cat) $< | xargs segA test/seq.fna 2>&1 | diff -q $(word 2, $?) - > /dev/null || (echo Test $@ failed && exit 1)

