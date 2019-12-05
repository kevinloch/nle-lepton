CC = gcc
CFLAGS = -Ofast -I.
LIBS = -lm
OBJ = cscanner.o reference.o initMultiplierArray.o phase1.o phase2.o util.o verifyMatches.o usage.o nle-lepton.o
DEPS = cscanner.h reference.h initMultiplierArray.h phase1.h phase2.h util.h verifyMatches.h usage.h nle-lepton.h

.PHONY: all clean

all: nle-lepton

clean:
	rm -f nle-lepton $(OBJ)

$(OBJ): %.o : %.c $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $<

nle-lepton: $(OBJ)
	$(CC) $(CFLAGS) $(LIBS) -o nle-lepton $^
