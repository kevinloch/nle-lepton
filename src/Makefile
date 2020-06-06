CC = gcc
CFLAGS = -I. -Wall -Ofast
LIBS = -lm
OBJ = cscanner.o getFormulaStr.o initDynamicfactorArray.o initInfactorArray.o initOutfactorArray.o initSmrfactorArray.o phase1.o phase2.o util.o selectOutputs.o verifyMatches.o usage.o nle-config.o generateExponents.o nle-lepton.o
DEPS = cscanner.h getFormulaStr.h initDynamicfactorArray.h initInfactorArray.h initOutfactorArray.h initSmrfactorArray.h phase1.h phase2.h util.h selectOutputs.h verifyMatches.h usage.h nle-config.h generateExponents.h nle-lepton.h

.PHONY: all clean

all: nle-lepton

clean:
	rm -f nle-lepton $(OBJ)

$(OBJ): %.o : %.c $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $<

nle-lepton: $(OBJ)
	$(CC) $(CFLAGS) $(LIBS) -o nle-lepton $^
