CFLAGS = -funroll-all-loops -O2
CC = g++
INC = -I$(HOME)/include -I/usr/local/dmp/nifti/include
LIBS = -L$(HOME)/lib -L/usr/local/dmp/nifti/lib -lbabak_lib_linux -lniftiio -lznz -lpthread -lm -lz -lc

all: brainwash brainwash_merge avgcurve jaccardindex diceindex booster

diceindex: diceindex.cxx

	$(CC) $(CFLAGS) -o diceindex diceindex.cxx $(LIBS) $(INC)

jaccardindex: jaccardindex.cxx
	$(CC) $(CFLAGS) -o jaccardindex jaccardindex.cxx $(LIBS) $(INC)

brainwash: brainwash.cxx
	$(CC) $(CFLAGS) -o brainwash brainwash.cxx $(LIBS) $(INC)

brainwash_merge: brainwash_merge.cxx
	$(CC) $(CFLAGS) -o brainwash_merge brainwash_merge.cxx $(LIBS) $(INC)

booster: booster.cxx
	$(CC) $(CFLAGS) -o booster booster.cxx $(LIBS) $(INC)

avgcurve: avgcurve.cxx
	$(CC) $(CFLAGS) -o avgcurve avgcurve.cxx $(LIBS) $(INC)
