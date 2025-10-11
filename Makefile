CC	:= g++
INCDIRS	:=
LIBDIRS	:=
LIBS	:= -lcsdr++ -lfftw3f
CFLAGS	:= -O1 $(INCDIRS)
OBJECTS	:= skimmer.o bufmodule.o

all: csdr-rttyskimmer

csdr-rttyskimmer: $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBDIRS) $(LIBS)

clean:
	rm -f $(OBJECTS) csdr-rttyskimmer
