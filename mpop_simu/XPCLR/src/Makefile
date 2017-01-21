# "make xp-clr" to make XP-CLR program
# "make install" to copy all programs to ../bin
# "make clean" to clean up extra files in the directory


HOMEL=$(PWD)
TOP=../bin

CC=cl
CFLAGS=-c 
LIBDIR=


all:	XPCLR $(LIBDIR)p2sweep.obj $(LIBDIR)hualib.obj

XPCLR:	xpclr.obj $(LIBDIR)p2sweep.obj $(LIBDIR)hualib.obj xpclr.h
	$(CC) $(LIBDIR)p2sweep.obj $(LIBDIR)hualib.obj xpclr.obj -o XPCLR

xpclr.obj:	xpclr.cpp xpclr.h
	$(CC) $(CFLAGS) xpclr.cpp

p2sweep.obj:	$(LIBDIR)p2sweep.cpp $(LIBDIR)p2sweep.h
	$(CC) $(CFLAGS) $(LIBDIR)p2sweep.cpp 

hualib.obj:	$(LIBDIR)hualib.cpp $(LIBDIR)hualib.h
	$(CC) $(CFLAGS) $(LIBDIR)hualib.cpp

clean:
	rm -f *.obj

install:	XPCLR
	cp XPCLR $(TOP)
