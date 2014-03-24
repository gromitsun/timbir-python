# Edit the following variables as needed
HDF_INSTALL = /opt/local/include
EXTLIB = -L$(HDF_INSTALL)/lib
CC          = gcc
CFLAGS      = 
LIB         = -lsz -lz -lm

INCLUDE   = -I$(HDF_INSTALL)/include
LIBSHDF   = $(EXTLIB) $(HDF_INSTALL)/lib/libhdf5.a 

all: h5_crtdat \
     h5_rdwt \
     h5_crtatt \
     h5_crtgrp \
     h5_crtgrpar \
     h5_crtgrpd \
     XT_Engine \
     XT_AMatrix \
     XT_ICD_Update \
     XT_IOMisc \
     XT_Init \
     XT_Profile \
     XT_genSinogram \
     allocate \
     randlib \
     tiff \ 

h5_crtdat: h5_crtdat.c
	$(CC) $(CFLAGS) -o $@ h5_crtdat.c $(INCLUDE) $(LIBSHDF) $(LIB)

h5_rdwt: h5_rdwt.c
	$(CC) $(CFLAGS) -o $@ h5_rdwt.c $(INCLUDE) $(LIBSHDF) $(LIB)

h5_crtatt: h5_crtatt.c
	$(CC) $(CFLAGS) -o $@ h5_crtatt.c $(INCLUDE) $(LIBSHDF) $(LIB)

h5_crtgrp: h5_crtgrp.c
	$(CC) $(CFLAGS) -o $@ h5_crtgrp.c $(INCLUDE) $(LIBSHDF) $(LIB)

h5_crtgrpar: h5_crtgrpar.c
	$(CC) $(CFLAGS) -o $@ h5_crtgrpar.c $(INCLUDE) $(LIBSHDF) $(LIB)

h5_crtgrpd: h5_crtgrpd.c
	$(CC) $(CFLAGS) -o $@ h5_crtgrpd.c $(INCLUDE) $(LIBSHDF) $(LIB)

XT_Engine: XT_Engine.c
	$(CC) $(CFLAGS) -o $@ XT_Engine.c $(INCLUDE) $(LIBSHDF) $(LIB)

XT_AMatrix: XT_AMatrix.c
	$(CC) $(CFLAGS) -o $@ XT_AMatrix.c $(INCLUDE) $(LIBSHDF) $(LIB)

XT_ICD_Update: XT_ICD_Update.c
	$(CC) $(CFLAGS) -o $@ XT_ICD_Update.c $(INCLUDE) $(LIBSHDF) $(LIB)

XT_IOMisc: XT_IOMisc.c
	$(CC) $(CFLAGS) -o $@ XT_IOMisc.c $(INCLUDE) $(LIBSHDF) $(LIB)

XT_Init: XT_Init.c
	$(CC) $(CFLAGS) -o $@ XT_Init.c $(INCLUDE) $(LIBSHDF) $(LIB)

XT_Profile: XT_Profile.c
	$(CC) $(CFLAGS) -o $@ XT_Profile.c $(INCLUDE) $(LIBSHDF) $(LIB)

XT_genSinogram: XT_genSinogram.c
	$(CC) $(CFLAGS) -o $@ XT_genSinogram.c $(INCLUDE) $(LIBSHDF) $(LIB)

allocate: allocate.c
	$(CC) $(CFLAGS) -o $@ allocate.c $(INCLUDE) $(LIBSHDF) $(LIB)

randlib: randlib.c
	$(CC) $(CFLAGS) -o $@ randlib.c $(INCLUDE) $(LIBSHDF) $(LIB)

tiff: tiff.c
	$(CC) $(CFLAGS) -o $@ tiff.c $(INCLUDE) $(LIBSHDF) $(LIB)

clean: 
	rm -f *.h5 *.o \
        h5_crtdat \
        h5_rdwt \
        h5_crtatt \
        h5_crtgrp \
        h5_crtgrpar \
        h5_crtgrpd \
 

.SUFFIXES:.o.c

