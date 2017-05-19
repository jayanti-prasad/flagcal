AUTOMAKE_OPTIONS= subdir-objects
ACLOCAL_AMFLAGS= ${ACLOCAL_FLAGS}  

CC=gcc -fopenmp 


bin_PROGRAMS = flagcal liststruc readfitshdr readfits plotbadant  plotbadchan plotbadbase igplot ivplot igrplot ivsrplot

flagcal_SOURCES=  average.c flag_apply.c  median.c complex.c minmax.c mad.c\
	  peak_remove.c piksrt.c longest_block.c group_data.c printerror.c io.c\
	 read_fits.c initilise.c smooth_data.c  mad_filter.c flag_data.c gain_solve.c\
	 gain_compute.c remove_rfi.c gain_interpolate.c gain_apply.c  flux_calib.c\
	 write_fits.c  main.c flag_timerange.c find_bad_timerange.c gettime.c

flagcal_CPPFLAGS=$(CPPFLAGS)
flagcal_LDFLAGS=-lcfitsio  -lm

liststruc_SOURCES= liststruc.c
liststruc_CPPFLAGS=$(CPPFLAGS)
liststruc_LDFLAGS=-lcfitsio  -lm

readfitshdr_SOURCES=readfitshdr.c
readfitshdr_CPPFLAGS=$(CPPFLAGS)
readfitshdr_LDFLAGS=-lcfitsio  -lm

readfits_SOURCES=gettime.c complex.c printerror.c group_data.c jdcal.c readfits.c
readfits_CPPFLAGS=$(CPPFLAGS)
readfits_LDFLAGS=-lcfitsio  -lm

plotbadant_SOURCES=minmax.c plotbadant.c 
plotbadant_CPPFLAGS=$(CPPFLAGS)
plotbadant_LDFLAGS=-lcpgplot -lpgplot -lpng -lgfortran -lX11 -lm

plotbadchan_SOURCES=minmax.c  plotbadchan.c 
plotbadchan_CPPFLAGS=$(CPPFLAGS)
plotbadchan_LDFLAGS=-lcpgplot -lpgplot  -lpng -lgfortran -lX11 -lm

plotbadbase_SOURCES=minmax.c  plotbadbase.c
plotbadbase_CPPFLAGS=$(CPPFLAGS)
plotbadbase_LDFLAGS=-lcpgplot -lpgplot -lpng -lgfortran -lX11 -lm

igplot_SOURCES=jdcal.c complex.c  minmax.c igplot.c
igplot_CPPFLAGS=$(CPPFLAGS)
igplot_LDFLAGS=-lcpgplot -lpgplot -lpng -lgfortran -lX11 -lm

ivsrplot_SOURCES=jdcal.c io.c complex.c baselines.c  minmax.c ivsrplot.c
ivsrplot_CPPFLAGS=$(CPPFLAGS)
ivsrplot_LDFLAGS=-lcpgplot -lpgplot -lpng -lgfortran -lX11 -lm

ivplot_SOURCES=jdcal.c minmax.c baselines.c  ivplot.c
ivplot_CPPFLAGS=$(CPPFLAGS)
ivplot_LDFLAGS=-lcpgplot -lpgplot -lpng -lgfortran -lX11 -lm

igrplot_SOURCES=gettime.c complex.c printerror.c baselines.c minmax.c read_fits.c jdcal.c   group_data.c igrplot.c
igrplot_CPPFLAGS=$(CPPFLAGS)
igrplot_LDFLAGS=-lcpgplot -lpgplot -lpng -lgfortran -lX11 -lm


