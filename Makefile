obj = main.o acor.o data.o evidence.o exception.o exoplanetjd.o exoplanet_hyperpara.o exoplanet_init.o keplers_eqn.o int2str.o level.o max.o mean.o model.o quicksort.o randn.o rng.o sampling.o sampling_DNest.o uniformtest.o var.o
exe = MCMC
cmp = icpc

#cflag = -g #for debugging
cflag = -O3 #for running

$(exe): $(obj)
	$(cmp) $(cflag) -o $(exe) $(obj)  -llapack -lblas

%.o: %.cpp Makefile preprocessors.h
	$(cmp) $(cflag) -c $<

clean:
	rm $(exe) $(obj)
