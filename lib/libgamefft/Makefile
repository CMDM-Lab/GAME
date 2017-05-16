LIBNAME=libgamefft
VERSION=1.0
LIBDIR=/usr/local/lib

all:
	make gamefft
	
clean:
	rm *.o *.so.*

gamefft:
	g++ -fPIC -c DPKernelCpuCcOpt.cc 
	g++ -O3 -shared -o $(LIBNAME).so.$(VERSION) *.o
	
install:
	cp $(LIBNAME).so.$(VERSION) $(LIBDIR)
	ln -s $(LIBDIR)/$(LIBNAME).so.$(VERSION) $(LIBDIR)/$(LIBNAME).so
	rm *.o
