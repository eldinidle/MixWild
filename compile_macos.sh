!/bin/bash

Must delete libmath for static to compile


echo 'Creating base folder'
mkdir 'Mac Binaries'
echo 'Compiling and Archiving'

echo 'mixregls_both'
cd lsboth_random_mixblank
gfortran -c -o amod_mls_sstar.o -g -static -D _UNIX -cpp amod_mls_sstar.f90
gfortran -c -o hermite_rule.o -g -static -D _UNIX -cpp hermite_rule.f90
gfortran -c -o lsboth_random_mixblank.o -g -static -D _UNIX -cpp lsboth_random_mixblank.f90
mkdir modules
mv *.mod modules
ar -qs libmix.a amod_mls_sstar.o hermite_rule.o lsboth_random_mixblank.o
gfortran libmix.a -Jmodules -L/usr/local/gfortran/lib -static-libgcc -static-libgfortran -D _UNIX -cpp -o lsboth_random_mixblank
rm *.o
rm -rf modules
rm libmix.a
cd ..
mv lsboth_random_mixblank/lsboth_random_mixblank 'Mac Binaries'

echo 'mixreg'
cd mixreg
gfortran -c -o mixlib_mixreg.o -g -static -D _UNIX -cpp mixlib_mixreg.f90
gfortran -c -o matcal1.o -g -static -D _UNIX -cpp matcal1.f90
gfortran -c -o matcal2.o -g -static -D _UNIX -cpp matcal2.F90
gfortran -c -o rrmset.o -g -static -D _UNIX -cpp rrmset.f90
gfortran -c -o mixregc.o -g -static -D _UNIX -cpp mixregc.f90
gfortran -c -o dllstub.o -g -static -D _UNIX -cpp dllstub.f90
mkdir modules
mv *.mod modules
ar -qs libmix.a matcal1.o matcal2.o mixlib_mixreg.o rrmset.o mixregc.o dllstub.o
gfortran libmix.a -Jmodules -L/usr/local/gfortran/lib -static-libgcc -static-libgfortran -D _UNIX -cpp -o mixreg
rm *.o
rm -rf modules
rm libmix.a
cd ..
mv mixreg/mixreg 'Mac Binaries'

echo 'mixors_both'
cd mixors_random_mixblank
gfortran -c -o hermite_rule.o -g -static -D _UNIX -cpp hermite_rule.f90
gfortran -c -o mixors_random_mixblank.o -g -static -D _UNIX -cpp mixors_random_mixblank.f90
gfortran -c -o sstar2d.o -g -static -D _UNIX -cpp sstar2d.f90
ar -qs libmix.a hermite_rule.o mixors_random_mixblank.o sstar2d.o
mkdir modules
mv *.mod modules
gfortran libmix.a -Jmodules -L/usr/local/gfortran/lib -static-libgcc -static-libgfortran -D _UNIX -cpp -o mixors_random_mixblank
rm *.o
rm -rf modules
rm libmix.a
cd ..
mv mixors_random_mixblank/mixors_random_mixblank 'Mac Binaries'

echo 'mixors'
cd mixors
gfortran -c -o hermite_rule.o -g -static -D _UNIX -cpp hermite_rule.f90
gfortran -c -o mixors_simple.o -g -static -D _UNIX -cpp mixors_simple.f90
gfortran -c -o sstar2d.o -g -static -D _UNIX -cpp sstar2d.f90
ar -qs libmix.a hermite_rule.o mixors_simple.o sstar2d.o
mkdir modules
mv *.mod modules
gfortran libmix.a -Jmodules -L/usr/local/gfortran/lib -static-libgcc -static-libgfortran -D _UNIX -cpp -o mixors
rm *.o
rm -rf modules
rm libmix.a
cd ..
mv mixors/mixors 'Mac Binaries'

echo 'mixno'
cd mixno
gfortran -c -o MIXLIB.O -g -static -D _UNIX -cpp MIXLIB.F90
gfortran -c -o mixnoc.o -g -static -D _UNIX -cpp mixnoc.f90
ar -qs libmix.a MIXLIB.O mixnoc.o
mkdir modules
mv *.mod modules
gfortran libmix.a -Jmodules -L/usr/local/gfortran/lib -static-libgcc -static-libgfortran -D _UNIX -cpp -o mixno
rm *.o
rm *.O
rm -rf modules
rm libmix.a
cd ..
mv mixno/mixno 'Mac Binaries'

echo 'mixpreg'
cd mixpreg
gfortran -c -o MIXLIB.O -g -static -D _UNIX -cpp MIXLIB.F90
gfortran -c -o mixpreg.o -g -static -D _UNIX -cpp mixpreg.f90
ar -qs libmix.a MIXLIB.O mixpreg.o
mkdir modules
mv *.mod modules
gfortran libmix.a -Jmodules -L/usr/local/gfortran/lib -static-libgcc -static-libgfortran -D _UNIX -cpp -o mixpreg
rm *.o
rm *.O
rm -rf modules
rm libmix.a
cd ..
mv mixpreg/mixpreg 'Mac Binaries'

echo 'stage2only'
cd stage2only
gfortran -c -o amod_mls_sstar.o -g -static -D _UNIX -cpp amod_mls_sstar.f90
gfortran -c -o stage2only.o -g -static -D _UNIX -cpp stage2only.f90
ar -qs libmix.a amod_mls_sstar.o stage2only.o
mkdir modules
mv *.mod modules
gfortran libmix.a -Jmodules -L/usr/local/gfortran/lib -static-libgcc -static-libgfortran -D _UNIX -cpp -o stage2only
rm *.o
rm -rf modules
rm libmix.a
cd ..
mv stage2only/stage2only 'Mac Binaries'


echo 'all done'
