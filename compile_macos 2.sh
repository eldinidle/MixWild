#!/bin/bash
#
# Must delete libmath for static to compile
#
#
echo 'Creating base folder'
mkdir 'Mac Binaries'
echo 'Compiling and Archiving'
echo 'mixreg'
cd mixreg
rm libmix.a
gfortran -c -o mixlib.o -g -static -D _UNIX -cpp mixlib.f90
gfortran -c -o matcal1.o -g -static -D _UNIX -cpp matcal1.f90
gfortran -c -o matcal2.o -g -static -D _UNIX -cpp matcal2.f90
gfortran -c -o rrmset.o -g -static  -D _UNIX -cpp rrmset.f90
gfortran -c -o mixregc.o -g -static  -D _UNIX -cpp mixregc.f90
gfortran -c -o dllstub.o -g -static -D _UNIX -cpp dllstub.f90
mv *.mod modules
ar -qs libmix.a matcal1.o matcal2.o mixlib.o rrmset.o mixregc.o dllstub.o
gfortran libmix.a -Jmodules -L/usr/local/gfortran/lib -static-libgcc -static-libgfortran -D _UNIX -cpp -o mixreg
cd ..
mv mixreg/mixreg 'Mac Binaries'
echo 'repeat_mixreg'
cd repeat_mixreg
gfortran -c -o repeat_mixreg.o -g -static -D _UNIX -cpp repeat_mixreg.f90
gfortran repeat_mixreg.o  -L/usr/local/gfortran/lib  -static-libgcc -static-libgfortran  -D _UNIX -cpp -o repeat_mixreg
cd ..
mv repeat_mixreg/repeat_mixreg 'Mac Binaries'
echo 'mix_random'
cd mix_random
gfortran -c -o mix_random.o -g -static -D _UNIX -cpp mix_random.f90
gfortran mix_random.o  -L/usr/local/gfortran/lib  -static-libgcc -static-libgfortran -D _UNIX -cpp -o mix_random
cd ..
mv mix_random/mix_random 'Mac Binaries'
echo 'mixregls_random_mixreg'
cd mixregls_random_mixreg
rm libmix.a
gfortran -c -o hermite_rule.o -g -static -D _UNIX -cpp hermite_rule.f90
gfortran -c -o mixregls_random_mixreg.o -g -static -D _UNIX -cpp mixregls_random_mixreg.f90
gfortran -c -o sstar.o -g -static -D _UNIX -cpp sstar.f90
ar -qs libmix.a hermite_rule.o sstar.o mixregls_random_mixreg.o
gfortran libmix.a -L/usr/local/gfortran/lib -static-libgcc -static-libgfortran -D _UNIX -cpp -o mixregls_random_mixreg
cd ..
mv mixregls_random_mixreg/mixregls_random_mixreg 'Mac Binaries'
echo 'mixregls_random_mixor'
cd mixregls_random_mixor
rm libmix.a
gfortran -c -o hermite_rule.o -g -static -D _UNIX -cpp hermite_rule.f90
gfortran -c -o sstar.o -g -static -D _UNIX -cpp sstar.f90
gfortran -c -o mixregls_random_mixor.o -g -static -D _UNIX -cpp mixregls_random_mixor.f90
ar -qs libmix.a hermite_rule.o sstar.o mixregls_random_mixor.o
gfortran libmix.a -L/usr/local/gfortran/lib -static-libgcc -static-libgfortran -D _UNIX -cpp -o mixregls_random_mixor
cd ..
mv mixregls_random_mixor/mixregls_random_mixor 'Mac Binaries'
echo 'repeat_mixor'
cd repeat_mixor
gfortran -c -o repeat_mixor.o -g -static -D _UNIX -cpp repeat_mixor.f90
gfortran repeat_mixor.o  -L/usr/local/gfortran/lib  -static-libgcc -static-libgfortran -D _UNIX -cpp  -o repeat_mixor
cd ..
mv repeat_mixor/repeat_mixor 'Mac Binaries'
echo 'mixregmls_random_mixreg'
cd mixregmls_random_mixreg
rm libmix.a
gfortran -c -o hermite_rule.o -g -static -D _UNIX -cpp hermite_rule.f90
gfortran -c -o mixregmls_random_mixreg.o -g -static -D _UNIX -cpp mixregmls_random_mixreg.f90
gfortran -c -o sstar.o -g -static -D _UNIX -cpp sstar.f90
ar -qs libmix.a hermite_rule.o sstar.o mixregmls_random_mixreg.o
gfortran libmix.a -L/usr/local/gfortran/lib -static-libgcc -static-libgfortran -D _UNIX -cpp -o mixregmls_random_mixreg
cd ..
mv mixregmls_random_mixreg/mixregmls_random_mixreg 'Mac Binaries'
echo 'mixregmls_random_mixor'
cd mixregmls_random_mixor
rm libmix.a
gfortran -c -o hermite_rule.o -g -static -D _UNIX -cpp hermite_rule.f90
gfortran -c -o mixregmls_random_mixor.o -g -static -D _UNIX -cpp mixregmls_random_mixor.f90
gfortran -c -o sstar.o -g -static -D _UNIX -cpp sstar.f90
ar -qs libmix.a hermite_rule.o sstar.o mixregmls_random_mixor.o
gfortran libmix.a -L/usr/local/gfortran/lib -static-libgcc -static-libgfortran -D _UNIX -cpp -o mixregmls_random_mixor
cd ..
mv mixregmls_random_mixor/mixregmls_random_mixor 'Mac Binaries'
echo 'mixor'
cd mixor
rm libmix.a
gfortran -c -o mixlib.o -g -static -D _UNIX -cpp mixlib.f90
gfortran -c -o sstar.o -g -static -D _UNIX -cpp sstar.f90
gfortran -c -o mixord_chol.o -g -static -D _UNIX -cpp mixord_chol.f90
gfortran -c -o dllstub.o -g -static -D _UNIX -cpp dllstub.f90
cp *.mod modules
ar -qs libmix.a mixlib.o  sstar.o mixord_chol.o dllstub.o
gfortran libmix.a -Jmodules -L/usr/local/gfortran/lib -static-libgcc -static-libgfortran -D _UNIX -cpp -o mixor
cd ..
mv mixor/mixor 'Mac Binaries'
echo 'all done'
