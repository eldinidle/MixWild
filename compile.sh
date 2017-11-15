#!/bin/bash
#
# Must delete libmath for static to compile 
#
#
echo 'mixreg'
cd mixreg
rm libmix.a
gfortran -c -o mixlib.o -g -static mixlib.f90
gfortran -c -o matcal1.o -g -static matcal1.f90
gfortran -c -o matcal2.o -g -static matcal2.f90
gfortran -c -o rrmset.o -g -static rrmset.f90
gfortran -c -o mixregc.o -g -static  mixregc.f90
gfortran -c -o dllstub.o -g -static dllstub.f90
mv *.mod modules 
ar -qs libmix.a matcal1.o matcal2.o mixlib.o rrmset.o mixregc.o dllstub.o
gfortran libmix.a -Jmodules -L/usr/local/gfortran/lib -static-libgcc -static-libgfortran -o mixreg
cd ..
echo 'repeat_mixreg'
cd repeat_mixreg 
gfortran -c -o repeat_mixreg.o -g -static repeat_mixreg.f90
gfortran repeat_mixreg.o  -L/usr/local/gfortran/lib  -static-libgcc -static-libgfortran  -o repeat_mixreg
cd ..
echo 'mix_random'
cd mix_random
gfortran -c -o mix_random.o -g -static mix_random.f90
gfortran mix_random.o  -L/usr/local/gfortran/lib  -static-libgcc -static-libgfortran  -o mix_random
cd ..
echo 'mixregls_random_mixreg'
cd mixregls_random_mixreg 
rm libmix.a 
gfortran -c -o hermite_rule.o -g -static hermite_rule.f90
gfortran -c -o mixregls_random_mixreg.o -g -static mixregls_random_mixreg.f90
gfortran -c -o sstar.o -g -static sstar.f90
ar -qs libmix.a hermite_rule.o sstar.o mixregls_random_mixreg.o 
gfortran libmix.a -L/usr/local/gfortran/lib -static-libgcc -static-libgfortran -o mixregls_random_mixreg
cd ..
echo 'mixregls_random_mixor'
cd mixregls_random_mixor 
rm libmix.a 
gfortran -c -o hermite_rule.o -g -static hermite_rule.f90
gfortran -c -o sstar.o -g -static sstar.f90
gfortran -c -o mixregls_random_mixor.o -g -static mixregls_random_mixor.f90
ar -qs libmix.a hermite_rule.o sstar.o mixregls_random_mixor.o 
gfortran libmix.a -L/usr/local/gfortran/lib -static-libgcc -static-libgfortran -o mixregls_random_mixor
cd..
echo 'repeat_mixor'
cd repeat_mixor 
gfortran -c -o repeat_mixor.o -g -static repeat_mixor.f90
gfortran repeat_mixor.o  -L/usr/local/gfortran/lib  -static-libgcc -static-libgfortran  -o repeat_mixor
cd ..
echo 'mixregmls_random_mixreg'
cd mixregmls_random_mixreg 
rm libmix.a 
gfortran -c -o hermite_rule.o -g -static hermite_rule.f90
gfortran -c -o mixregmls_random_mixreg.o -g -static mixregmls_random_mixreg.f90
gfortran -c -o sstar.o -g -static sstar.f90
ar -qs libmix.a hermite_rule.o sstar.o mixregmls_random_mixreg.o 
gfortran libmix.a -L/usr/local/gfortran/lib -static-libgcc -static-libgfortran -o mixregmls_random_mixreg
cd ..
echo 'mixregls_random_mixor'
cd mixregmls_random_mixor 
rm libmix.a 
gfortran -c -o hermite_rule.o -g -static hermite_rule.f90
gfortran -c -o mixregmls_random_mixor.o -g -static mixregmls_random_mixor.f90
gfortran -c -o sstar.o -g -static sstar.f90
ar -qs libmix.a hermite_rule.o sstar.o mixregmls_random_mixor.o 
gfortran libmix.a -L/usr/local/gfortran/lib -static-libgcc -static-libgfortran -o mixregmls_random_mixor
cd ..
echo 'mixor'
cd mixor 
rm libmix.a
gfortran -c -o mixlib.o -g -static mixlib.f90
gfortran -c -o sstar.o -g -static sstar.f90
gfortran -c -o mixord_chol.o -g -static mixord_chol.f90
gfortran -c -o dllstub.o -g -static dllstub.f90
cp *.mod modules 
ar -qs libmix.a mixlib.o  sstar.o mixord_chol.o dllstub.o
gfortran libmix.a -Jmodules -L/usr/local/gfortran/lib -static-libgcc -static-libgfortran -o mixor
cd .. 
echo 'all done'