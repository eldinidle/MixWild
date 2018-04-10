echo 'Creating base folder'
mkdir "Windows32 Binaries"
echo 'Compiling and Archiving'
echo 'mixreg'
cd mixreg
del libmix.a
i686-w64-mingw32-gfortran -c -o mixlib.o -g -static -m32 -D _WIN32 -cpp mixlib.f90
i686-w64-mingw32-gfortran -c -o matcal1.o -g -static -m32 -D _WIN32 -cpp matcal1.f90
i686-w64-mingw32-gfortran -c -o matcal2.o -g -static -m32 -D _WIN32 -cpp matcal2.f90
i686-w64-mingw32-gfortran -c -o rrmset.o -g -static -m32 -D _WIN32 -cpp rrmset.f90
i686-w64-mingw32-gfortran -c -o mixregc.o -g -static -m32 -D _WIN32 -cpp mixregc.f90
i686-w64-mingw32-gfortran -c -o dllstub.o -g -static -m32 -D _WIN32 -cpp dllstub.f90
mkdir modules
move *.mod modules 
ar -qs libmix.a matcal1.o matcal2.o mixlib.o rrmset.o mixregc.o dllstub.o
i686-w64-mingw32-gfortran libmix.a -Jmodules -LC:\msys64\mingw32\lib -static -static -static-libgcc -static-libgfortran -m32 -D _WIN32 -cpp -o mixreg
cd ..
move "mixreg\mixreg.exe" "Windows32 Binaries"
echo 'repeat_mixreg'
cd repeat_mixreg 
i686-w64-mingw32-gfortran -c -o repeat_mixreg.o -g -static -m32 -D _WIN32 -cpp repeat_mixreg.f90
i686-w64-mingw32-gfortran repeat_mixreg.o  -LC:\msys64\mingw32\lib -static -static -static-libgcc -static-libgfortran  -m32 -D _WIN32 -cpp -o repeat_mixreg
cd ..
move repeat_mixreg\repeat_mixreg.exe "Windows32 Binaries"
echo 'mix_random'
cd mix_random
i686-w64-mingw32-gfortran -c -o mix_random.o -g -static -m32 -D _WIN32 -cpp mix_random.f90
i686-w64-mingw32-gfortran mix_random.o  -LC:\msys64\mingw32\lib  -static -static-libgcc -static-libgfortran  -m32 -D _WIN32 -cpp -o mix_random
cd ..
move mix_random\mix_random.exe "Windows32 Binaries"
echo 'mixregls_random_mixreg'
cd mixregls_random_mixreg 
del libmix.a 
i686-w64-mingw32-gfortran -c -o hermite_rule.o -g -static -m32 -D _WIN32 -cpp hermite_rule.f90
i686-w64-mingw32-gfortran -c -o mixregls_random_mixreg.o -g -static -m32 -D _WIN32 -cpp mixregls_random_mixreg.f90
i686-w64-mingw32-gfortran -c -o sstar.o -g -static -m32 -D _WIN32 -cpp sstar.f90
ar -qs libmix.a hermite_rule.o sstar.o mixregls_random_mixreg.o 
i686-w64-mingw32-gfortran libmix.a -LC:\msys64\mingw32\lib -static -static-libgcc -static-libgfortran -m32 -D _WIN32 -cpp -o mixregls_random_mixreg
cd ..
move mixregls_random_mixreg\mixregls_random_mixreg.exe "Windows32 Binaries"
echo 'mixregls_random_mixor'
cd mixregls_random_mixor 
del libmix.a 
i686-w64-mingw32-gfortran -c -o hermite_rule.o -g -static -m32 -D _WIN32 -cpp hermite_rule.f90
i686-w64-mingw32-gfortran -c -o sstar.o -g -static -m32 -D _WIN32 -cpp sstar.f90
i686-w64-mingw32-gfortran -c -o mixregls_random_mixor.o -g -static -m32 -D _WIN32 -cpp mixregls_random_mixor.f90
ar -qs libmix.a hermite_rule.o sstar.o mixregls_random_mixor.o 
i686-w64-mingw32-gfortran libmix.a -LC:\msys64\mingw32\lib -static -static-libgcc -static-libgfortran -m32 -D _WIN32 -cpp -o mixregls_random_mixor
cd..
move mixregls_random_mixor\mixregls_random_mixor.exe "Windows32 Binaries"
echo 'repeat_mixor'
cd repeat_mixor 
i686-w64-mingw32-gfortran -c -o repeat_mixor.o -g -static -m32 -D _WIN32 -cpp repeat_mixor.f90
i686-w64-mingw32-gfortran repeat_mixor.o  -LC:\msys64\mingw32\lib  -static -static-libgcc -static-libgfortran  -m32 -D _WIN32 -cpp -o repeat_mixor
cd ..
move repeat_mixor\repeat_mixor.exe "Windows32 Binaries"
echo 'mixregmls_random_mixreg'
cd mixregmls_random_mixreg 
del libmix.a 
i686-w64-mingw32-gfortran -c -o hermite_rule.o -g -static -m32 -D _WIN32 -cpp hermite_rule.f90
i686-w64-mingw32-gfortran -c -o mixregmls_random_mixreg.o -g -static -m32 -D _WIN32 -cpp mixregmls_random_mixreg.f90
i686-w64-mingw32-gfortran -c -o sstar.o -g -static -m32 -D _WIN32 -cpp sstar.f90
ar -qs libmix.a hermite_rule.o sstar.o mixregmls_random_mixreg.o 
i686-w64-mingw32-gfortran libmix.a -LC:\msys64\mingw32\lib -static -static-libgcc -static-libgfortran -m32 -D _WIN32 -cpp -o mixregmls_random_mixreg
cd ..
move mixregmls_random_mixreg\mixregmls_random_mixreg.exe "Windows32 Binaries"
echo 'mixregmls_random_mixor'
cd mixregmls_random_mixor 
del libmix.a 
i686-w64-mingw32-gfortran -c -o hermite_rule.o -g -static -m32 -D _WIN32 -cpp hermite_rule.f90
i686-w64-mingw32-gfortran -c -o mixregmls_random_mixor.o -g -static -m32 -D _WIN32 -cpp mixregmls_random_mixor.f90
i686-w64-mingw32-gfortran -c -o sstar.o -g -static -m32 -D _WIN32 -cpp sstar.f90
ar -qs libmix.a hermite_rule.o sstar.o mixregmls_random_mixor.o 
i686-w64-mingw32-gfortran libmix.a -LC:\msys64\mingw32\lib -static -static-libgcc -static-libgfortran -m32 -D _WIN32 -cpp -o mixregmls_random_mixor
cd ..
move mixregmls_random_mixor\mixregmls_random_mixor.exe "Windows32 Binaries"
echo 'mixor'
cd mixor 
del libmix.a
i686-w64-mingw32-gfortran -c -o mixlib.o -g -static -m32 -D _WIN32 -cpp mixlib.f90
i686-w64-mingw32-gfortran -c -o sstar.o -g -static -m32 -D _WIN32 -cpp sstar.f90
i686-w64-mingw32-gfortran -c -o mixord_chol.o -g -static -m32 -D _WIN32 -cpp mixord_chol.f90
i686-w64-mingw32-gfortran -c -o dllstub.o -g -static -m32 -D _WIN32 -cpp dllstub.f90
mkdir modules
move *.mod modules 
ar -qs libmix.a mixlib.o  sstar.o mixord_chol.o dllstub.o
i686-w64-mingw32-gfortran libmix.a -Jmodules -LC:\msys64\mingw32\lib -static -static-libgcc -static-libgfortran -m32 -D _WIN32 -cpp -o mixor
cd .. 
move mixor\mixor.exe "Windows32 Binaries"
echo 'all done'