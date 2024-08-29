system('mkoctfile --mex ./src/*.c -O3 -o bensolve -lglpk');
system('rm *.o');