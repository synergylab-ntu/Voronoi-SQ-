system('mkoctfile --mex src/bslv_algs.c src/bslv_lists.c src/bslv_lp.c src/bslv_main.c src/bslv_poly.c src/bslv_vlp.c -O3 -o bensolve -lglpk');
system('rm *.o');
