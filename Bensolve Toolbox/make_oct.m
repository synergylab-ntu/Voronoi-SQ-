% call this script to genereate a mex file for Octace on a Windows platform, see the manual for more detailed instructions
system('mkoctfile --mex src/bslv_algs.c src/bslv_lists.c src/bslv_lp.c src/bslv_main.c src/bslv_poly.c src/bslv_vlp.c -lglpk -O3 -o bensolve');
system('rm *.o');