
 rm fs.bin
 echo 'Compile...'
 gfortran -O3 constants.f90 fs_sch.f90 main_fs.f90 -o fs
 echo 'Solving Falkner-Skan.'
 ./fs 
 echo 'Saving data. Finished.'
 rm fs

