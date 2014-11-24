autoPar-demo
============

Genenerated files for autoPar

The corresponding input files are located in rose/projects/autoParallelization/tests 

https://github.com/rose-compiler/edg4x-rose/tree/master/projects/autoParallelization/tests

Several types of files are generated during "make check" of autoPar.

rose_?*.c/C: the generated file, with OpenMP directives inserted if feasible.

*.out: the screen output, explaining why some loops cannot be parallized when possible.

*.patch: the patch files generated

*.diff: the diff output, comparing generated files to reference files. This should be empty.

