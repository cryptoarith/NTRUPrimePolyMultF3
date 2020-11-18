# NTRUPrimePolyMultF3

This repository contains C implementation of the new efficient hybrid polynomial multiplication algorithms over Z<sub>3</sub> that are first described  "Faster Characteristic Three Polynomial Multiplication and Its Application to NTRU Prime Decapsulation" by E. Yeniaras and M. Cenk. One can refer to the paper: https://eprint.iacr.org/2020/1336.pdf for the details of the algorithms. We add sample input files corresponding to the polynomial sizes n=761, n=765, and n=768.

Benchmark tests for all implementations are performed on an Intel (R) Core (TM) i7-9750H processor running at 2600 MHz. The operating system is Ubuntu 20.04.1 LTS and Linux Kernel 5.4.0. All software is compiled with gcc-9.3.0.

Note that we execute each algorithm 100.000 times and get the median of the executions to compare the cycle counts. 


# Contributors

Esra Yeniaras <esramath@gmail.com>




# License and Copyright 

Licensed under the [GNU AFFERO GENERAL PUBLIC LICENSE](LICENSE) 


