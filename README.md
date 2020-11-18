# NTRUPrimePolyMultF3

This repository contains C implementation of the new efficient hybrid polynomial multiplication algorithms over Z<sub>3</sub> that are first described in: "Faster Characteristic Three Polynomial Multiplication and Its Application to NTRU Prime Decapsulation" by E. Yeniaras and M. Cenk.  The repository contains the C implementation of the following nine hybrid polynomial multiplication algorithms: Hybrid-1, Hybrid-2, N1-Hybrid, N1-Hybrid2, N2-Hybrid, A3-Hybrid, A3-Hybrid2, V1-Hybrid, LT-Hybrid. The hybrid algorithms are constructed as different combinations of the newly proposed N1, N2, and V1 algorithms along with the other well-known algorithms such as KA2, A2, A3, UB, LT, SB. One can refer to the paper: https://eprint.iacr.org/2020/1336.pdf for the details of these new algorithms. We add sample input files corresponding to the polynomial sizes n=761, n=765, and n=768. 

Note that we execute each algorithm 100.000 times and get the median of the executions to compare the cycle counts. Benchmark tests for all implementations are performed on an Intel (R) Core (TM) i7-9750H processor running at 2600 MHz. The operating system is Ubuntu 20.04.1 LTS and Linux Kernel 5.4.0. All software is compiled with gcc-9.3.0. We apply our new efficient hybrid algorithms to the NTRU Prime Protocol which is a quantum-resistant, lattice-based Key Encapsulation Mechanism (KEM) submitted to NIST PQC Competition and advanced to the third round as an alternative candidate by Bernstein et al.





# Contributors

Esra Yeniaras <esramath@gmail.com>




# License and Copyright 

Licensed under the [GNU AFFERO GENERAL PUBLIC LICENSE](LICENSE) 


