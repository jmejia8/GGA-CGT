# GGA-CGT

## GROUPING GENETIC ALGORITHM WITH CONTROLLED GENE TRANSMISSION FOR THE BIN PACKING PROBLEM

The program excecutes GGA-CGT over a set instances using different configurations
given by the user. Each configuration represents an independent execution of the GA.

Author: 

- Marcela Quiroz-Castellanos
- qc.marcela@gmail.com
- Tecnológico Nacional de México
- Instituto Tecnológico de Ciudad Madero
- División de Estudios de Posgrado e Investigación
- Depto. de Sistemas y Computación



### Input:
File `instances.txt` including the name of the BPP instances to be solved;
Files including the standard instances to be solve;
File `configurations.txt` including the parameter values for each experiment;

### Output:
A set of files `GGA-CGT_(i).txt` including the experimental results for each
configuration i in the input, stored in directory: Solutions_GGA-CGT;
If(`save_bestSolution = 1`) a set of files `HGGA_S_(i)_instance.txt` including the
obtained solution for each instance, for each configuration i, stored in directory: `Details_GGA-CGT`.

### Reference:

```
Quiroz-Castellanos, M., Cruz-Reyes, L., Torres-Jimenez, J.,
Gómez, C., Huacuja, H. J. F., & Alvim, A. C. (2015).
  A grouping genetic algorithm with controlled gene transmission for
  the bin packing problem. Computers & Operations Research, 55, 52-64.
```
