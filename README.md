# GGA-CGT

## Grouping Genetic Algorithm with Controlled Gene Transmission for the Bin Packing Problem

The program executes GGA-CGT over a set instances using different configurations
given by the user. Each configuration represents an independent execution of the GA.

Author: 

- Marcela Quiroz-Castellanos
- qc.marcela@gmail.com
- Tecnológico Nacional de México
- Instituto Tecnológico de Ciudad Madero
- División de Estudios de Posgrado e Investigación
- Depto. de Sistemas y Computación

Uploaded and documented by [jmejia8](https://github.com/jmejia8/).


### My first time (Linux):
To build this program just compile by running `make` . Use the following command to
run the instance `instances/BPP1.txt`:

Running algorithm with default parameters (settings):

```shell
./GGA-CGT --instance instances/BPP1.txt
```


Running algorithm with different optional parameters:

```shell
./GGA-CGT --P_size 100 --max_gen 500 --p_m 0.83 --p_c 0.2 --k_ncs 1.3 --k_cs 4 --B_size 0.1 --life_span 10 --seed 1 --save_bestSolution 0 --instance instances/BPP1.txt
```


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
