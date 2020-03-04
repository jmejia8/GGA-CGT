/************************************************************************************************************************
GGA-CGT                                                                                                       *
GROUPING GENETIC ALGORITHM WITH CONTROLLED GENE TRANSMISSION FOR THE BIN PACKING PROBLEM                      *
************************************************************************************************************************/
/************************************************************************************************************************
 Author:  Marcela Quiroz-Castellanos                                                          *
      qc.marcela@gmail.com                                                                            *
         Tecnológico Nacional de México                                                                 *
         Instituto Tecnológico de Ciudad Madero                                                         *
         División de Estudios de Posgrado e Investigación                                               *
         Depto. de Sistemas y Computación                                                               *
************************************************************************************************************************/
/************************************************************************************************************************
 The program excecutes GGA-CGT over a set instances using different configurations                          *
 given by the user. Each configuration represents an independent execution of the GA.                       *
                                                                                *
 Reference:                                                                                                   *
    Quiroz-Castellanos, M., Cruz-Reyes, L., Torres-Jimenez, J.,                                     *
    Gómez, C., Huacuja, H. J. F., & Alvim, A. C. (2015).                                          *
      A grouping genetic algorithm with controlled gene transmission for                                *
      the bin packing problem. Computers & Operations Research, 55, 52-64.                                *
                                                                                                          *
 Input:                                                                                                       *
  File "instances.txt" including the name of the BPP instances to be solved;                            *
  Files including the standard instances to be solve;                                           *
  File "configurations.txt" including the parameter values for each experiment;                           *
                                                                                                              *
 Output:                                                                            *
  A set of files "GGA-CGT_(i).txt" including the experimental results for each                            *
  configuration i in the input, stored in directory: Solutions_GGA-CGT;                               *
   If(save_bestSolution = 1) a set of files HGGA_S_(i)_instance.txt including the                         *
   obtained solution for each instance, for each configuration i, stored in directory: Details_GGA-CGT;           *
************************************************************************************************************************/
#include "linked_list.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>


//CONSTANTS DEFINING THE SIZE OF THE PROBLEM
#define ATTRIBUTES          5000
#define P_size_MAX            500

char
  file[300],
  nameC[300];

int
  is_optimal_solution,
   save_bestSolution,
  generation,
   repeated_fitness,
   max_gen,
   life_span,
   P_size,
   seed;

long int
  i,
   j,
   conf,
  number_items,
   bin_capacity,
  best_solution,
   n_,
   L2,
   bin_i,
   higher_weight,
   lighter_weight,
  ordered_weight[ATTRIBUTES],
   weight[ATTRIBUTES],
   permutation[ATTRIBUTES],
   items_auxiliary[ATTRIBUTES],
   ordered_population[P_size_MAX],
   best_individuals[P_size_MAX],
   random_individuals[P_size_MAX];

float
  p_m,
  p_c,
  k_ncs,
  k_cs,
  B_size,
  TotalTime;

long double
  total_accumulated_weight,
   weight1[ATTRIBUTES],
  _p_;

clock_t
  start,
  end;

FILE  *output,
    *input_Configurations,
      *input_Instances;

struct SOLUTION
{
  linked_list L;
  double Bin_Fullness;
};

SOLUTION  global_best_solution[ATTRIBUTES],
      population[P_size_MAX][ATTRIBUTES],
         children[P_size_MAX][ATTRIBUTES];

//Initial seeds for the random number generation
int   seed_emptybin,
    seed_permutation;

//GA COMPONENTS
long int Generate_Initial_Population();
long int Generation();
void  Gene_Level_Crossover_FFD(long int, long int, long int);
void  Adaptive_Mutation_RP(long int, float, int);
void  FF_n_(int); //First Fit with ñ pre-allocated-items (FF-ñ)
void  RP(long int, long int &, long int[], long int); //Rearrangement by Pairs

//BPP Procedures
void  FF(long int, SOLUTION[], long int&, long int, int);
void  LowerBound();

//Auxiliary functions
void  Find_Best_Solution();
void  Sort_Ascending_IndividualsFitness();
void  Sort_Descending_Weights(long int[], long int);
void  Sort_Ascending_BinFullness(long int[], long int);
void  Sort_Descending_BinFullness(long int [], long int);
void  Sort_Random(long int[], long int, int);
void  Copy_Solution(SOLUTION[], SOLUTION[], int);
void  Clean_population();
long int  Used_Items(long int, long int, long int[]);
void  Adjust_Solution(long int);
long int  LoadData();
void  WriteOutput();
void  sendtofile(SOLUTION[]);

//Pseudo-random number generator functions
int get_rand_ij(int *, int, int);
int   get_rand(int *, int);
float randp(int *);
int trand();

char* get_parm(int argc, char *argv[], char* parm){
  for (int i = 0; i < argc-1 ; ++i) {
    if (strcmp(parm, argv[i]) == 0)
    {
      return argv[i+1];
    }
  }

  char* str = (char*) malloc(2*sizeof(char));
  sprintf(str, "%s", "-1");
  return str;
}

int main(int argc, char *argv[])
{
  char  aux[300], nombreC[300], string[300];
   system("mkdir -p Solutions_GGA-CGT");
   system("mkdir -p Details_GGA-CGT");


   // printf(">>> %s\n", get_parm(argc, argv, "--mN"));

   //READING EACH CONFIGURATION IN FILE "configurations.txt", CONTAINING THE PARAMETER VALUES FOR EACH EXPERIMENT
   // if((input_Configurations = fopen("configurations.txt","rt")) == NULL)
   // {  printf("\n INVALID FILE");
   //    //getch();
   //    exit(1);
   // }
   // fscanf(input_Configurations,"%[^\n]",string);
  // while(!feof(input_Configurations))
  // {
    // fscanf(input_Configurations,"%ld",&conf);
      conf = 1;
      strcpy(nameC, "Solutions_GGA-CGT/GGA-CGT_(");
      //itoa(conf, aux, 10);
      sprintf(aux,"%ld", 0);
      strcat(nameC, aux);
      strcat(nameC, ").txt");
      // output = fopen(nameC,"w+");
      char str_tmp[64];
      sprintf(str_tmp, "%s", "--P_size");
      sscanf(get_parm(argc, argv, str_tmp), "%d",    &P_size);
      sprintf(str_tmp, "%s", "--max_gen");
      sscanf(get_parm(argc, argv, str_tmp), "%d",   &max_gen);
      sprintf(str_tmp, "%s", "--p_m");
      sscanf(get_parm(argc, argv, str_tmp), "%f",       &p_m);
      sprintf(str_tmp, "%s", "--p_c");
      sscanf(get_parm(argc, argv, str_tmp), "%f",       &p_c);
      sprintf(str_tmp, "%s", "--k_ncs");
      sscanf(get_parm(argc, argv, str_tmp), "%f",     &k_ncs);
      sprintf(str_tmp, "%s", "--k_cs");
      sscanf(get_parm(argc, argv, str_tmp), "%f",      &k_cs);
      sprintf(str_tmp, "%s", "--B_size");
      sscanf(get_parm(argc, argv, str_tmp), "%f",    &B_size);
      sprintf(str_tmp, "%s", "--life_span");
      sscanf(get_parm(argc, argv, str_tmp), "%d", &life_span);
      sprintf(str_tmp, "%s", "--seed");
      sscanf(get_parm(argc, argv, str_tmp), "%d",      &seed);
      sprintf(str_tmp, "%s", "--save_bestSolution");
      sscanf(get_parm(argc, argv, str_tmp), "%d",      &save_bestSolution);
      sprintf(str_tmp, "%s", "--instance");
      sscanf(get_parm(argc, argv, str_tmp), "%s",      file);
      // sscanf(get_parm(argc, argv, "--save_bestSolution"), "%d",     &save_bestSolution);      

      P_size = P_size >= 0 ? P_size : 100;
      max_gen = max_gen >= 0 ? max_gen : 500;
      p_m = p_m >= 0 ? p_m : 0.83;
      p_c = p_c >= 0 ? p_c : 0.2;
      k_ncs = k_ncs >= 0 ? k_ncs : 1.3;
      k_cs = k_cs >= 0 ? k_cs : 4;
      B_size = B_size >= 0 ? B_size : 0.1;
      life_span = life_span >= 0 ? life_span : 10;
      seed = seed >= 0 ? seed : 1;
      save_bestSolution = save_bestSolution > 0 ? save_bestSolution : 0;

      // ======================================
      // ======================================
      // ======================================
      // printf("P_size %d \n",    P_size);
      // printf("max_gen %d \n",   max_gen);
      // printf("p_m %f \n",       p_m);
      // printf("p_c %f \n",       p_c);
      // printf("k_ncs %f \n",     k_ncs);
      // printf("k_cs %f \n",      k_cs);
      // printf("B_size %f \n",    B_size);
      // printf("life_span %d \n", life_span);
      // printf("seed %d \n",      seed);
      // printf("save_bestSolution %d \n",     save_bestSolution);
      // printf("instance %s \n",     file);
      // ======================================
      // ======================================
      // ======================================

    //   fprintf(output, "CONF\t|P|\tmax_gen\tn_m\tn_c\tk1(non-cloned_solutions)\tk2(cloned_solutions)\t|B|\tlife_span\tseed");
    //   fprintf(output, "\n%ld\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%d\t%d",conf,P_size,max_gen,p_m,p_c,k_ncs,k_cs,B_size,life_span,seed);
    //  fprintf(output, "\nInstancias \t L2 \t Bins \t Gen \t Time");
    // fclose(output);

    //READING FILE "instances.txt" CONTAINING THE NAME OF BPP INSTANCES TO PROCESS
     //  if((input_Instances = fopen("instances.txt","rt")) == NULL)
     // {  printf("\n INVALID FILE");
     //   //getch();
     //    exit(1);
     // }
    // while(!feof(input_Instances))
    // {

      // fscanf(input_Instances,"%s",file);
      LoadData();
      // return 0;
      for(i = 0; i < number_items; i ++)
        ordered_weight[i] = i;
      Sort_Descending_Weights(ordered_weight, number_items);
      LowerBound();
        seed_permutation = seed;
         seed_emptybin = seed;
         for(i = 0; i < P_size; i++)
         {  ordered_population[i] = i;
          random_individuals[i] = i;
            best_individuals[i] = i;
         }
      Clean_population();
      is_optimal_solution = 0;
         generation = 0;
         for(i = 0, j = n_; j < number_items;i++)
          permutation[i] = ordered_weight[j++];
         repeated_fitness = 0;
         //procedure GGA-CGT
         start = clock();
         if(!Generate_Initial_Population())  //Generate_Initial_Population() returns 1 if an optimal solution was found
           for(generation = 0; generation < max_gen; generation++)
        {  if(Generation()) //Generation() returns 1 if an optimal solution was found
              break;
          Find_Best_Solution();
        }
         if(!is_optimal_solution) //is_optimal_solution is 1 if an optimal solution was printed before
      { end = clock();
        TotalTime = (end - start) / (CLOCKS_PER_SEC * 1.0);
            Find_Best_Solution();
            WriteOutput();
         }
    // }
    // fclose(input_Instances);
  // }
  // fclose(input_Configurations);
  // printf("\n\n\tEnd of process");
  // getchar();
  // printf("1\n");
  // printf("2\n");
  // printf("ll666ppk\n");
  printf("%f\n", global_best_solution[number_items + 1].Bin_Fullness);
// [number_items + 1].Bin_Fullness

        return 0;
}



/************************************************************************************************************************
 To generate an initial population P of individuals with FF-ñ packing heuristic.                            *
  population[i][number_items].Bin_Fullness: Saves the fitness of the solution i                                         *
  population[i][number_items + 1].Bin_Fullness:                                                                         *
  Saves the total number of bins in the solution i                                                                      *
  population[i][number_items + 2].Bin_Fullness:                                                                         *
  Saves the generation in which the solution i was generated                                                          *
  population[i][number_items + 3].Bin_Fullness:                                                                         *
  Saves the number of bins in the solution i that are fully at 100%                                                     *
  population[i][number_items + 4].Bin_Fullness:                                                   *
  Saves the fullness of the bin with the highest avaliable capacity in the solution i                                 *
 Output:                                                                                                        *
  (1) when it finds a solution for which the size matches the L2 lower bound                              *
   (0) otherwise                                                                      *
************************************************************************************************************************/
long int Generate_Initial_Population()
{ for(i = 0; i < P_size; i++)
  { FF_n_(i);
      population[i][number_items + 2].Bin_Fullness = generation;
    population[i][number_items].Bin_Fullness /= population[i][number_items + 1].Bin_Fullness;
      if(population[i][number_items + 1].Bin_Fullness == L2)
     {  end = clock();
        Copy_Solution(global_best_solution, population[i], 0);
        global_best_solution[number_items].Bin_Fullness = population[i][number_items].Bin_Fullness;;
      global_best_solution[number_items + 2].Bin_Fullness = generation;
      global_best_solution[number_items + 1].Bin_Fullness = population[i][number_items + 1].Bin_Fullness;
        global_best_solution[number_items + 3].Bin_Fullness = population[i][number_items + 3].Bin_Fullness;
      TotalTime = (end - start) / (CLOCKS_PER_SEC * 1.0);
      WriteOutput();
        is_optimal_solution = 1;
      return (1);
     }
  }
   return 0;
}



/************************************************************************************************************************
 To apply the reproduction technique: Controlled selection and Controlled replacement.                          *
 Output:                                                                                                        *
  (1) when it finds a solution for which the size matches the L2 lower bound                              *
   (2) if more than 0.1*P_size individuals (solutions) have duplicated-fitness                            *
   (0) otherwise                                                                      *
************************************************************************************************************************/
long int Generation()
{
   long int
      f1,
      f2,
         h,
      k;

   /*-----------------------------------------------------------------------------------------------------
   ---------------------------------Controlled selection for crossover------------------------------------
   -----------------------------------------------------------------------------------------------------*/
  Sort_Ascending_IndividualsFitness();
   //if(generation > 1 && repeated_fitness > 0.1*P_size)
  //  return (2);
   Sort_Random(random_individuals,0,(int)(P_size-(int)(P_size*B_size)));
   Sort_Random(best_individuals,(1-p_c)*P_size,P_size);
   k = 0;
   h = P_size - 1;
   for(i = P_size - 1, j = 0; i > P_size - (p_c/2*P_size); i--, j+=2)
    { f1 = best_individuals[h--];
      f2 = random_individuals[k++];
      if(f2 == f1)
    { f1 = best_individuals[h--];
      }
      Gene_Level_Crossover_FFD(ordered_population[f1], ordered_population[f2], j);
       children[j][number_items + 2].Bin_Fullness = generation + 1;
    children[j][number_items].Bin_Fullness /= children[j][number_items+1].Bin_Fullness;
      if(children[j][number_items + 1].Bin_Fullness == L2)
      {  end = clock();
      Copy_Solution(global_best_solution, children[j], 0);
        global_best_solution[number_items].Bin_Fullness = children[j][number_items].Bin_Fullness;;
      global_best_solution[number_items + 2].Bin_Fullness = generation+1;
      global_best_solution[number_items + 1].Bin_Fullness = children[j][number_items + 1].Bin_Fullness;
        global_best_solution[number_items + 3].Bin_Fullness = children[j][number_items + 3].Bin_Fullness;
      TotalTime = (end - start) / (CLOCKS_PER_SEC * 1.0);
      WriteOutput();
        is_optimal_solution = 1;
      return (1);
     }
      Gene_Level_Crossover_FFD(ordered_population[f2], ordered_population[f1], j+1);
      children[j+1][number_items + 2].Bin_Fullness = generation + 1;
    children[j+1][number_items].Bin_Fullness /= children[j+1][number_items+1].Bin_Fullness;
      if(children[j+1][number_items + 1].Bin_Fullness == L2)
      {  end = clock();
      Copy_Solution(global_best_solution, children[j+1], 0);
        global_best_solution[number_items].Bin_Fullness = children[j+1][number_items].Bin_Fullness;;
      global_best_solution[number_items + 2].Bin_Fullness = generation+1;
      global_best_solution[number_items + 1].Bin_Fullness = children[j+1][number_items + 1].Bin_Fullness;
        global_best_solution[number_items + 3].Bin_Fullness = children[j+1][number_items + 3].Bin_Fullness;
      TotalTime = (end - start) / (CLOCKS_PER_SEC * 1.0);
      WriteOutput();
        is_optimal_solution = 1;
        return (1);
     }
   }

   /*-----------------------------------------------------------------------------------------------------
   ---------------------------------Controlled replacement for crossover----------------------------------
   -----------------------------------------------------------------------------------------------------*/
   k = 0;
   for(j = 0; j < p_c/2*P_size - 1; j++)
    Copy_Solution(population[ordered_population[random_individuals[k++]]], children[j], 1);
   k = 0;
   for(i = P_size - 1; i > P_size - (p_c/2*P_size); i--, j++)
   {  while(population[ordered_population[k]][number_items + 2].Bin_Fullness == generation + 1)
        k++;
    Copy_Solution(population[ordered_population[k++]], children[j], 1);
   }
    /*-----------------------------------------------------------------------------------------------------
   --------------------------------Controlled selection for mutation--------------------------------------
   -----------------------------------------------------------------------------------------------------*/
  Sort_Ascending_IndividualsFitness();
   //if(generation > 1 && repeated_fitness > 0.1*P_size)
  //  return (2);
   j = 0;
   for(i = P_size - 1; i > P_size - (p_m*P_size); i--)
  {
      if(j < (int)(P_size*B_size) && generation+ 1 - population[ordered_population[i]][number_items + 2].Bin_Fullness < life_span)
      { /*-----------------------------------------------------------------------------------------------------
      ----------------------------------Controlled replacement for mutation----------------------------------
      -----------------------------------------------------------------------------------------------------*/
        Copy_Solution(population[ordered_population[j]], population[ordered_population[i]], 0);
      Adaptive_Mutation_RP(ordered_population[j], k_cs, 1);
        population[ordered_population[j]][number_items + 2].Bin_Fullness = generation + 1;
      population[ordered_population[j]][number_items].Bin_Fullness /= population[ordered_population[j]][number_items+1].Bin_Fullness;
        if(population[ordered_population[j]][number_items + 1].Bin_Fullness == L2)
      {  end = clock();
        Copy_Solution(global_best_solution, population[ordered_population[j]], 0);
          global_best_solution[number_items].Bin_Fullness = population[ordered_population[j]][number_items].Bin_Fullness;;
        global_best_solution[number_items + 2].Bin_Fullness = generation+1;
        global_best_solution[number_items + 1].Bin_Fullness = population[ordered_population[j]][number_items + 1].Bin_Fullness;
          global_best_solution[number_items + 3].Bin_Fullness = population[ordered_population[j]][number_items + 3].Bin_Fullness;
        TotalTime = (end - start) / (CLOCKS_PER_SEC * 1.0);
        WriteOutput();
          is_optimal_solution = 1;
        return (1);
       }
         j++;
    }
    else
    { Adaptive_Mutation_RP(ordered_population[i], k_ncs, 0);
      population[ordered_population[i]][number_items + 2].Bin_Fullness = generation + 1;
      population[ordered_population[i]][number_items].Bin_Fullness /= population[ordered_population[i]][number_items+1].Bin_Fullness;
       if(population[ordered_population[i]][number_items+1].Bin_Fullness == L2)
       {  end = clock();
          Copy_Solution(global_best_solution, population[ordered_population[i]], 0);
          global_best_solution[number_items].Bin_Fullness = population[ordered_population[i]][number_items].Bin_Fullness;;
        global_best_solution[number_items + 2].Bin_Fullness = generation+1;
        global_best_solution[number_items + 1].Bin_Fullness = population[ordered_population[i]][number_items + 1].Bin_Fullness;
        global_best_solution[number_items + 3].Bin_Fullness = population[ordered_population[i]][number_items + 3].Bin_Fullness;
        TotalTime = (end - start) / (CLOCKS_PER_SEC * 1.0);
        WriteOutput();
          is_optimal_solution = 1;
        return (1);
       }
    }
   }
   return 0;
}



/************************************************************************************************************************
 To recombine two parent solutions producing a child solution.                                          *
 Input:                                                                                                         *
  The positions in the population of the two parent solutions: father_1 and father_2                        *
  The position in the set of children of the child solution: child                                    *
************************************************************************************************************************/
void Gene_Level_Crossover_FFD(long int father_1, long int father_2, long int child)
{
  long int  k,
      counter,
      k2 = 0,
      ban = 1,
         items[ATTRIBUTES] = {0},
         free_items[ATTRIBUTES] = {0};
   children[child][number_items + 4].Bin_Fullness = bin_capacity;

   if(population[father_1][number_items + 1].Bin_Fullness > population[father_2][number_items + 1].Bin_Fullness)
    counter = population[father_1][number_items + 1].Bin_Fullness;
   else
    counter = population[father_2][number_items + 1].Bin_Fullness;

   long int *random_order1 = new long int [counter];
   long int *random_order2 = new long int [counter];


   for(k = 0; k < counter; k++)
   {  random_order1[k] = k;
     random_order2[k] = k;
   }

   Sort_Random(random_order1,0, population[father_1][number_items + 1].Bin_Fullness);
    Sort_Random(random_order2,0, population[father_2][number_items + 1].Bin_Fullness);
   Sort_Descending_BinFullness(random_order1, father_1);
    Sort_Descending_BinFullness(random_order2, father_2);

   for(k = 0; k < population[father_1][number_items + 1].Bin_Fullness; k++)
   {
    if(population[father_1][random_order1[k]].Bin_Fullness >= population[father_2][random_order2[k]].Bin_Fullness)
      { ban = Used_Items(father_1, random_order1[k], items);
      if (ban == 1)
      {
          children[child][k2].L.clone_linked_list(population[father_1][random_order1[k]].L);
        children[child][k2++].Bin_Fullness = population[father_1][random_order1[k]].Bin_Fullness;
            if(children[child][k2-1].Bin_Fullness < children[child][number_items + 4].Bin_Fullness)
              children[child][number_items + 4].Bin_Fullness = children[child][k2-1].Bin_Fullness;
      }
         if(population[father_2][random_order2[k]].Bin_Fullness > 0)
         {
          ban = Used_Items(father_2, random_order2[k], items);
        if (ban == 1)
          {
            children[child][k2].L.clone_linked_list(population[father_2][random_order2[k]].L);
          children[child][k2++].Bin_Fullness = population[father_2][random_order2[k]].Bin_Fullness;
               if(children[child][k2-1].Bin_Fullness < children[child][number_items + 4].Bin_Fullness)
                children[child][number_items + 4].Bin_Fullness = children[child][k2-1].Bin_Fullness;
         }
         }

      }
      else
    {
        if(population[father_2][random_order2[k]].Bin_Fullness > 0)
        { ban = Used_Items(father_2, random_order2[k], items);
        if (ban == 1)
          {
            children[child][k2].L.clone_linked_list(population[father_2][random_order2[k]].L);
          children[child][k2++].Bin_Fullness = population[father_2][random_order2[k]].Bin_Fullness;
               if(children[child][k2-1].Bin_Fullness < children[child][number_items + 4].Bin_Fullness)
                children[child][number_items + 4].Bin_Fullness = children[child][k2-1].Bin_Fullness;
          }
         }
         ban = Used_Items(father_1, random_order1[k], items);
      if (ban == 1)
      {
          children[child][k2].L.clone_linked_list(population[father_1][random_order1[k]].L);
        children[child][k2++].Bin_Fullness = population[father_1][random_order1[k]].Bin_Fullness;
            if(children[child][k2-1].Bin_Fullness < children[child][number_items + 4].Bin_Fullness)
              children[child][number_items + 4].Bin_Fullness = children[child][k2-1].Bin_Fullness;
      }
      }
   }
  k = 0;
   for(counter = 0; counter < number_items; counter++)
  {
      if(items[ordered_weight[counter]] == 0)
        free_items [k++] = ordered_weight[counter];
   }
   if(k > 0)
   {  bin_i = 0;
      for(counter = 0; counter < k-1; counter++)
      FF(free_items[counter], children[child], k2, bin_i,0);
    FF(free_items[counter], children[child], k2, bin_i,1);
   }
   else
    for(k = 0; k < k2; k++)
        children[child][number_items].Bin_Fullness += pow((children[child][k].Bin_Fullness / bin_capacity), 2);
   children[child][number_items+1].Bin_Fullness = k2;
   free(random_order1);
  free(random_order2);
}



/************************************************************************************************************************
 To produce a small modification in a solution.                                                 *
 Input:                                                                                                         *
  The position in the population of the solution to mutate: individual                                *
  The rate of change to calculate the number of bins to eliminate: k                                  *
  A value that indicates if the solution was cloned: is_cloned                                      *
************************************************************************************************************************/
void Adaptive_Mutation_RP(long int individual, float k, int is_cloned)
{
  long int
         number_bins,
         i,
         lightest_bin = 0,
         number_free_items = 0,
         free_items[ATTRIBUTES] = {0},
         ordered_BinFullness[ATTRIBUTES] = {0};
   node *p;
   for(i = 0; i < population[individual][number_items + 1].Bin_Fullness; i++)
    ordered_BinFullness[i] = i;
   if(is_cloned)
      Sort_Random(ordered_BinFullness,0, population[individual][number_items + 1].Bin_Fullness);
   Sort_Ascending_BinFullness(ordered_BinFullness, individual);
   i = 1;
   while(population[individual][ordered_BinFullness[i]].Bin_Fullness < bin_capacity && i < population[individual][number_items + 1].Bin_Fullness)
    i++;
   _p_ = 1 / (float)(k);
   number_bins = (long int)ceil(i*((2 - i/population[individual][number_items + 1].Bin_Fullness) / pow(i,_p_))*(1 - ((double)get_rand(&seed_emptybin,(long int)ceil((1/pow(i,_p_))*100))/100)));
  for(i = 0; i < number_bins; i++)
  {  p = population[individual][ordered_BinFullness[lightest_bin]].L.first;
    while(p != NULL)
      { free_items[number_free_items++] = p->data;
        p = p->next;
      }
    population[individual][ordered_BinFullness[lightest_bin]].L.free_linked_list();
    population[individual][ordered_BinFullness[lightest_bin]].Bin_Fullness = 0;
      lightest_bin++;
  }
  population[individual][number_items + 1].Bin_Fullness -= number_bins;
  number_bins = population[individual][number_items + 1].Bin_Fullness;
  Adjust_Solution(individual);
   RP(individual, number_bins, free_items, number_free_items);
  population[individual][number_items+1].Bin_Fullness = number_bins;
}



/************************************************************************************************************************
 To generate a random BPP solution with the ñ large items packed in separate bins.                          *
 Input:                                                                                                         *
  The position in the population of the new solution: individual                                    *
************************************************************************************************************************/
void FF_n_(int individual)
{
  long int
      i,
      j = 0,
      total_bins = 0;
   bin_i = 0;
   population[individual][number_items + 3].Bin_Fullness = 0;
   population[individual][number_items + 4].Bin_Fullness = bin_capacity;
   if(n_ > 0)
  { for(i = 0; i < n_; i++)
    {  population[individual][i].Bin_Fullness = weight[ordered_weight[i]];
      population[individual][i].L.insert(ordered_weight[i]);
        total_bins++;
       if(population[individual][i].Bin_Fullness < population[individual][number_items + 4].Bin_Fullness)
          population[individual][number_items + 4].Bin_Fullness = population[individual][i].Bin_Fullness;
    }
     i = number_items - i;
    Sort_Random(permutation,0, i);
    for(j = 0; j < i-1; j++)
        FF(permutation[j], population[individual], total_bins, bin_i,0);
     FF(permutation[j], population[individual], total_bins, bin_i,1);
   }
   else
   {  Sort_Random(permutation,0, number_items);
    for(j = 0; j < number_items-1; j++)
        FF(permutation[j], population[individual], total_bins, bin_i,0);
    FF(permutation[j], population[individual], total_bins, bin_i,1);
   }
  population[individual][number_items + 1].Bin_Fullness = total_bins;
}



/************************************************************************************************************************
 To reinsert free items into an incomplete BPP solution.                                            *
 Input:                                                                                                                 *
  The position in the population of the incomplete solution where the free items must be reinserted: individual     *
   The number of bins of the partial_solution: b                                                *
   A set of free items to be reinserted into the partial_solution: F                                    *
   The number of free items of F: number_free_items                                             *
************************************************************************************************************************/
void RP(long int individual, long int &b, long int F[], long int number_free_items)
{
  long int
      i,
         k,
         k2,
         ban,
         sum = 0,
         total_free = 0,
         ordered_BinFullness[ATTRIBUTES] = {0},
      *new_free_items = new long int [2];

   node   *p,
      *s,
         *aux;

   higher_weight = weight[F[0]];
   lighter_weight = weight[F[0]];
   bin_i = b;
   population[individual][number_items].Bin_Fullness = 0;
   population[individual][number_items + 3].Bin_Fullness = 0;
   population[individual][number_items + 4].Bin_Fullness = bin_capacity;

   for (i = 0; i < b; i++ )
    ordered_BinFullness[i] = i;
   Sort_Random(ordered_BinFullness,0, b);
  Sort_Random(F,0, number_free_items);

  for (i = 0; i < b; i++)
  {  sum = (long int)population[individual][ordered_BinFullness[i]].Bin_Fullness;
    p = population[individual][ordered_BinFullness[i]].L.first;
    while(p->next != NULL)
    {  ban = 0;
         aux = p;
         s = p->next;
        while(s != NULL)
      {  for (k = 0; k < number_free_items - 1; k++)
        {  if(i == b-1)
                if(weight[F[k]] > higher_weight)
                    higher_weight = weight[F[k]];
              for (k2 = k + 1; k2 < number_free_items; k2++)
               {  if(weight[F[k]] >= weight[p->data] + weight[s->data] && ((sum - (weight[p->data] + weight[s->data]) + (weight[F[k]]) <= bin_capacity)))
                  {  sum = sum - (weight[p->data] + weight[s->data]) + (weight[F[k]]);
                    new_free_items[0] = p->data;
                    new_free_items[1] = s->data;
                     p->data = F[k];
                     aux->next = s->next;
                     free(s);
                     if(population[individual][ordered_BinFullness[i]].L.last == s)
                      population[individual][ordered_BinFullness[i]].L.last = aux;
                     population[individual][ordered_BinFullness[i]].L.num--;
                     F[k] = new_free_items[0];
                     F[number_free_items + total_free] = new_free_items[1];
                     total_free++;
                     ban = 1;
                     break;
                  }
                  if(weight[F[k2]] >= weight[p->data] + weight[s->data] && ((sum - (weight[p->data] + weight[s->data]) + (weight[F[k2]]) <= bin_capacity)))
                  {  sum = sum - (weight[p->data] + weight[s->data]) + (weight[F[k2]]);
                    new_free_items[0] = p->data;
                    new_free_items[1] = s->data;
                     p->data = F[k2];
                     aux->next = s->next;
                     free(s);
                     if(population[individual][ordered_BinFullness[i]].L.last == s)
                      population[individual][ordered_BinFullness[i]].L.last = aux;
                     population[individual][ordered_BinFullness[i]].L.num--;
                     F[k2] = new_free_items[0];
                     F[number_free_items + total_free] = new_free_items[1];
                     total_free++;
                     ban = 1;
                     break;
                  }
                  if((weight[F[k]] + weight[F[k2]] > weight[p->data] + weight[s->data]) || ((weight[F[k]] + weight[F[k2]] == weight[p->data] + weight[s->data]) && !(weight[F[k]] == weight[p->data] || weight[F[k]] == weight[s->data])))
                  {  if(sum - (weight[p->data] + weight[s->data]) + (weight[F[k]] + weight[F[k2]]) > bin_capacity)
                        break;
                    sum = sum - (weight[p->data] + weight[s->data]) + (weight[F[k]] + weight[F[k2]]);
                    new_free_items[0] = p->data;
                    new_free_items[1] = s->data;
                     p->data = F[k];
                     s->data = F[k2];
                     F[k] = new_free_items[0];
                     F[k2] = new_free_items[1];
                     if(sum == bin_capacity)
                     {  ban = 1;
                        break;
                     }
                  }
              }
               if(ban)
                break;
        }
            if(ban)
                break;
            aux = s;
            s = s->next;
      }
         if(ban)
            break;
         p = p->next;
      }
      population[individual][ordered_BinFullness[i]].Bin_Fullness = sum;
      if(population[individual][ordered_BinFullness[i]].Bin_Fullness < population[individual][number_items + 4].Bin_Fullness)
        population[individual][number_items + 4].Bin_Fullness = population[individual][ordered_BinFullness[i]].Bin_Fullness;
      if(population[individual][ordered_BinFullness[i]].Bin_Fullness == bin_capacity)
        population[individual][number_items + 3].Bin_Fullness++;
      else if(population[individual][ordered_BinFullness[i]].Bin_Fullness + weight[ordered_weight[number_items-1]] <= bin_capacity)
      {  if(ordered_BinFullness[i] < bin_i)
          bin_i = ordered_BinFullness[i];
      }
  }
   for(i = 0; i < bin_i; i++)
    population[individual][number_items].Bin_Fullness += pow((population[individual][i].Bin_Fullness / bin_capacity), 2);

   free(new_free_items);
   number_free_items += total_free;

   if(higher_weight < .5*bin_capacity)
    Sort_Random(F,0, number_free_items);
   else
   {  Sort_Descending_Weights(F, number_free_items);
      lighter_weight = weight[F[number_free_items-1]];
   }

   if(lighter_weight > bin_capacity - population[individual][number_items + 4].Bin_Fullness)
   {  for(i = bin_i; i < b; i++)
        population[individual][number_items].Bin_Fullness += pow((population[individual][i].Bin_Fullness / bin_capacity), 2);
    bin_i = b;
   }
   for(i = 0; i < number_free_items-1; i++)
    FF(F[i], population[individual], b, bin_i, 0);
   FF(F[i], population[individual], b, bin_i, 1);
}



/************************************************************************************************************************
 To insert an item into an incomplete BPP solution.                                             *
 Input:                                                                                                         *
   An item to be inserted into the individual: item                                             *
  An incomplete chromosome where the item must be inserted: individual                                      *
   The number of bins of the individual: total_bins                                             *
   The first bin that could have sufficient available capacity to store the item: beginning                   *
   A value that indicates if it is the last item to be stored into the individual: is_last                    *
************************************************************************************************************************/
void FF(long int item, SOLUTION individual[], long int &total_bins, long int beginning, int is_last)
{
  long int  i;

  if(!is_last && weight[item] > (bin_capacity - (long int)individual[number_items + 4].Bin_Fullness))
      i = total_bins;
  else
    for(i = beginning; i < total_bins; i++)
    {
      if((long int)individual[i].Bin_Fullness + weight[item] <= bin_capacity)
      {
        individual[i].Bin_Fullness += weight[item];
        individual[i].L.insert(item);
            if((long int)individual[i].Bin_Fullness == bin_capacity)
              individual[number_items + 3].Bin_Fullness++;
           if(is_last)
          { for(i=i; i < total_bins; i++)
              individual[number_items].Bin_Fullness += pow((individual[i].Bin_Fullness / bin_capacity), 2);
            return;
          }
           if((long int)individual[i].Bin_Fullness + weight[ordered_weight[number_items-1]] > bin_capacity && i == bin_i)
           {  bin_i++;
            individual[number_items].Bin_Fullness += pow((individual[i].Bin_Fullness / bin_capacity), 2);
           }
        return;
      }
        if(is_last)
          individual[number_items].Bin_Fullness += pow((individual[i].Bin_Fullness / bin_capacity), 2);
    }
  individual[i].Bin_Fullness += weight[item];
  individual[i].L.insert(item);
   if(individual[i].Bin_Fullness < individual[number_items + 4].Bin_Fullness)
      individual[number_items + 4].Bin_Fullness = individual[i].Bin_Fullness;
   if(is_last)
      individual[number_items].Bin_Fullness += pow((individual[i].Bin_Fullness / bin_capacity), 2);
  total_bins++;
}



/************************************************************************************************************************
 To calculate the lower bound L2 of Martello and Toth and the ñ large items n_                            *
************************************************************************************************************************/
void LowerBound()
{
  long int k, m, i, j, aux1, aux2;
   long double sjx=0, sj2=0, sj3=0;
   long int jx=0, cj12, jp=0, jpp=0, cj2;

   while(weight[ordered_weight[jx]] > bin_capacity/2 && jx < number_items)
    jx++;
   n_ = jx;
   if(jx == number_items)
   {  L2 = jx;
    return;
   }
   if(jx == 0)
  {  if(fmod(total_accumulated_weight,bin_capacity) >= 1)
      L2 = (long int)ceil(total_accumulated_weight / bin_capacity);
    else
      L2 = (long int)(total_accumulated_weight / bin_capacity);
      return;
   }
   else
   {  cj12 = jx;
      for(i=jx; i < number_items; i++)
      sjx += weight[ordered_weight[i]];
      jp = jx;
      for(i = 0; i < jx; i++)
      { if(weight[ordered_weight[i]] <= bin_capacity - weight[ordered_weight[jx]])
         {  jp = i;
          break;
         }
      }

      cj2 = jx - jp;
      for(i=jp; i <= jx-1; i++)
        sj2 += weight[ordered_weight[i]];
      jpp = jx;
      sj3 = weight[ordered_weight[jpp]];
      ordered_weight[number_items] = number_items;
      weight[number_items]=0;
      while(weight[ordered_weight[jpp+1]]==weight[ordered_weight[jpp]])
      {
        jpp++;
        sj3 += weight[ordered_weight[jpp]];
      }
      L2 = cj12;

      do
      {  if(fmod((sj3 + sj2),bin_capacity) >= 1)
          aux1 = (long int)ceil((sj3 + sj2)/bin_capacity - cj2);
        else
          aux1 = (long int)((sj3 + sj2)/bin_capacity - cj2);

        if(L2 < (cj12 + aux1))
          L2 = cj12 + aux1;
        jpp++;
         if(jpp < number_items)
         {  sj3 += weight[ordered_weight[jpp]];
          while(weight[ordered_weight[jpp+1]] == weight[ordered_weight[jpp]])
            {
              jpp++;
              sj3 += weight[ordered_weight[jpp]];
            }
            while(jp > 0 && weight[ordered_weight[jp-1]] <= bin_capacity - weight[ordered_weight[jpp]])
            { jp--;
              cj2++;
               sj2 += weight[ordered_weight[jp]];
            }
         }
         if(fmod((sjx + sj2),bin_capacity) >= 1)
          aux2 = (long int)ceil((sjx + sj2) / bin_capacity - cj2 );
         else
          aux2 = (long int)((sjx + sj2) / bin_capacity - cj2 );
      }while(jpp <= number_items || (cj12 + aux2) > L2);
   }
}



/************************************************************************************************************************
 To find the solution with the highest fitness of the population and update the global_best_solution              *
************************************************************************************************************************/
void Find_Best_Solution()
{
  long int  i,
      best_individual = 0;
  for(i = 0; i < P_size; i++)
  { if(population[i][number_items].Bin_Fullness > population[best_individual][number_items].Bin_Fullness)
      best_individual = i;
  }
  if(generation + 1 > 1)
  {
    if(population[best_individual][number_items].Bin_Fullness > global_best_solution[number_items].Bin_Fullness)
    {
      Copy_Solution(global_best_solution, population[best_individual], 0);
    }
  }
  else
  {
      Copy_Solution(global_best_solution, population[best_individual], 0);
  }
}



/************************************************************************************************************************
 To sort the individuals of the population in ascending order of their fitness                              *
************************************************************************************************************************/
void Sort_Ascending_IndividualsFitness()
{
  long int i,
       k = P_size - 1,
       i2 = 0,
     aux,
       ban = 1;

  while(ban)
  {
    ban = 0;
    for(i = i2; i < k; i++)
    {
         if(population[ordered_population[i]][number_items].Bin_Fullness > population[ordered_population[i+1]][number_items].Bin_Fullness)
      {
        aux = ordered_population[i];
        ordered_population[i] = ordered_population[i+1];
        ordered_population[i+1] = aux;
            ban = 1;
      }
         else if(population[ordered_population[i]][number_items].Bin_Fullness == population[ordered_population[i+1]][number_items].Bin_Fullness)
         {
        aux = ordered_population[i+1];
        ordered_population[i+1] = ordered_population[i2];
        ordered_population[i2] = aux;
            i2++;
         }
    }
    k--;
  }
   repeated_fitness = i2;
}



/************************************************************************************************************************
 To sort the bins of a solution in ascending order of their filling                                   *
 Input:                                                                                                                 *
  An array to save the order of the bins: ordered_BinFullness                                       *                                                                                                                 *
  The position in the population of the solution: individual                                        *
************************************************************************************************************************/
void Sort_Ascending_BinFullness(long int ordered_BinFullness[], long int individual)
{
  long int  m,
      k,
         temporary_variable,
         ban = 1;

  k = population[individual][number_items + 1].Bin_Fullness - 1;
  while(ban)
  {
    ban = 0;
    for(m = 0; m < k; m++)
    {
      if(population[individual][ordered_BinFullness[m]].Bin_Fullness > population[individual][ordered_BinFullness[m+1]].Bin_Fullness)
      {
        temporary_variable = ordered_BinFullness[m];
        ordered_BinFullness[m] = ordered_BinFullness[m + 1];
        ordered_BinFullness[m + 1] = temporary_variable;
        ban = 1;
      }
    }
    k--;
  }
}



/************************************************************************************************************************
 To sort the bins of a solution in descending order of their filling                                    *
 Input:                                                                                                                 *
  An array to save the order of the bins: ordered_BinFullness                                       *                                                                                                                 *
  The position in the population of the solution: individual                                        *
************************************************************************************************************************/
void Sort_Descending_BinFullness(long int ordered_BinFullness[], long int individual)
{
  long int  m,
      k,
         temporary_variable,
         ban = 1;

  k = population[individual][number_items + 1].Bin_Fullness - 1;
  while(ban)
  {
    ban = 0;
    for(m = 0; m < k; m++)
    {
      if(population[individual][ordered_BinFullness[m]].Bin_Fullness < population[individual][ordered_BinFullness[m+1]].Bin_Fullness)
      {
        temporary_variable = ordered_BinFullness[m];
        ordered_BinFullness[m] = ordered_BinFullness[m + 1];
        ordered_BinFullness[m + 1] = temporary_variable;
        ban = 1;
      }
    }
    k--;
  }
}



/************************************************************************************************************************
 To sort the elements between the positions [k] and [n] of an array in random order                         *
 Input:                                                                                                                 *
  The array to be randomized: random_array                                                    *                                                                                                                 *
  The initial random position: k                                                          *
    The final random position: n                                                            *
************************************************************************************************************************/
void Sort_Random(long int random_array[], long int k, int n)
{
   long int i,
        aux,
        random_number;

   for(i = n - 1; i >= k; i--)
   {
      random_number = k + get_rand(&seed_permutation,n-k) - 1;
      aux = random_array[random_number];
      random_array[random_number] = random_array[i];
      random_array[i] = aux;
      if(weight[random_array[i]] < lighter_weight)
        lighter_weight = weight[random_array[i]];
    if(weight[random_array[random_number]] < lighter_weight)
        lighter_weight = weight[random_array[random_number]];
   }

}



/************************************************************************************************************************
 To sort a set of items in descending order of their weights                                        *
 Input:                                                                                                                 *
  An array to save the order of the items: ordered_weight                                         *                                                                                                                 *
  The number of items in the set: n                                                       *
************************************************************************************************************************/
void Sort_Descending_Weights(long int ordered_weight[], long int n)
{
  long int  m,
      k,
         temporary_variable,
         ban = 1;

  k = n - 1;
  while(ban)
  {
    ban = 0;
    for(m = 0; m < k; m++)
    {
      if(weight[ordered_weight[m]] < weight[ordered_weight[m+1]])
      {
        temporary_variable = ordered_weight[m];
        ordered_weight[m] = ordered_weight[m + 1];
        ordered_weight[m + 1] = temporary_variable;
        ban = 1;
      }
    }
    k--;
  }
}



/************************************************************************************************************************
 To copy solution2 in solution                                                            *
 Input:                                                                                                                 *
  A solution to save the copied solution: solution                                              *                                                                                                                 *
  The solution to be copied: solution2                                                        *
  A value that indicates if the copied solution must be deleted: delete_solution2                         *
************************************************************************************************************************/
void Copy_Solution(SOLUTION solution[], SOLUTION solution2[], int delete_solution2)
{
  long int  j;

   for(j = 0; j < solution2[number_items + 1].Bin_Fullness; j++)
  {
    solution[j].Bin_Fullness = solution2[j].Bin_Fullness;
    solution[j].L.clone_linked_list(solution2[j].L);
      if(delete_solution2)
      { solution2[j].Bin_Fullness = 0;
      solution2[j].L.free_linked_list();
      }
  }
   while(j < solution[number_items + 1].Bin_Fullness)
   {  solution[j].Bin_Fullness = 0;
    solution[j++].L.free_linked_list();
   }
   solution[number_items].Bin_Fullness = solution2[number_items].Bin_Fullness;
   solution[number_items + 1].Bin_Fullness = solution2[number_items + 1].Bin_Fullness;
   solution[number_items + 2].Bin_Fullness = solution2[number_items + 2].Bin_Fullness;
   solution[number_items + 3].Bin_Fullness = solution2[number_items + 3].Bin_Fullness;
   solution[number_items + 4].Bin_Fullness = solution2[number_items + 4].Bin_Fullness;
   if(delete_solution2)
   {  solution2[number_items].Bin_Fullness = 0;
    solution2[number_items + 1].Bin_Fullness = 0;
    solution2[number_items + 2].Bin_Fullness = 0;
    solution2[number_items + 3].Bin_Fullness = 0;
    solution2[number_items + 4].Bin_Fullness = 0;
   }

}



/************************************************************************************************************************
 To free the memory of the individuals of the population                                            *
************************************************************************************************************************/
void Clean_population()
{
  long int  i,
      j;

  for(i = 0; i < P_size; i++)
  {
    for(j = 0; j < number_items + 5; j++)
    {
      population[i][j].L.free_linked_list();
      population[i][j].Bin_Fullness = 0;
         children[i][j].L.free_linked_list();
      children[i][j].Bin_Fullness = 0;
    }
   }
}



/************************************************************************************************************************
 To check if any of the items of the current bin is already in the solution                             *
 Input:                                                                                                                 *
  The position in the population of the solution: individual                                        *
  A new bin that could be added to the solution: bin                                              *
   An array that indicates the items that are already in the solution: items                                            *
 Output:                                                                                                                *
  (1) when none of the items in the current bin is already in the solution                                *
   (0) otherwise                                                                      *
************************************************************************************************************************/
long int Used_Items(long int individual, long int bin, long int items[])
{
  long int   item,
      i,
      counter = 0;
   node *p;

  p = population[individual][bin].L.first;
    while(p != NULL)
    { item = p->data;
      p = p->next;
      if(items [item] != 1)
      {
      items_auxiliary[counter++] = item;
      items[item] = 1;
    }
      else
      {
      for(i = 0; i < counter ; i++)
            items[items_auxiliary[i]] = 0;
        return 0;
      }
    }
    return (1);
}



/************************************************************************************************************************
 To put together all the used bins of the solution                                                *
 Input:                                                                                                                 *
  The position in the population of the solution: individual                                        *
************************************************************************************************************************/
void Adjust_Solution(long int individual)
{
  long int  i = 0,
      j = 0,
         k;
   while(population[individual][i].Bin_Fullness > 0)
    i++;
   for(j = i, k = i; j < number_items; j++, k++)
  { if(j < population[individual][number_items + 1].Bin_Fullness)
    { while(population[individual][k].Bin_Fullness == 0)
        k++;
        population[individual][j].L.first = NULL;
         population[individual][j].L.last = NULL;
      population[individual][j].Bin_Fullness = population[individual][k].Bin_Fullness;
        population[individual][j].L.get_linked_list(population[individual][k].L);

      }
      else
      { population[individual][j].Bin_Fullness = 0;
        population[individual][j].L.first = NULL;
         population[individual][j].L.last = NULL;
         population[individual][j].L.num = 0;
      }
   }
}



/************************************************************************************************************************
 To read the data defining a BPP instance                                                     *
************************************************************************************************************************/
long int LoadData()
{
    char  string[300];
  long  k;
   long int ban = 0;
   long double bin_capacity1;
   long double total_accumulated_aux = 0;

    FILE  *data_file;

   string[0] = '\0';
    strcpy(string, file);
  if((data_file = fopen(string, "rt")) == NULL)
    {
      printf("\nThere is no data file ==> [%s]%c", string, 7);
      return 0;
    }
   // printf("\nThe file is %s", string);
    fgets(string, 300, data_file);
    fscanf(data_file, "%ld\n", &number_items);
    bin_capacity = 0;
    fgets(string, 300, data_file);
    fscanf(data_file, "%Lf\n", &bin_capacity1);
    best_solution = 0;
    fgets(string, 300, data_file);
    fscanf(data_file, "%ld\n", &best_solution);
    fgets(string, 300, data_file);
    total_accumulated_weight = 0;
    for(k = 0; k < number_items; k++)
    {   fscanf(data_file, "%Lf", &weight1[k]);
      weight[k] = (long int)weight1[k];
      total_accumulated_weight = (total_accumulated_weight + weight[k]);
      total_accumulated_aux += weight1[k];
      if(ban == 0)
    {
        if(weight1[k] / weight[k] > 1)
          ban = 1;
      }
    }
   if(ban)
   {  total_accumulated_weight = 0;
    for(k = 0; k < number_items; k++)
    { weight[k] = (long int)(ceil(weight1[k]*bin_capacity1 - .5) );
         total_accumulated_weight = (total_accumulated_weight + weight[k]);
      }
      bin_capacity1 *= bin_capacity1;
   }
   bin_capacity = (long int)bin_capacity1;
   fclose(data_file);
   if(ban)
   {  if((long int)total_accumulated_weight != (long int)(ceil(total_accumulated_aux*sqrt(bin_capacity) - .5) ))
    { printf("\t Error loading weights");
      //getch();
        exit(1);
    }
   }
   return 1;
}



/************************************************************************************************************************
 To print the performance of the procedure on a BPP instance in a data file                             *
************************************************************************************************************************/
void WriteOutput()
{

  //  output = fopen(nameC, "a");
  //  fprintf(output, "\n%s \t %d \t %d \t %d \t %f", file, (int)L2, (int)global_best_solution[number_items + 1].Bin_Fullness, generation, TotalTime);
  // if(save_bestSolution == 1)
  //   sendtofile(global_best_solution);
  //  fclose(output);
}



/************************************************************************************************************************
 To print the global best solution in a data file                                               *
************************************************************************************************************************/
void sendtofile(SOLUTION best[])
{
  char  string1[300],
      fil[300],
         aux[300];

  long double accumulated = 0;
  long int   bin,
    ban = 1,
    item = 0,
    position = 0,
      bins[ATTRIBUTES] = {0},
    n_bins = best[number_items + 1].Bin_Fullness;

   int  binError = -1,
         banError = 0;
   long int j;



    FILE *output;
   node *p;
    strcpy(fil, "Details_GGA-CGT/GGA-CGT_S_(");
    strcpy(string1, file);
   //itoa(conf, aux, 10);
   sprintf(aux,"%ld", conf);
   strcat(fil,aux);
   strcat(fil,")_");
  strcat(fil,string1);
  if((output = fopen(fil, "w+")) == NULL)
    {
      //printf("\nThere is no data file ==> [%s]%c", file, 7);
    getchar();
    exit(1);
    }
  fprintf(output,"Instance:\t%s\n", file);
  fprintf(output,"Number of items:\t%ld\n", number_items);
  fprintf(output,"Bin capacity:\t%ld\n", bin_capacity);
  fprintf(output,"L2:\t%ld\n", L2);
   fprintf(output,"\n****************************GGA-CGT global best solution******************************\n");
  fprintf(output,"Number of bins:\n%ld\n", n_bins);
    fprintf(output,"Fitness:\n%f\n", best[number_items].Bin_Fullness);
  fprintf(output,"Optimal order of the weights:\n");

  for (bin = 0; bin < n_bins; bin++)
  {  bins[bin] = 0;
    p = best[bin].L.first;
    while(true)
    {  if(p == NULL)
          break;
         item = p->data;
         p = p->next;
      bins[bin] += weight[item];
         accumulated += weight[item];
         fprintf(output, "%ld\n",weight[item]);
      if(bins[bin] > bin_capacity)
      {
        printf("ERROR the capacity of bin %ld was exceeded", bin);
            binError = bin;
        getchar();
            banError = 1;
      }

    }
  }

   if(accumulated != total_accumulated_weight)
   {   printf("ERROR inconsistent sum of weights");
       getchar();

   }

  fprintf(output,"\nDetailed solution:");
  for (j=0; j < n_bins; j++)
  {
    if((long int)bins[j] > (long int)bin_capacity)
        fprintf(output, " \n ********************ERROR the capacity of the bin was exceeded******************");
    fprintf(output, "\n\nBIN %ld\nFullness: %ld Gap: %ld\nStored items:\t ",j + 1, bins[j], bin_capacity-bins[j]);
      p = best[j].L.first;
    for (position=0; ; position++)
    {  if(p == NULL)
          break;
      item = p->data;
      p = p->next;
      fprintf(output, "[Item: %ld, Weight: %ld]\t", item + 1,weight[item]);
    }

  }

  fclose(output);

   if(banError)
    exit(1);
}



/******************************************************
 Author: Adriana Cesario de Faria Alvim               *
         (alvim@inf.puc-rio.br, adriana@pep.ufrj.br ) *
*******************************************************/
/***************************************************************************
 Portable pseudo-random number generator                                   *
 Machine independent as long as the machine can represent all the integers *
         in the interval [- 2**31 + 1, 2**31 - 1].                         *
                                                                           *
 Reference: L. Schrage, "A more Portable Fortran Random Number Generator", *
            ACM Transactions on Mathematical Software, Vol. 2, No. 2,      *
            (June, 1979).                                                  *
                                                                           *
 The generator produces a sequence of positive integers, "ix",             *
      by the recursion: ix(i+1) = A * ix(i) mod P, where "P" is Mersenne   *
      prime number (2**31)-1 = 2147483647 and "A" = 7**5 = 16807. Thus all *
      integers "ix" produced will satisfy ( 0 < ix < 2147483647 ).         *
                                                                           *
 The generator is full cycle, every integer from 1 to (2**31)-2 =          *
      2147483646 is generated exactly once in the cycle.                   *
                                                                           *
 Input: integer "ix", ( 0 < ix < 2147483647 )                              *
 Ouput: real "xrand", ( 0 < xrand < 1 )                                    *
***************************************************************************/
float
randp(int *ix)
{
  int xhi, xalo, leftlo, fhi, k;

  const int A   = 16807;  /* = 7**5          */
  const int P   = 2147483647; /* = Mersenne prime (2**31)-1  */
  const int b15   = 32768;  /* = 2**15                 */
  const int b16 = 65536;  /* = 2**16                 */

  /* get 15 hi order bits of ix */
  xhi = *ix/b16;

  /* get 16 lo bits of ix and form lo product */
  xalo  = (*ix-xhi*b16)*A;

  /* get 15 hi order bits of lo product */
  leftlo  = xalo/b16;

  /* from the 31 highest bits of full product */
  fhi = xhi*A+leftlo;

  /* get overflo past 31st bit of full product */
  k = fhi/b15;

  /* assemble all the parts and presubtract P */
        /* the parentheses are essential            */
  *ix = (((xalo-leftlo*b16)-P)+(fhi-k*b15)*b16)+k;

  /* add P back in if necessary  */
  if (*ix < 0) *ix = *ix + P;

  /* multiply by 1/(2**31-1) */
  return (float)(*ix*4.656612875e-10);
}

/***********************************************************
 To generate a random integer "k" in [i,j].                *
 Input: integers "seed", "i" and "j" in [1,2147483646]     *
 Ouput: integer in [i,j]                                   *
************************************************************/
int
get_rand_ij(int *seed, int i, int j )
{

   randp(seed);

   return ((double)*seed/((double)2147483647/((double)(j-i+1))))+i;

}

/**************************************************************
 To generate a random integer "k" in [1,size].                *
 Input: integers "seed" and "size" in [1,2147483646]          *
 Ouput: integer in [1,size]
                                   *
**************************************************************/
int
get_rand(int *seed, int size )
{
   randp(seed);

   return ((double)*seed/((double)2147483647/((double)size)))+1;

}


/*************************************************
 To check the correctness of the implementation. *
 Output: correct (1) or wrong (0)                 *
**************************************************/
int
trand()
{
  int i, seed;

  seed = 1;

  for (i=0;i<1000;i++)
    randp(&seed);

  if ( seed == 522329230 )
      return 1;
  else
      return 0;

}