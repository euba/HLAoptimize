# HLAoptimize
A Python module to test ranking and optimizing HLA coverage of vaccine elements. 

## Introduction
This repository contains an example project to find optimal vaccine elements according to specific criteria involving HLA coverage. In particular, following analyses are considered for this problem:
- **Maximized HLA coverage**: Here, vaccine elements are ranked based on their HLA coverage, which is defined by a immunogenicity cut-off.
- **Minimal sets of vaccine elements**: Here, an optimization with genetic algorithms is performed that is finding a minimal set of vaccine elements, that optimize HLA coverage, immunigenicity, binding, error rate and coverage across viruses.

## Installation
- Clone this repository to get all the necessary code
- Navigate to the directory with:
```
cd HLAoptimize
```
- To install the CLI and associated functions run:
```
python setup.py develop
```

## Usage

### Maximized HLA coverage
To perform a ranking of vaccine elements that maximize the HLA coverage run:
```
HLAopti_cli 'exercise_data.csv' 'Zika Virus' 0.7 'test_output'
```
For help on the input run 
```
HLAopti_cli -h
```
which will give you the following information:

```
HLAoptimize euba$ HLAopti_cli -h
usage: HLAopti_cli [-h] input virus immunocut outdir

positional arguments:
  input       Specify input csv file with HLA data
  virus       Specify the virus you are interested in
  immunocut   Specify the cut-off value of Immunogenicity
  outdir      Specify an output directory for results

optional arguments:
  -h, --help  show this help message and exit
```

### Minimal sets of vaccine elements
To get an optimized minimum set of vaccine elements using a genetic algorithms with fixed weights run:
```
HLAoptiGA_cli 'exercise_data.csv' 0.7 'test_output'
```
For help on the input run 
```
HLAoptiGA_cli -h
```
which will give you the following information:
```
 HLAoptiGA_cli [-h] input immunocut outdir

positional arguments:
  input       Specify input csv file with HLA data
  immunocut   Specify the cut-off value of Immunogenicity
  outdir      Specify an output directory for results

optional arguments:
  -h, --help  show this help message and exit
```
In case you want to use a genetic algorithm with randomized weights you can use a similar command: 
```
HLAoptiGArand_cli 'exercise_data.csv' 0.7 'test_output'
```

## Details

### Maximized HLA coverage
To find suitable vaccine elements, the HLA coverage is calculated by selected HLAs with an user defined immunogenicity cut-off for a specific user defined virus. The resulting table is then ranked by HLA coverage (number of covered HLAs per vaccine element) and median immunigenicity score of the vaccine element.

The resulting table is then exported as .tsv file to an user defined output directory along with additional plots. The plot includes a bar chart with the HLA coverage per vaccine element of the seleted virus, in which top vaccine elements are highlighted in red. In addition, a scatter plot shows the median immunogenicity over the HLA coverage, in which top vaccine elements are again highlighted in red and points are scaled by the error rate:
<p align="center">
  <img width="460" src="score2d.png">
</p>
It is recommended to investigate to output plots and the ranked .tsv table to make a decision on which vaccine element might be best suited for further investigation.

### Minimal sets of vaccine elements

To find optimal sets of vaccine elements that satisfy different objective functions, a genetic algorithm is used for optimization. The following objective functions need to be satisfied:

- $f_1(x)$: Minimize the set of selected vaccine elements as a maximization of the inverse size
- $f_2(x)$: Maximize the sum of immogenicity score over all HLAs in a vaccine element over all viruses
- $f_3(x)$: Maximize the sum of binding over all HLAs in a vaccine element over all viruses
- $f_4(x)$: Minimize the sum of overall error rate by maximizing the inverse
- $f_5(x)$: Maximize the HLA coverage over all vaccine elements over all viruses

In principle, multi-objective optimization can be performed, but for simplicity all objectives are scalarized and combined into one weighted objective fitness function:

```math
f(x) = \frac{w_1 * f_1(x) + w_2 * f_2(x) + w_3 * f_3(x) + w_4 * f_4(x) + w_5 * f_5(x)}{\sum_{k=1}^{n=5} w_k}
```
this objective function can then be maximized in each iteration of the genetic algorithm.

#### Data transformation

To make the values of the different objective functions comparable to each other, the data is scaled and transformed to lay within a range of floats between 0 and 1 for each objective. This effectively means, that for each individual objective function (column of the transformed matrix) there is a vaccine element (rows of the transformed matrix) that has a value of 0 and a vaccine element that has a value of 1. This makes the vaccine elements and also the objectives numerically comparable to each other. In addition, this has the benefit, that the objective (fitness) function is as well scaled from 0 (theoretically worst solution) to 1 (theoretically best solution), which aids in the interpretation of the results.

The scaled as well as the raw transformed data can be found as .tsv tables in the user specified ouput folder after each run.

#### Genetic algorithm

The genetic algorithm is implemented as a binary genotype with the length of the total number of vaccine elements, where a 1 indicates the presence of the vaccine element and 0 denotes the absence of a vaccine element. The phenotype (search) space of this genotype denotation would then be all possible sets and combinations of vaccine elements. Based on the phenotype, the overall objective function $f(x)$ can be calculated for each individual (solution) in the population.

The exact steps of the genetic algorithm can be summarized as following:
1. Initialize a random set of individual solutions, e.g., sets of varying vaccine elements
2. Calculate $f(x)$ for each individual $x$ in the population
3. Replace the worst (lowest $f(x)$ values) one third of the population with the best solution (highest $f(x)$ value)
4. Mutate each individual based on a mutation rate (number of vaccine elements taken in or out of the solution)
5. Cross-over the individual with the highest fitness with a random individual of the population
6. Keep the genotype of the individual with the highest fitness in the population
7. Repeat from step 2 until the number of simulated generations is reached
8. Find the individual with the highest fitness and report it as the solution

#### Benchmark

The mutation rate has one of the biggest influence on solution finding in genetic algorithms. To optimize the mutation rate, the genetic algorithm has therefore been simulated with different mutation rates. In addition, to benchmark the objective function, the above described weights $w$ were randomized in each generation.

<p align="center">
  <img height="300" src="GA_fixed.png"><img height="300" src="GA_random.png">
</p>

Mutation rates below $m=2$ seem to yield no improvement in the value of the maximum fitness and mutation rates above $m=4$ seem to yield a sharp decline in the value of the maximum fitness. Based on the results, a mutation rate of $m=2$ has been chosen in the implementation, since this yielded the highest fitness values.

In addition, the objective function utilizing randomized weights yields a more stochastic increase in maximum fitness values in comparison to the steady increase of maximum fitness with static weights. However, both approaches can be valueable in finding different possible solutions, as the randomized weights approach is less sensitive towards biases in individual objective functions.

#### Output

In addition to the above described tables with raw and scaled data, the functions for the genetic algorithm return a plot with the maximum fitness increase and a .tsv table with the combined solutions of independent runs. The vaccine elements in the solution are scored based on the number of times they occur in the independent runs.

## Discussion

- As an alternative to using a scalarization approach for the individual objective functions, a multi-objective optimization could be performed with different solution populations that are assigned differing objective functions. Cross-over events across populations could ensure diversity of solutions.
- Test scripts based on the example data could be included in the code to automatically test different scenarios
- Cross-over events in the current code are not strongly implemented and might not have an direct effect on the solution finding.
- Additional benchmarks including the population size, immunogenicity cut-off, or publicly available tools for optimization could be included.
