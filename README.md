# HLAoptimize
A module to test ranking and optimizing HLA coverage of vaccine elements. 

## Introduction
This repository contains an example project to find optimal vaccine elements according to specific criteria involving HLA coverage. In particular, following analyses are considered for this problem:
- **Minimal HLA coverage**: Here, vaccine elements are ranked based on their HLA coverage, which is defined by a immunogenicity cut-off.
- **Minimal sets of vaccine elements**: Here, an optimization with genetic algorithms is performed that is finding a minimal set of vaccine elements, that optimize HLA coverage, immunigenicity, binding, error rate and coverage across viruses.

## Installation
- Clone this repository to get all the necessary code
- Navigate to the directory:
```
cd HLAoptimize
```
- To install the CLI and associated functions run:
```
python setup.py develop
```

## Usage

### Minimal HLA coverage
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
HLAoptiGA_cli 'exercise_data.csv' 0.7 'test_output'
```