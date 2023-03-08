# KmerGO ver 2.0

**KmerGO is a user-friendly tool to identify the group-specific sequences on two groups of high throughput sequencing datasets. A sequence that is present, or rich, in one group, but absent, or scarce, in another group is considered “group-specific” here. Furthermore, KmerGO can also be applied to capture trait-associated sequences for continuous-trait dataset.**

KmerGO2 is almost 2x faster than KmerGO.

Github for KmerGO(v1.x.x) is [here](https://github.com/ChnMasterOG/KmerGO).

Please cite: Wang Y, Chen Q, Deng C, Zheng Y and Sun F (2020) KmerGO: A Tool to Identify Group-Specific Sequences With k-mers. Front. Microbiol. 11:2067. doi: 10.3389/fmicb.2020.02067

### Running on Linux (CommandLine)

Download the Linux command version of KmerGO2(x86 64bit): [download link](https://github.com/ChnMasterOG/KmerGO2/releases/download/v2.0.1/KmerGO2_for_linux_x64.zip).

Decompress the file and enter the software path.

Type these commands if you first use KmerGO2:

> sudo chmod +x KmerGO2

> sudo chmod +x ./bin/*

You can type *`./KmerGO2 -h`* for help message.

### Test data

We prepared some FASTA format files which are stored in "test_data/samples".

And trait files can be also found in folder "test_data".

Because of the small size of testing data, **"k=25 MinValue=1" and other parameters as defaults are recommended**.

## Usage

CommandLine version of KmerGO2 supports END-TO-END mode(*`tool=all`*) and FOUR-TOOLS mode(*`tool=kmc3/union/filtering/cap3`*). The main running command is *`./KmerGO2 <tool>`* with following options:

option | required | description 
:----: | :------: | :---------:  
-i     | Yes      | input files path 
-o     | Yes      | output files path 
-t     | No       | a csv file path of trait information (required when tool is filtering or all) 
-e     | No       | temp files path (required when tool is kmc3 or all, default: ./) 
-m     | No       | mode, 0-categorical, 1-continuous (default: 0) 
-k     | No       | k-mer length (k from 14 to 256, default: 40) 
-ci    | No       | minimal K-mer occurring times (default: 2) 
-cs    | No       | maximal K-mer occurring times (default: 65535) 
-n     | No       | number of processes (default: 24) 
-assl  | No       | when mode = 0, logical features ASS value (default: 0.8) 
-p     | No       | numeric(mode=0) or logical(mode=1) features rank sum test p threshold value (default: 0.01) 
-assn  | No       | when mode = 0, numeric features logistic regression ASS value (default: 0.8) 
-corr  | No       | when mode = 1, numeric features coefficient of association ρ threshold value (default: 0.8) 

Example: *`./KmerGO2 all -k 25 -ci 1 -n 2 -i ./test_data/samples/ -t ./test_data/categorical\ trait\ information.csv`*

## Trait information file format

**Categorical Trait**

Example1:

id,trait

1.fasta,A

2.fasta,A

3.fasta,B

...

Example2:

id,trait

SRR1.fastq,Case

SRR2.fastq,Control

SRR3.fastq,Case

...

**Continuous Trait**

Example:

id,trait

SRR1.fasta,1.02

SRR2.fasta,2.35

SRR3.fasta,5.22

...

**Templates can be also found in the folder "test_data".**

## Tips

1) KmerGO2 must be called from its own directory now because otherwise it cannot find the bin/ directory. (refer to [issue6](https://github.com/ChnMasterOG/KmerGO2/issues/6))

2) The directory path reference by -i must have ONLY the relevant .fasta (or .fastq) files in that directory and nothing else.

## Contributions

Thank meshpv for suggestions.

## Contacts and bug reports

Please send bug reports, comments, or questions to [here](https://github.com/ChnMasterOG/KmerGO2/issues).

----------

Last update: 2023-03-08
