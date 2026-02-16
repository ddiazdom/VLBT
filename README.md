# VLBT: an adaptive encoding for BWTs and compressed suffix arrays

This repository provides implementations of run-length BWTs and BWT-based compressed suffix arrays 
(CSA) leveraging *variable-length blocking* (VLB), a novel technique that exploits the skew distribution of BWT 
runs to balance space usage and query speed.

## TL;DR

State of the art in practical CSAs:

* The $sr$-index is the most space-efficient BWT-based CSA.
* The move data structure is the most query-efficient BWT-based CSA.

When comparing their time and space tradeoffs, they are at opposite ends of the Pareto frontier.

This repository implements two data structures:

* A run-length BWT with support for $count$ queries.
* The fast variant of the $sr$-index (valid areas). 

The key takeaway is this:

Our VLB framework strikes a balance: its space usage is comparable to the $sr$-index,
but it is substantially faster. On the other hand, it remains slower than move-based structures while consuming 
significantly less space. This tradeoff makes VLBT practical for pangenomics and similar applications, where the input
text is highly repetitive, but sequence variation affects BWT compressibility.

VLBT is a promising alternative, as it can effectively handle variation to produce compact representations—essential for
terabyte-scale inputs—while still supporting fast pattern-matching queries.

## Motivation

A compressed suffix array (CSA) is a data structure that stores a text in compressed form and allows counting and
locating occurrences of a given pattern in the text.

This idea takes many forms, but the most popular are those based on the Burrows-Wheeler Transform (BWT).
Combining the BWT of the text with some samples of the suffix array allows building the so-called FM index, the
algorithmic workhorse behind popular bioinformatics tools such as [BWA-MEM](https://github.com/lh3/bwa) and
[Bowtie2](https://github.com/BenLangmead/bowtie2).

In practical implementations, the FM index usually uses space proportional to the plain text, which is not ideal for
large texts. So, Gagie et al. created the $r$-index, a compressed version of the FM index whose space cost is
proportional to the number of runs in the BWT. When the text is highly compressible (say, copies of nearly identical
sequences), the $r$-index can be much smaller than the FM index, thus motivating its use in pangenomics applications.

However, the $r$-index has an important limitation: the number of runs in the BWT is sensitive to edits, meaning
that as soon as we introduce variation to the text, the size of the $r$-index inflates quickly to reach the
size of the FM index. The $sr$-index addresses this problem by partially removing information that is later recomputed
on the fly during query time, showing a significant improvement in space under not-so-repetitive scenarios.

Another relevant problem is that performing pattern matching on the $r$-index (and $sr$-index) requires the
interplay of multiple internal data structures that lack spatial locality, making the process relatively slow
compared to plain alternatives. More recent encodings (the move data structure) partially alleviate the locality problem,
using a much more straightforward layout that sacrifices space efficiency for speed. Overall, state-of-the-art
encodings for BWT-based CSAs either prioritize space efficiency or speed, limiting their applicability in terabyte-scale
applications.

## Design principle

A way to deal with the space issue is to group BWT runs into blocks, storing global indexing information about the 
blocks, and recomputing the missing information on the fly during query time. This idea, in principle, should keep 
variation-related space overhead controlled. The challenge is to find a suitable way to distribute the indexing 
information across the runs such that we still achieve good query performance.

VLB constructs an *unbalanced shallow* tree over the run-length BWT (the VLB-tree), where the leaves encode
variable-length BWT blocks and the internal nodes store indexing information that speeds up access to those blocks.
Answering rank and successor queries (core operations in pattern matching) involves descending the tree until a leaf is
reached and then performing a cache-friendly scan of a bounded number of runs in the leaf.

The key feature of our design is that compressible areas are placed near the tree root and are fast to access 
(close to one cache miss), while incompressible areas are placed at deeper levels. Incompressible areas contain many 
shorts that are more costly to access with a linear scan. However, the indexing information in the 
internal nodes of the path allows skipping many runs, improving the query performance. Deeper nodes trigger more cache 
misses, but they are still faster than a linear scan. You can think of the VLB-tree as a data structure that 
relocates space from compressible BWT areas to incompressible ones. 

The VLB-tree can also place suffix array samples near their corresponding BWT runs, so you can efficiently get
the lexicographically smallest occurrence of the queried pattern.

Additionally, the VLB-tree can be used to encode the function $\phi^{-1}(SA[j])=SA[j+1]$, necessary to decode
the rest of the occurrences. Both trees (from the BWT and $\phi^{-1}$) form a fully functional CSA. We also provide an
implementation of the $sr$-index, with the fast variant that speeds up $locate$ queries.

## Dependencies

* C++17 compiler
* CMake

So far, we have tested VLBT on Linux and macOS, using GCC x and Clang x. We do not guarantee that VLBT will work on 
other platforms, yet.

## External repositories

* [CLI](https://github.com/CLIUtils/CLI11) (already included in the repository)

## How to build

Enter the repository directory and run:

```
$ mkdir build
$ cd build
$ cmake ..
$ make
```

This process will generate a binary called `vlbt-cli` (among other things) in the `build` directory, which allows
building and querying indexes. This binary is **not** a full-fledged sequence aligner. It is meant to be used for
indexing and testing the performance of $count$ and $locate$ queries. 

## Block size

We use a block size $\ell$ that guides the shape of the VLB-tree. In general, small values should 
increase the space but improve query speed, while large values should have the opposite effect. 
However, this behavior is not strict, as $\ell$ is a reference value that the construction algorithm changes with the
local run structure in the BWT. The performance should not vary substantially as we change $\ell$, assuming it
is large enough. 

We limited $\ell$ in the implementation to powers of $4$ in $4^{5}–4^{10}$. This range is fairly wide to cover 
repetitive and non-repetitive texts, even at a large scale.

Here is a general rule of thumb to decide its value:

* Text has near-identical sequences: X
* Text highly repetitive but with more variation (e.g., metagenomes): X
* Text is highly repetitive and has a large alphabet: X

These values are rough and based on our experiments. You can explore others. We would like to devise a 
mechanism to recommend a suitable $\ell$ based on the distribution of BWT runs, but that is future work.

## Building VLBT data structures

For the moment, we do not provide a mechanism to compute the BWT and/or the $2r$ suffix array samples. These components 
are necessary but have to be computed externally. In the meantime, you can use [BigBWT](https://gitlab.com/manzai/Big-BWT)
to produce these files. Our cli expects input files in their format. 

### Run-length BWT

To create the run-length BWT, you have to run

```
./vlbt-cli build mytext.txt.bwt -b 4096 -d 0 
```

Where `mytext.txt.bwt` is the BWT of `mytext.txt` in one-byte-per-symbol encoding (i.e., plain). The `-b` option 
specifies the block size, while the `-d` option specifies the structure we are building (0 means run-lenth BWT).
The option `--help` gives more details.

### BWT-based CSA 

In this case, you also have to have `mytext.txt.bwt` beforehand, but also the files 

* `mytext.txt.ssa`
* `mytext.txt.esa`

The first (`ssa` extension) stores the suffix array samples corresponding to BWT run heads, and the second (`esa` 
extension) stores the suffix array samples corresponding to BWT run tails. Both files must store the samples using
five bytes per symbol and in suffix array order. Notice that if a BWT run has length $1$, it is simultaneously a head 
and a tail. In this case, the corresponding suffix array sample has to be in both files. 

The command to build the CSA is 

```
./vlbt-cli build mytext.txt -b 4096 -d 2 -s 5 
```

Where `-d 2` indicates that we are building the CSA and `s` is the subsampling parameter of the $sr$-index. The CLI 
will look for files `mytext.txt.bwt`, `mytext.txt.ssa`, and `mytext.txt.esa` in the same directory as `mytext.txt`.
Our VLB-based CSA for the moment uses the same block size $\ell$ for both the BWT and $\phi^{-1}$. This may change in
the future.  

## Querying an index:

To $count$ the occurrences of a pattern in an indexed text, use the command 

```
./vlbt-cli count index.vlbt pat_file
```

where `index.vlbt` is the VLBT index (run-length BWT or CSA) and `pat_file` is the pattern file in
[Pizza&Chilli](https://pizzachili.dcc.uchile.cl/utils/genpatterns.c) 
format.

To $locate$ the occurrences of a pattern, use 

```
./vlbt-cli locate index.vlbt pat_file
```
The input index in this case must be a VLBT CSA.

**Note:** this interface performs the queries, but it only reports statistics (speed, number of occurrences, 
etc.). See below how to actually get the $locate$ results. 

## Including VLBT in your project: 

It is also possible to include VLBT as a library in your own project. Copy the `include` folder and add the following 
lines to your source files:

Run-length BWT:
```C++
#include "include/vlbt_build_bwt.h"
#include "include/vlbt_bwt.h"

int main() {
    
    //build the index
    vlbt_rlbwt<4096> bwt;//block size as a template parameter
    build_bwt(bwt, input_bwt_file, PLAIN, "/tmp/folder");//PLAIN means BWT format
    
    //store the index to disk
    store_to_file("/path/to/bwt_index", bwt);
    
    //count occurrences
    std::string pattern = "atggagag";
    size_t count = bwt.count(pattern); 
   
    //load from disk
    vlbt_rlbwt<4096> bwt2;
    load_from_file("/path/to/bwt_index", bwt2);//make sure template parameters match 
}
```

The construction algorithm will place temporary files in `/tmp/folder`. It then will delete them. 

CSA:
```C++
#include "include/vlbt_sr_index.h"
#include "include/vlbt_build_sr_index.h"

int main() {
    
    size_t s = 5;//subsampling parameter
    vlbt_sri_va<4096, 4096> csa;//left is block size for the BWT and right for $\phi^{-1}$
    build_sr_index<uint64_t>(csa, bwt_file, PLAIN, s, "/tmp/folder");//SA samples are stored in uint64_t cells
    
    //store to disk
    store_to_file("/output/csa/index", csa);

    std::string pattern = "atggagag";
    
    //count
    size_t count = csa.count(pattern); 
    //locate
    std::vector<size_t> positions = csa.locate(pattern);
    
    //load from disk
    vlbt_sri_va<4096, 4096> csa2;
    load_from_file("/output/csa/index", csa2);//make sure template parameters match 
}
```

We have not tested using VLBT as a library yet, but it should work. If not, please open an issue. The 
fix should be straightforward.

## Experimental results

### Datasets:

* BAC: genome assemblies of 30 bacterial species from the [AllTheBacteria]() collection. Strings 
belonging to the same species are highly repetitive, while strings from different species are dissimilar.
 
* COVID: $4{,}494{,}508$ SARS-CoV-2 genomes downloaded from the [NCBI](https://uud.ncbi.nlm.nih.gov/home/genomes) genome
  portal. These sequences are short and near-identical, with an average length of $29{,}748$.

* HUM: genome assemblies of 40 individuals from the Human Pangenome Reference Consortium (HPRC).
  Individual genomes are near-identical, but assembly differences introduced variability.

* KERNEL: $2{,}609{,}417$ versions of the [Linux kernel](https://github.com/torvalds/linux) repository. This dataset is
  highly repetitive and has a large alphabet.

In DNA collections (BAC, COVID, and HUM), we also considered the DNA reverse complement of each string (as is standard in
bioinformatics). The numbers presented in the table below already consider these extra sequences.

| Dataset | Size (GB) | Alphabet | $n/r$  | Longest BWT run (MB) |
|---------|-----------|----------|--------|----------------------|
| BAC     | 133.12    | 7        | 116.77 | 0.23                 |
| COVID   | 267.41    | 17       | 940.49 | 4.11                 |
| HUM     | 241.24    | 7        | 61.82  | 6.66                 |
| KERNEL  | 54.45     | 190      | 263.11 | 70.6                 |

### Competitor tools

 * [mn](https://github.com/simongog/sdsl-lite/blob/master/include/sdsl/wt_rlmn.hpp) (release 2.1.1): 
 the run-length BWT of Mäkinen and Navarro as implemented in the [SDSL library](https://github.com/simongog/sdsl-lite).
 * [fbb](https://github.com/dominikkempa/faster-minuter) (commit 9238178): BWT encoding using 
   fixed block boosting.
 * [movc]() and [movl](https://github.com/LukasNalbach/Move-r) (commit bed2fe9): optimized 
   implementations of the move data structure. The variant [movc]() encodes only the BWT, while [movl]() also 
   includes $r$-suffix array samples.
 * [ri](https://github.com/nicolaprezza/r-index) (commit 7009b53): the original $r$-index.
 * [sri-va](https://github.com/duscob/sr-index) (commit f99b54a): the original $sr$-index with 
   valid $\phi^{-1}$ areas. We varied the sampling $s$ across values $8,12,16,20$.

### Count queries in run-length BWTs:

Random patterns of length 105 were generated using [Pizza&Chilli](https://pizzachili.dcc.uchile.cl/utils/genpatterns.c).
The table shows query speed in microseconds per pattern (μs/pat) and index space usage in bits per symbol (bps). 

| run-length BWT    | 30bac |   30bac | 40hum |  40hum | covid |  covid | kernel | kernel |
|:------------------|------:|--------:|------:|-------:|------:|-------:|-------:|-------:|
|                   |   bps |  μs/pat |   bps | μs/pat |   bps | μs/pat |    bps | μs/pat |
| vlbt-bwt_b_4096   | 0.142 |   45.83 | 0.332 |  52.99 | 0.024 |  33.27 |  0.172 |  44.85 |
| vlbt-bwt_b_16384  | 0.133 |   53.72 | 0.317 |  75.52 | 0.017 |  34.14 |  0.127 |  39.62 |
| vlbt-bwt_b_65536  | 0.129 |   72.25 | 0.313 |  91.02 | 0.016 |  41.08 |  0.109 |  49.93 |
| vlbt-bwt_b_262144 | 0.128 |   83.05 | 0.313 | 101.44 | 0.015 |  49.85 |  0.104 |  60.65 |
| fbb               | 0.232 |  102.33 | 0.273 | 121.73 | 0.069 |  71.73 |  0.179 |  69.77 |
| mn                | 0.192 |  190.25 | 0.331 | 204.25 | 0.031 | 184.80 |  0.115 | 204.06 |
| movc              | 0.962 |   17.99 |    NA |     NA | 0.125 |  11.96 |  0.438 |  15.53 |

### Locate queries:

The following table shows the average time (in seconds) to count and locate $10^{10}$ occurrences of a pattern in the 
datasets. 

| CSA                | 30bac |  30bac |   hum |    hum | covid |  covid | kernel |  kernel |
|:-------------------|------:|-------:|------:|-------:|------:|-------:|-------:|--------:|
|                    |   bps | μs/occ |   bps | μs/occ |   bps | μs/occ |    bps |  μs/occ |
| movloc             | 2.038 |  0.092 |    NA |     NA | 0.255 |  0.148 |  0.946 |   0.105 |
| ri                 | 0.824 |  1.406 | 1.552 |  1.545 | 0.109 |  0.741 |  0.381 |   2.544 |
| sri-s8             | 0.363 |  0.765 | 0.627 |  1.086 | 0.096 |  0.600 |  0.232 |   1.543 |
| sri-vm-s8          | 0.365 |  0.771 | 0.632 |  1.089 | 0.097 |  0.612 |  0.233 |   1.543 |
| sri-va-s8          | 0.397 |  0.754 | 0.685 |  1.028 | 0.098 |  0.609 |  0.258 |   1.490 |
| vlbt-sri-va-b6-s8  | 0.392 |  0.337 | 0.762 |  0.344 | 0.117 |  0.415 |  0.341 |   0.342 |
| vlbt-sri-va-b7-s8  | 0.378 |  0.250 | 0.742 |  0.387 | 0.104 |  0.362 |  0.291 |   0.273 |
| vlbt-sri-va-b8-s8  | 0.373 |  0.262 | 0.737 |  0.454 | 0.101 |  0.357 |  0.271 |   0.283 |
| vlbt-sri-va-b9-s8  | 0.372 |  0.292 | 0.736 |  0.510 | 0.100 |  0.473 |  0.266 |   0.335 |
| sri-s16            | 0.300 |  0.691 | 0.492 |  1.307 | 0.091 |  0.647 |  0.182 |   1.675 |
| sri-vm-s16         | 0.301 |  0.703 | 0.495 |  1.093 | 0.092 |  0.627 |  0.183 |   1.674 |
| sri-va-s16         | 0.322 |  0.648 | 0.528 |  1.044 | 0.094 |  0.616 |  0.203 |   1.461 |
| vlbt-sri-va-b6-s16 | 0.325 |  0.317 | 0.607 |  0.345 | 0.112 |  0.411 |  0.285 |   0.340 |
| vlbt-sri-va-b7-s16 | 0.311 |  0.237 | 0.587 |  0.387 | 0.099 |  0.365 |  0.234 |   0.276 |
| vlbt-sri-va-b8-s16 | 0.306 |  0.248 | 0.583 |  0.457 | 0.096 |  0.354 |  0.215 |   0.282 |
| vlbt-sri-va-b9-s16 | 0.305 |  0.279 | 0.582 |  0.513 | 0.095 |  0.467 |  0.210 |   0.342 |

## How to cite