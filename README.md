# VLBT: An adaptive encoding for BWTs and compressed suffix arrays 

This repository implements run-length BWTs and BWT-based compressed suffix arrays using *variable-length blocking* (VLB),
a novel technique that exploits the skew distribution of BWT runs to balance space usage and query speed.

VLB constructs an unbalance tree over the run-length BWT (the VLB-tree). Compressible BWT areas (i.e., 
few runs spanning a large segment) are fast to operate, so the tree stores little indexing information for them. The 
tree "reallocates" this spared space to add more indexing information in incompressible BWT areas (many short runs 
covering a shoft segment), where operating is more costly. In the tree, answering rank and successor queries (core 
operations in pattern matching) involves descending the tree until a leaf is reached and then performing a 
cache-friendly scan of a bounded number of BWT runs in the leaf.

Compressible areas are placed near the tree root and are fast to access (close to one cache miss), 
while incompressible areas are placed at deeper levels. Although deeper levels may trigger more cache misses during a 
descent, they contain information that speeds up the access to specific BWT positions in incompressible areas. 
Besides, cache locality improves at deeper tree levels.

The VLB-tree is a shallow tree, and its maximum height is bounded by a logarithmic factor (large base). You can 
think of the VLB-tree as a B-tree-like structure tailored for the skew distribution of BWT runs.

The VLB-tree can also place suffix array samples near to their corresponding BWT runs so you can efficiently get 
the lexicographically smallest occurrence of the queried pattern.  

Additionally, the VLB-tree can be used to encode the function $\phi^{-1}(S\!A[j])=S\!A[j+1]$, necessary to decode 
the rest of occurrences. Both trees (from the BWT and $\phi^{-1}$) conform a fully-functionally compressed suffix 
array. We also provide an implementation the sr-index, the fast variant storing valid index areas to speed up the
query time. 

## TL;DR

## Motivation

A compressed suffix array is a data structure that stores a text in compressed form at the time it allows counting and 
locating occurrences of a given pattern in the text. 

This idea has many forms, but probably the most popular ones are those based on the Burrows-Wheeler Transform (BWT). 
Combining the BWT of the text with some samples of the suffix array allows building the so-called FM index, the 
algorithmic workhorse behind popular bioinformatics tools such as BWA-MEM and Bowtie2. 

In practical implementations, the FM index usually uses space proportional to the plain text, which is not ideal for 
large texts. So, Gagie et al. created the r-index, a compressed version of the FM index whose space cost is 
proportional to the number of runs in the BWT. When the text is highly compressible (say, copies of nearly identical
sequences), the r-index can be much smaller than the FM index, thus motivating its use in pangenomics applications. 

However, the r-index has an important limitation: the number of runs in the BWT is highly sensitive to edits, meaning 
that as soon as we introduce more variation to the text, the size of the r-index inflates quickly to reach the 
size of the FM index. The sr-index addresses this problem by partially removing information from the index that is 
later recomputed on the fly during query time, showing a significant improvement in space under not-so-repetitive 
scenarios. 

Another relevant problem is that performing pattern matching on the r-index (and sr-index) requires the 
interplay of multiple composition data structures that lack spatial locality, making the process relatively slow 
compared to plain alternatives. More recent encodings (the move data structure) partially alleviate the locality problem,
using a much more straightforward layout that sacrifices space efficiency for speed. Overall, state-of-the-art 
encodings for BWT-based compressed suffix arrays either prioritize space efficiency or speed, limiting their 
applicability in terabyte-scale applications. 

