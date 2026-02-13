# VLBT: an adaptive encoding for BWTs and compressed suffix arrays

This repository provides implementations of encodings for the run-length BWT and BWT-based compressed suffix array, leveraging
variable-length blocking (VLB), a novel technique that exploits the skew distribution of BWT runs to balance
space usage and query speed.

VLB constructs an unbalanced tree over the run-length BWT (the VLB-tree). Compressible BWT areas (i.e.,
few runs spanning a large segment) are fast to operate, so the tree stores little indexing information for them. The
tree “reallocates” this spared space to add more indexing information in incompressible BWT areas (many short runs
covering a short segment), where operating is more costly. In the tree, answering rank and successor queries (core
operations in pattern matching) involves descending the tree until a leaf is reached and then performing a
cache-friendly scan of a bounded number of BWT runs in the leaf.

Compressible areas are placed near the tree root and are fast to access (close to one cache miss),
while incompressible areas are placed at deeper levels. Although deeper levels may trigger more cache misses during a
descent, they contain information that speeds up the access to specific BWT positions in incompressible areas.
Additionally, cache locality improves at deeper levels of the tree.

The VLB-tree is a shallow tree, and its maximum height is bounded by a logarithmic factor (large base). You can
think of the VLB-tree as a B-tree-like structure tailored for the skew distribution of BWT runs.

The VLB-tree can also place suffix array samples near their corresponding BWT runs, so you can efficiently get
the lexicographically smallest occurrence of the queried pattern.

Additionally, the VLB-tree can be used to encode the function $\phi^{-1}(S!A[j])=S!A[j+1]$, necessary to decode
the rest of the occurrences. Both trees (from the BWT and $\phi^{-1}$) form a fully-functional compressed suffix
array. We also provide an implementation of the sr-index, the fast variant storing valid index areas to speed up the
query time.

## TL;DR

We implement two data structures: the run-length BWT and the fast variant of the $sr$-index.

State of the art in practical CSAs:

* The $sr$-index is the most space-efficient BWT-based CSA variant in the literature.
* The move data structure is the most query-efficient BWT-based CSA variant.

When comparing their tradeoffs, they are at opposite ends of the Pareto frontier.

The key takeaway is this:

Our VLBT-based CSA implementation strikes a balance between these methods: its space usage is comparable to the $sr$-index, but it is substantially faster. While the move data structure remains faster, it consumes significantly more space. This tradeoff makes VLBT practical for pangenomics and similar applications. In such scenarios, BWT-based CSAs remain the most efficient option for pattern matching in lossless compressed space.
However, current data structures are still too large because pangenomes and metagenomes contain substantial variation, such as misassemblies and genetic diversity.

VLBT is a promising alternative, as it can effectively handle variation to produce compact representations—essential for terabyte-scale inputs—while still supporting fast pattern-matching queries.

## Motivation

A compressed suffix array is a data structure that stores a text in compressed form, which allows counting and
locating occurrences of a given pattern in the text.

This idea takes many forms, but the most popular are those based on the Burrows-Wheeler Transform (BWT).
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