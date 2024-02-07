# fbrl-bwt: an efficient data structure to encode BWTs with small alphabets

This data structure is designed to encode massive and repetitive text collections with tiny alphabets (16 symbols or
less),
such as DNA data. The main advantages of fbrl-bwt are:

* Efficient construction (fast and low memory footprint)
* High compression ratios
* Fast query performance

## Design

Our data structure follows the scheme of (x et al., ESA 2023): we keep the input BWT[1..n] in run-length format.
Additionally, we choose a parameter b to divide the BWT into ceil(n/b) blocks of length b (those BWT runs crossing two
or more blocks are split). For each block, we store precomputed rank answers and store them as a header for the block.
If the block contains more than s runs (s being a parameter), we divide the block into ceil(b/m) mini blocks of length
m (also a parameter), also storing precomputed rank answers for the mini block.

The BWT runs of each block (or mini block) use one or two bytes, depending on the block content. The maximum length a
byte can encode is defined by the BWT alphabet. For instance, an alphabet=16 requires 4 bits, so one byte can encode
runs of length up to 16 symbols (because we have 8-4 bits (2^4=16) available in the byte). For two-byte cells,
we can encode run lengths up to 4096 (log2(4096)=16-4=12). Notices that the alphabet also sets a cap on the maximum
block size b we can choose. The smaller the alphabet, the longer are the runs we can fit in one byte and the longer are
the blocks we can set. On the other hand, by choosing a low value for s, say 128, the blocks tend to have fewer but
longer runs, so we always use two bytes for the blocks. Mini blocks are built over fragmented blocks (i.e., they contain
several runs). In this case, we check if its better two use two bytes or one byte and break the longer runs to fit one
byte. We store in the header of the mini block which encoding we are using.

## Usage

## Benchmarks

We compared our data structure against different data structures aimed to encode the BWT with text indexing queries:

* **wt_huff_bv** : the huffman-shaped wavelet tree whose underlying scheme uses plain bit vectors
* **wt_huff_bv** : the run-length-encoded BWT of Mäkinen and Navarro (2005?)
* **wt_fbb_hyb** : the wavelet tree that uses the fix-block technique to boost compression (Gog et al., DCC 2016), with
  the hybrid vector as the underlying bit vector representation.

We used the implementations of those data structures available in the SDSL-library. We used a machine x and x. All the
data structures used the compiler flags ''.

## Construction and space usage

### **Dataset**: sampled PacBio reads from the E.coli genome

The size of this collection is 342.02 MB and has r=9547033 equal-symbols runs in its BCR BWT (n/r=35.8249). The alphabet
is 9 as the equal-symbol runs of length >1 in the text (not in the BWT) were transformed into metasymbols. The '*'
character denotes the best performance in the tables above.

|                               | fb_rl   | wt_huff_bv | wt_rlmn | wt_fbb_hyb | 
|-------------------------------|---------|------------|---------|------------|
| build time (ss.ms)            | *00.491 | 20.655     | 32.364  | 32.88      |
| space usage (bits per symbol) | 0.560   | 4.355      | 0.596   | *0.461     |

## Query performance

The average time in nanoseconds of standard BWT queries. The inputs for those queries
were chosen at random. Data structures with a "-" symbol do not implement that query

|                  | fb_rl    | wt_huff_bv | wt_rlmn | wt_fbb_hyb | 
|------------------|----------|------------|---------|------------|
| interval symbols | *732.601 | 1643.26    | -       | -          |
| select           | *1158.08 | 1869.65    | 2097.77 | -          |
| inverse select   | *361.209 | 602.398    | 1135.21 | 817.013    |
| access           | *275.389 | 598.574    | 549.113 | 738.682    |
| rank             | *358.273 | 695.783    | 954.587 | 625.148    |

## Citation