# Fast and small BWT data structure for small alphabets 


|                  | fb_rl    | wt_huff_bv | wt_rlmn | wt_fbb_hyb | 
|------------------|----------|------------|---------|------------|
| bps              | *0.56    | 4.35       | 0.59    | 2.76       |
| interval symbols | *727.394 | 1648.79    | -       | -          |
| select           | *1140.51 | 1741.59    | 2071.4  | -          |
| inverse select   | *341.331 | 593.935    | 1121.01 | 597.573    |
| access           | *264.998 | 588.346    | 533.264 | 556.793    |
| rank             | *339.992 | 687.775    | 951.653 | 533.144    |