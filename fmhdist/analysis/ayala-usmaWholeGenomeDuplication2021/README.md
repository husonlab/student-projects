# Summary

The content of this directory is not directly included in the thesis. I did the
same analysis with dataset **C** (`../ref-db`). However, I initially planned to
do the following:

Here, I want to analyze the genomes of _P. infestans_ T30-4, _P. infestans_
RC10-1, _P. ramorum_ and _P. sojae_ in terms of any correlation between hash
density in a given window and effector gene density in that window.

The genomes used are

- GCA_020800215.1 (_P. ramorum_)
- GCA_000149755.2 (_P. sojae_)
- GCA_000142945.1 (_P. infestans_ T30-4)
- GCA_011316315.1 (_P. infestans_ RC1-10)

The genomes are sketched using the MurMur hash function and different hashing
seeds, each sketch comes with the corresponding hash coordinates.

Also, for each genome, the sequence complexity was calculated using
[macle](https://github.com/evolbioinf/macle) by splitting each genome in the
individual sequences stored in the fasta file.

The last ingredient for each genome is the annotation file in GTF format. I use
this to obtain the locations of the effector genes.

For each genome, the density of the hashes is analysed using the
`../../misc/coordinates_complexity_correlation.py` script with non-overlapping
windows of size $w=10000$ and the computed sequence complexities.

It turns out that there are regions with a lower density as one would expect
(expected density is calculated using the scaling parameter $s$ of FracMinHash,
the window size and a user defined range, i.e. 0.5 for $s=2000$ and a window
size of $10000$ implies a expected density range of $5 \pm (0.5 \times 5)$) and
that those regions are linked to low-complexity regions (where low complexity is
defined as everything $\leq 0.5$).

There is, however, no link between the unusual density windows and any coding
gene. Thus, I didn't dive into the effector genes at all. For this, I was
planning to use the list compiled by [Ayala-Usma _et
al._](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-021-08079-y).
