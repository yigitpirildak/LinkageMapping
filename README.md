# Problem Description
DNA sequences that are close to one another on a chromosome are more likely to be inherited together compared to genes that are far apart. This is called genetic linkage. This property is an exception to Mendel's Law of Independent assortment which states that the possibility of each gene being inherited is completely independent of each other.

By crossing a trihybrid (heterozygous at all three loci) individual to a purebreed homozygous individual, one can observe the genotype frequencies of each progeny to determine gene order and calculate recombination frequencies. These recombination frequencies can then be used to estimate distance between each gene.

![](https://bio.libretexts.org/@api/deki/files/5339/Fig7.12.png?revision=1&size=bestfit&width=459&height=380)

# How to use
The script takes a file path as an argument and attempts to classify genotypes, gene order and distances between them, coefficient of coincidence and interference. File format is expected to be:

**Gamete that came from the heterozygous *parent* / **freqeuncy**
ABC 416
abc 426
AbC 72
aBc 78
Abc 3
aBC 4
abC 0
ABc 1

For reference, there are some test samples that can be found under **test** folder.
