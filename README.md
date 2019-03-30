# Requirements

* bedtools (`sudo apt install bedtools`)
* jellyfish (if genome counts will be calculated) (`sudo apt install jellyfish`)
* `ncbi-blast+` package for `dustmasker` utility (`sudo apt install ncbi-blast+`)
* R packages: `optparse`, `tidyverse`, `lsa`

# Usage

Before using the `Makefile` please update `VCF` and `GENOME` arguments according to your own project.

Then test if coordinates match or not and finally `make`

```bash
make test_coordinate
make analysis
make report
```
