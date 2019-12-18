# Requirements

* bedtools (`sudo apt install bedtools`)
* kmc (if genome counts will be calculated) (`sudo apt install kmc`)
* `ncbi-blast+` package for `dustmasker` utility (`sudo apt install ncbi-blast+`)
* R packages: `optparse`, `tidyverse`, `lsa`, `DT`, `plotly`
* perl script [SeqComplex](https://github.com/caballero/SeqComplex) (**note**: update Makefile to use `ce` entropy column from file)

# Usage

Before using the `Makefile` please update `VCF` and `GENOME` arguments according to your own project.

Then test if coordinates match or not and finally `make`

```bash
make test_coordinate
make analysis
make report
```

## Docker usage

If running locally the genome files can be located at any folder in your computer. However, if you are using Docker image to run the analysis the genome file (and all other files) should reside under mounted directory.

```
docker run --rm -it --user 1000 -v $(pwd):/app alperyilmaz/snp-dna2vec make analysis
```

> Docker image does not contain SeqComplex script yet
