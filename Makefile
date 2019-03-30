# complexity calculation is done by dustmasker, to install: sudo apt install ncbi-blast+
# word counts is done by jellyfish, to install: sudo apt install jellyfish

GENOME = ../GRCh37.p13/Homo_sapiens.GRCh37.dna.primary_assembly.fa
VCF = input/human-common_20180423_test.vcf.gz
WORDVECTOR = input/dna2vec-20161219-0153-k3to8-100d-10c-29320Mbp-sliding-Xat.w2v.gz

# all: snp151_hg37_8mer_overlap_count_complexity_genomcount.csv
# .PHONY: all

snp151_hg37_8mer_overlap_count_complexity_genomcount.csv: snp151_hg37_8mer_overlap_count hg37_08mer_counts
	@echo "adding complexity and genome counts to overlap counts.."
	@echo
	@awk -f scripts/generate_fasta_from_kmer $< | \
	dustmasker -in - -window 8 -level 8 -outfmt acclist | \
	awk -f scripts/add_dust_output - $< | \
	awk -f scripts/add_genome_count hg37_08mer_counts - > $@

snp151_hg37_8mer_overlap_count:
	@echo "generating overlap count.. might take a while.."
	@echo
	@zcat $(VCF) | \
	egrep -v "^#" | \
	awk -f scripts/generate_windows_vcf | \
	bedtools getfasta  -name -tab -fi $(GENOME) -bed - | \
	tr ":-" "\t\t" | \
	awk -f scripts/count_snp_per_coordinate | \
	awk -f scripts/add_neighbor_to_variant | \
	grep -v N | \
	sort > $@

hg37_08mer_counts:
	@echo "genome counts file is missing.. generating counts with jellyfish.."
	@echo
	jellyfish count -m 8 -s 100000000 -t 7 $(GENOME)
	jellyfish dump mer_counts.jf | awk '/>/{count=$$0; getline; gsub(/>/,"",count); printf "%s\t%s\n",$$0,count}' > $@
	rm mer_counts.jf

.PHONY: clean clean_all test_coordinate

test_coordinate:
	@echo "testing if coordinates match between vcf file and genome"
	@zgrep -v "^#" ${VCF} | \
	awk 'length($$4)==1 && length($$5)==1' | \
	shuf | \
	head -1000 | \
	awk '{printf "%s\t%s\t%s\t%s-%s-%s\t.\t.\n",$$1,$$2-1,$$2,$$3,$$4,$$5}' >| _tmp_coordinates

	@bedtools getfasta -fi ${GENOME} -bed _tmp_coordinates -name -tab | tr "-" "\t" | awk '$$2!=$$5' > _tmp_coordinates2

	@echo "VCF file coordinates and genome coordinates"
	@if [ -s _tmp_coordinates2 ]; then echo "DO NOT MATCH"; else echo "MATCH"; fi 
	@rm _tmp_coordinates _tmp_coordinates2

clean:
	rm snp151_hg37_8mer_overlap_count
	rm snp151_hg37_8mer_overlap_count_complexity_genomcount.csv

clean_all: clean
	rm hg37_08mer_counts
