# complexity calculation is done by dustmasker, to install: sudo apt install ncbi-blast+
# word counts is done by jellyfish, to install: sudo apt install jellyfish

GENOME = ../GRCh37.p13/Homo_sapiens.GRCh37.dna.primary_assembly.fa
VCF = ../human-common_20180423.vcf.gz

# all: snp151_hg37_8mer_overlap_count_complexity_genomcount.csv
# .PHONY: all

snp151_hg37_8mer_overlap_count_complexity_genomcount.csv: snp151_hg37_8mer_overlap_count hg37_08mer_counts
	@echo "adding complexity and genome counts to overlap counts.."
	@echo
	@awk -f generate_fasta_from_kmer $< | \
	dustmasker -in - -window 8 -level 8 -outfmt acclist | \
	awk -f add_dust_output - $< | \
	awk -f add_genome_count hg37_08mer_counts - > $@

snp151_hg37_8mer_overlap_count:
	@echo "generating overlap count.. might take a while.."
	@echo
	@zcat $(VCF) | \
	egrep -v "^#" | \
	awk -f generate_windows_vcf | \
	bedtools getfasta  -name -tab -fi $(GENOME) -bed - | \
	tr ":-" "\t\t" | \
	awk -f count_snp_per_coordinate | \
	awk -f add_neighbor_to_variant | \
	grep -v N | \
	sort > $@

hg37_08mer_counts:
	@echo "genome counts file is missing.. generating counts with jellyfish.."
	@echo
	jellyfish count -m 8 -s 100000000 -t 7 $(GENOME)
	jellyfish dump mer_counts.jf | awk '/>/{count=$$0; getline; gsub(/>/,"",count); printf "%s\t%s\n",$$0,count}' > $@
        rm mer_counts.jf

.PHONY: clean clean_all
clean:
	rm snp151_hg37_8mer_overlap_count
	rm snp151_hg37_8mer_overlap_count_complexity_genomcount.csv

clean_all: clean
	rm hg37_08mer_counts
