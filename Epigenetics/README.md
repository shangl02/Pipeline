# Filter
This pipeline predicts variants effect on chromatine openness. It uses [scEpiLock](https://github.com/aicb-ZhangLabs/scEpiLock) to train the model using single cell ATACSeq data. One good source of scATACseq can be download from this paper: [A single-cell atlas of chromatin accessibility in the human genome](https://www.sciencedirect.com/science/article/pii/S0092867421012794).

## How to run the pipeline

	python variant_predict/variant_impact.py \
		-m model.pt \
		-r hg38.fasta \
		-s snp_file.txt \
		-o outputfile.txt

* -m: model file produced using pytorch.
* -r: reference fasta file
* -s: snp file, each line has one snp, it needs to have 9 columns. [peak_chrom, peak_start, peak_end, snp_chrom, snp_start, snp_end, snp_id, ref_allele, alt_allele]
* -o: output file that store the results.
