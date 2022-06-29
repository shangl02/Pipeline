class Config:
    def __init__(self,
                 #peak_bed,
                 #snp_bed,
                 peak_snp_bed,
                 ref_fasta,
                 model_name,
                 model_weight_path,
                 n_class,
                 result_path):
        #self.peak_bed = peak_bed
        #self.snp_bed = snp_bed
        self.peak_snp_bed = peak_snp_bed
        self.ref_fasta = ref_fasta
        self.model_name = model_name
        self.model_weight_path = model_weight_path
        self.n_class = n_class
        self.result_path = result_path

