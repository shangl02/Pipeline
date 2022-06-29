from pybedtools import BedTool
import pandas as pd
import numpy as np
import os


class FaGenerator():
    def __init__(self, peak_snp_bed, ref_fasta, result_path):
        self.peak_snp_bed = peak_snp_bed
        self.ref_fasta = ref_fasta
        self.result_path = result_path

    # read in the bed
    def read_bed(self):
        combined = pd.read_csv(self.peak_snp_bed, sep="\t", header=None)
        # snp_bed_df = pd.read_csv(self.snp_bed, sep="\t", header=None)
        # combined = pd.concat([peak_bed_df, snp_bed_df], axis=1)

        combined.columns = ['peak_chr', 'peak_start', 'peak_end',
                            'snp_chr', 'snp_start', 'snp_end',
                            'snp_id', 'Effect_Allele', 'Allele']

        # check if the the snp is within peak
        if combined['peak_chr'].equals(combined['snp_chr']) and \
                all(combined['peak_start'] <= combined['snp_start'])and \
                all(combined['peak_end'] >= combined['snp_end']):
            pass
        else:
            raise ValueError('SNP is not within peak')

        return combined

    # compute the location index of snp within the peak
    def compute_pos_index(self, combined):
        return list(combined['snp_start'] - combined['peak_start'])

    # generate wt and mutant fasta
    # TODO check why the snp info and the returned nucliec does not match - returned one pos further
    def get_fasta(self, snp_pos, combined):
        peak_bed_bt = BedTool(self.peak_snp_bed)
        ref_bt = BedTool(self.ref_fasta)
        peak_fa = peak_bed_bt.sequence(fi=ref_bt)
        peak_fa = open(peak_fa.seqfn).read()
        peak_fa_list = self._process_seq(peak_fa, combined)
        mutant_peak_fa_list = self._generate_mutant_fa(combined, peak_fa_list)

        return np.array(peak_fa_list), np.array(mutant_peak_fa_list)

    def _process_seq(self, peak_fa, combined):
        n_peak = combined.shape[0]
        index_list = [i*2+1 for i in range(n_peak)]
        peak = peak_fa.split("\n")
        peak_list = [peak[i].upper() for i in index_list] # TODO: change to upper
        #print(peak_list)
        return peak_list

    def _generate_mutant_fa(self, combined, peak_fa_list):

        snp_index = combined['snp_start'] - combined['peak_start']
        snp_index_list = []
        for i in snp_index:
            snp_index_list.append(i)

        mutant_peak_fa_list = []
        # TODO report which one is effect, generate a new file
        # 0 means WT is effect allele, 1 means mutant is effect allele,
        # 2 means Given WT nucleic acid does not match with the one retrieve in fasta'
        effect_allele_index = []

        for i in range(len(peak_fa_list)):
            fa = peak_fa_list[i]
            snp_index = snp_index_list[i]

            ## WT is Effect_Allele
            if fa[snp_index] == combined.iloc[i, 7]:
                fa = list(fa)
                fa[snp_index] = combined.iloc[i, 8].upper()
                mutant_fa = ''.join(fa)
                mutant_peak_fa_list.append(mutant_fa)
                effect_allele_index.append(0)
            ## WT is Noneffect_Allele
            elif fa[snp_index] == combined.iloc[i, 8]:
                fa = list(fa)
                fa[snp_index] = combined.iloc[i, 7].upper()
                mutant_fa = ''.join(fa)
                mutant_peak_fa_list.append(mutant_fa)
                effect_allele_index.append(1)

            else:
                print(fa)
                print(fa[snp_index])
                print(combined.iloc[i, :])
                effect_allele_index.append(2)

                raise ValueError('Given WT nucleic acid does not match with the one retrieve in fasta')

        pd.DataFrame(np.array(effect_allele_index)).to_csv(
            os.path.join(self.result_path, "effect_allele_index.csv"), header = False, index = True
        )
        return mutant_peak_fa_list