import glob, os
import argparse
import pandas as pd
import gzip
from shutil import copyfile

data_path = '/hpc/grid/hgcb/workspace/projects/P049_Schizophrenia_Sequencing/data'
target_path = '/hpc/grid/hgcb/workspace/users/shangzhong/GATK/schizo'
param_file = '/home/lis262/Code/Pipeline/GATK/p01_GATK_Parameters.config'
ppl_fn = '/home/lis262/Code/Pipeline/GATK/p01_GATK_get_gvcf.nf'
# 1. find samples that are sequenced
fq_path = data_path + '/fastq/MY1610072_R1/fastq_Samples_1-50'
fq_path1 = data_path + '/fastq/MY1610072_R1/fastq_Samples_51-61'

fq_files = glob.glob(fq_path + '/*_R1_*.fastq.gz') + glob.glob(fq_path1 + '/*_R1_*.fastq.gz')
samples = ['PS-' + f.split('/')[-1].split('_')[0] for f in fq_files]

# 2. build dictionary {sample: family_personID}
def get_dict(covari_fn, samples):
	covari_df = pd.read_csv(covari_fn)
	covari_df = covari_df[covari_df['PS#'].isin(samples)]
	covari_df['NInd'] = covari_df['NInd'].astype('str')
	covari_df['Mom'] = covari_df['Mom'].astype('str')
	covari_df['Dad'] = covari_df['Dad'].astype('str')
	# change the id of the individuales
	covari_df['NInd'] = covari_df['Pedi#'] + '_' + covari_df['NInd']
	covari_df['Mom'] = covari_df['Pedi#'] + '_' + covari_df['Mom']
	covari_df['Dad'] = covari_df['Pedi#'] + '_' + covari_df['Dad']
	# build dictionary
	sample_id_dict = covari_df.set_index('PS#')['NInd'].to_dict()
	return sample_id_dict

covari_fn = data_path + '/scz_phe_covar.csv'
sample_id_dict = get_dict(covari_fn, samples)

# 3. copy file and change file names, and submit jobs.
from shutil import copyfile
for f in fq_files[45:61]:
    fq_id = 'PS-' + f.split('/')[-1].split('_')[0]
    fam_id = sample_id_dict[fq_id]
    if fam_id == 'AV-27_3': continue
    fam = fam_id.split('_')[-1]
    des_path = target_path + '/' + fam_id
    if not os.path.exists(des_path):
        os.mkdir(des_path)
    des_fn = des_path + '/' + fam_id + '_1.fq.gz'
    copyfile(f, des_fn)
    des_fn = des_path + '/' + fam_id + '_2.fq.gz'
    ori2 = f.replace('_R1_','_R2_')
    copyfile(ori2, des_fn)
    # change the parameter file
    des_param = des_path + '/p01_GATK_Parameters.config'
    cmd = ('sed \"/fq_path/c\\\tfq_path = \'{path}\'\" {old_param} > '
           '{target}').format(path=des_path, old_param=param_file,
                             target=des_param)
    os.system(cmd)
    # submit jobs
    os.chdir(des_path)
    cmd = ('bsub -n 8 -M 31457280 -R \"span[ptile=8]\" '
          '-o log.txt -q long -e log.err.txt '
          '\"nextflow run {ppl} -c {config} -resume\"').format(
            ppl=ppl_fn,config=des_param)
    print('run',fam_id)
    print(cmd)
    os.system(cmd)
