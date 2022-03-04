import os, glob, gzip
import argparse
import pandas as pd
import tabix


parser = argparse.ArgumentParser(description='extract significant QTL results and annotate them')


parser.add_argument('-p','--path',action='store',dest='path',help='work path')
parser.add_argument('--pval',type=float,action='store',dest='pval',help='pvalue threshold',default=5e-8)
parser.add_argument('-a','--anno',action='store',dest='anno',help='annotate results',default='/lustre/workspace/projects/BLCA/Results/snp_indel/anno.vep.tsv.gz')

args      = parser.parse_args()
work_path = args.path
pval      = args.pval
anno_fn   = args.anno

clump_path = f'{work_path}/matrixQTL/clump'
clump_fns  = sorted(glob.glob(f'{clump_path}/*.clumped'))

tb = tabix.open(anno_fn)
with gzip.open(anno_fn, 'rt') as in_f:
    columns = in_f.readline().strip().split('\t')

def anno_clump(row, tb, columns):
    '''tb: the object of tabix
    columns: header of the annotation file
    '''
    chrom, pos, ref, alt = row['ID'].split(':')
    genes = []
    conse = []
    cadd = ''
    gnomad = ''
    rsid = ''
    nearest = ''
    records = tb.query(f'chr{chrom}',int(pos),int(pos)+1)

    for record in records:
        try:
            dic = {k:v for k,v in zip(columns,record)}
        except:
            continue
        if (dic['ref'] == ref and dic['alt'] == alt) or \
           (dic['ref'] == alt and dic['alt'] == ref):
            cadd = dic['CADD_PHRED']
            gnomad = dic['gnomAD_AF']
            # some rsids have multiple other ids. eg: rs138213197&CM120253&COSV51705978
            rsid = [r for r in dic['Existing_variation'].split('&') if r.startswith('rs')]
            if rsid == []:
                rsid = ''
            else:
                rsid = rsid[0]
            nearest = dic['near_50Kgene']
            genes.append(dic['SYMBOL'])
            conse.append(dic['Consequence'])
    return pd.Series([','.join(genes),nearest,','.join(conse),rsid,gnomad,cadd])

# generates {rsid : chr:pos:ref:alt}
af_fn = f'{work_path}/plink/germ_newSp38.afreq'
af_df = pd.read_csv(af_fn, sep='\t', header=0)
af_df['#CHROM'] = af_df['#CHROM'].astype('str')
af_df['POS'] = af_df['POS'].astype('str')

af_df['snp'] = af_df['#CHROM']+':'+af_df['POS']+':' + \
               af_df['REF']+':'+af_df['ALT']
rsid_snp_dic = af_df.set_index('ID')['snp'].to_dict()

anno_clump_fn = f'{work_path}/matrixQTL/clump_sig_anno.tsv'
head = True
for clump in clump_fns:
    clump_df = pd.read_csv(clump, delim_whitespace=True, header=0)
    clump_df['cell'] = '.'.join(clump.split('/')[-1].split('.')[:-1])
    clump_df['ID'] = clump_df['SNP'].map(lambda x: rsid_snp_dic[x])
    clump_df[['Gene','near_50kgene','Consequence','rsid','gnomad','CADD']] = \
         clump_df.apply(lambda row:anno_clump(row, tb, columns), axis=1)
    if head:
        clump_df.to_csv(anno_clump_fn,sep='\t',index=False)
        head = False
    else:
        clump_df.to_csv(anno_clump_fn,sep='\t',index=False,header=False,mode='a')
    
