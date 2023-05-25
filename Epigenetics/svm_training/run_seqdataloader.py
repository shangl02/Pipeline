import os
import pybedtools
from splits import *
# from labelgen.__init__ import *

import sys
from seqdataloader import *
from  seqdataloader.labelgen import *

cluster = 1 #int(sys.argv[1])
split = int(sys.argv[1])

sc_tasks = '/hpc/grid/wip_drm_targetsciences/users/shangzhong/epigenetics/gkm_svm/seqdataloader/tasks.tsv'
test_chroms = hg19_splits[split]['test']

train_params={
    'task_list':sc_tasks,
    'outf':'/hpc/grid/wip_drm_targetsciences/users/shangzhong/epigenetics/gkm_svm/model/fold'+str(split)+'/train.inputs.bed.gz',
    'output_type':'gzip',
    'chrom_sizes':'/hpc/grid/wip_drm_targetsciences/users/shangzhong/epigenetics/gkm_svm/seqdataloader/hg38.chrom.sizes.txt',
    'bin_stride':50,
    'left_flank':400,
    'right_flank':400,
    'bin_size':200,
    'task_threads':4,
    'chroms_to_exclude':test_chroms,
    'chrom_threads':4,
    'labeling_approach':'peak_summit_in_bin_classification',
    'store_positives_only':False,
    'allow_ambiguous':True,
    }

test_params={
    'task_list':sc_tasks,
    'outf':'/hpc/grid/wip_drm_targetsciences/users/shangzhong/epigenetics/gkm_svm/model/fold'+str(split)+'/test.inputs.bed.gz',
    'output_type':'gzip',
    'chrom_sizes':'/hpc/grid/wip_drm_targetsciences/users/shangzhong/epigenetics/gkm_svm/seqdataloader/hg38.chrom.sizes.txt',
    'bin_stride':50,
    'left_flank':400,
    'right_flank':400,
    'bin_size':200,
    'task_threads':4,
    'chroms_to_keep':test_chroms,
    'chrom_threads':4,
    'labeling_approach':'peak_summit_in_bin_classification',
    'store_positives_only':False,
    'allow_ambiguous':True,
    }

genomewide_labels(train_params)
genomewide_labels(test_params)
