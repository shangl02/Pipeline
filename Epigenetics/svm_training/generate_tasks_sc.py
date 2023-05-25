import pandas as pd


bed_fn = '/hpc/grid/wip_drm_targetsciences/users/shangzhong/epigenetics/gkm_svm/Alveolar_Cap_Endo.bed.gz'

sc_tasks = '/hpc/grid/wip_drm_targetsciences/users/shangzhong/epigenetics/gkm_svm/seqdataloader/tasks.tsv'
task_dict = {}

for task in range(1,2):
    task = 'Cluster' + str(task)
    idr_peaks = bed_fn
    # idr_peaks = basedir_mnt + 'idr_peaks/' + task + '.idr.optimal.narrowPeak'
    # bigwigs = basedir_oak + task + '/signal/rep1/' + task + '.fc.signal.bigwig'
    # ambig_peaks = basedir_mnt + 'ambiguous/' + task + '.ambiguous.bed'
    # task_dict[task] = [idr_peaks, bigwigs, ambig_peaks]
    task_dict[task] = [idr_peaks]

outf = open(sc_tasks, 'w')
for task in task_dict:
    outf.write(task + '\t' + '\t'.join(task_dict[task]) + '\n')

