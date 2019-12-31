import os,sys
import pysam
insam=sys.argv[1]
outsam=sys.argv[2]
#m_path='/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target'

mapq_threshold=0
n_primary=0
n_not_xa=0
n_not_xa_mapq_gt_20=0
bam_f=insam#os.path.join(m_path,'2-0-0BWA_map',insam)
out_bam=outsam#os.path.join(m_path,'2-0-0BWA_map',outsam)
bam = pysam.AlignmentFile(bam_f, 'rb')
bam_unique_mis6=pysam.AlignmentFile(out_bam, "w", template=bam)
for read in bam:
    if read.is_secondary:  # not the primary alignment
        continue
    n_primary+=1
    #if read.has_tag('MD'):break
    if read.has_tag('XA'):continue
    n_not_xa+=1
    if read.mapping_quality<=mapq_threshold:continue
    n_not_xa_mapq_gt_20+=1
    md_tag = read.get_tag('MD').split()[0]
    flag=0
    mis=0
    for ab in md_tag :
        if ab =='^':
            flag=1
        elif ab in [str(x)for x in range(10)]:flag=0
        elif ab.lower() in  [chr(x) for x in range(97,123)]:
            if flag==0:
                mis+=1
    if mis>6:continue
    bam_unique_mis6.write(read)
bam_unique_mis6.close()
bam.close()
print "n_primary:",n_primary
print "n_not_xa:",n_not_xa
print "n_not_xa_mapq_gt_20:",n_not_xa_mapq_gt_20   

