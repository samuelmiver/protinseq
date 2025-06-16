#!/usr/bin/env python
import glob
import os, sys
def create_job(pcr_dup=1):
    fils = glob.glob('/users/lserrano/sequencing_data/Maria_Lluch/2019-02-15_AHWLNMBCX2/C_0_5*_read1.fastq.gz')
    fils += glob.glob('/users/lserrano/sequencing_data/Maria_Lluch/2019-02-15_AHWLNMBCX2/C_2*_read1.fastq.gz')
    fils += glob.glob('/users/lserrano/sequencing_data/Maria_Lluch/2019-02-15_AHWLNMBCX2/D*_read1.fastq.gz')
    tn_seqs = {'A':'CGAGGGGGGGCCCTTTTACACAATTATACGGACTTAATC', 'B':'TGATTTTTTTCTCTTTTACACAATTATACGGACTTAATC',
               'C':'CGAGGGGGGGCCCTTTTACACAGTTGTACGGACTTAATC', 'D':'TGATTTTTTTCTCTTTTACACAGTTGTACGGACTTAATC'}
    for fil in fils:
        ide = fil.split('/')[-1][:6]
        subsh = '/users/lserrano/smiravet/proteoseq/sub_'+ide+'.sh'
        r1 = fil
        r2 = fil.replace('_read1', '_read2')
        tn = tn_seqs[fil.split('users/lserrano/sequencing_data/Maria_Lluch/2019-02-15_AHWLNMBCX2/')[1][0]]
        
        if pcr_dup:
            director = '/users/lserrano/smiravet/proteoseq/datasets/pcr_duplicates_removed'
            pd = '1'
        else:
            director = '/users/lserrano/smiravet/proteoseq/datasets/with_pcr_duplicates'
            pd = '0'

        # cmd = "python /software/ls/fastqins/fastqins_web.py -i1 "+r1+" -i2 "+r2+" -t "+tn+" -v --rm_inter 0 -p 1 -o /users/lserrano/smiravet/proteoseq/datasets/real_dup_rm\n"
        # cmd = "python /software/ls/fastqins/fastqins_web.py -i1 "+r1+" -i2 "+r2+" -t "+tn+" -v --rm_inter 0 -p 0 -o /users/lserrano/smiravet/proteoseq/datasets/with_pcr_duplicates\n"
        cmd = "python ./protinseq_fastqins.py -i1 "+r1+" -i2 "+r2+" -t "+tn+" -v --rm_inter 0 -p "+pd+" -o "+director
        print fil+" being processed."
        os.system(cmd)

if __name__=='__main__':
    if sys.argv[1]=='0':
        print 'PCR duplicates not being removed'
        create_job(int(sys.argv[1]))
    elif sys.argv[1]=='1':
        print 'PCR duplicates being removed'
        create_job(int(sys.argv[1]))
    else:
        sys.exit('provide argument')
