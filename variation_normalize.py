import os,sys,json,re,gzip
import argparse
from pyfaidx import Fasta

###
#1. Left-align and normalize indels;
#2. check if REF alleles match the reference;
#3. split multiallelic sites into multiple rows;
#by wangtaifu 2022.7.6
###

ap = argparse.ArgumentParser()
ap.add_argument(
    '-ref',
    '--reference',
    type=str,
    default=None,
    required=True,
    metavar='FASTA',
    help='The genome sequencing'
)
ap.add_argument(
    '-f',
    '--input_file',
    type=str,
    default=None,
    required=True,
    metavar='info',
    help='the path of input file, should be Tab-separated file.'
)
ap.add_argument(
    '-c',
    '--cpra',
    type=str,
    default="1,2,4,5",
    required=True,
    metavar='column number',
    help='the column number of chr-pos-ref-alt(CPRA), for example: 1,2,4,5'
)
# ap.add_argument(
#     '-s',
#     '--split',
#     type=str,
#     default="yes",
#     required=True,
#     metavar='yes or no',
#     help='yes: split multiallelic sites into multiple rows. no: do not split'
# )
ap.add_argument(
    '-o',
    '--outpre',
    type=str,
    default="./var_norm",
    required=True,
    metavar='outpre',
    help='filename of output'
)
args = ap.parse_args()

def has_non_acgtn(seq):
    for c in seq.upper():
        if  c!='A' and c!='C' and c!='G' and c!='T' and c!='N' :
            return 1;
    return 0

def check_SymOrRef(fasta,chrom,pos,ref,alt):
    '''
ERR_SYMBOLIC
ERR_SPANNING_DELETION
ERR_REF_MISMATCH
    '''
    ERR_SYMBOLIC=1
    ERR_SPANNING_DELETION=2
    ERR_REF_MISMATCH=3
    if "<" in alt:
        sys.stderr.write("symbolic allele at {}:{}; alt: {} ".format(chrom,pos,alt))
        return ERR_SYMBOLIC
    if "*" in alt:
        sys.stderr.write("spanning deletion at {}:{}; alt: {} ".format(chrom,pos,alt))
        return ERR_SPANNING_DELETION
    try:
        RefInFasta=Fasta(fasta)[chrom][int(pos)-1:int(pos)+len(ref)-1].seq
    except:
        sys.stderr.write("The chromosome coordinates of VCF were inconsistent with fasta. {}:{}-{} in VCF \n".format(chrom,pos,str(int(pos)+len(ref))))
        return ERR_REF_MISMATCH
    if has_non_acgtn(ref):
        sys.stderr.write("Non-ACGTN reference allele at {}:{}; REF_SEQ:{} vs VCF:{}\n".format(chrom,pos,RefInFasta,ref))
        return ERR_REF_MISMATCH
    if RefInFasta.upper() != ref.upper():
        sys.stderr.write("Reference allele mismatch at {}:{}; REF_SEQ:{} vs VCF:{}\n".format(chrom,pos,RefInFasta,ref))
        return ERR_REF_MISMATCH
    return 0

def right_trim(fasta,chrom,pos,ref,alt):

    if ref[-1].upper() != alt[-1].upper():
        return chrom,pos,ref,alt

    #trim from right
    min_len=min(len(ref),len(alt))
    for i in range(1,min_len+1):
        flag=0
        if ref[len(ref)-i].upper() != alt[len(alt)-i].upper() : break
        flag=1

    if flag == 0 :
        new_ref=ref[0:len(ref)-i+1]
        new_alt=alt[0:len(alt)-i+1]
        return chrom,pos,new_ref,new_alt
    else:  
        extend_seq=Fasta(fasta)[chrom][int(pos)-2].seq   # extend from fasta
        new_ref=ref[0:len(ref)-i]
        new_alt=alt[0:len(alt)-i]
        #print(str(i),new_ref,new_alt,ref,alt,chrom,new_pos,extend_seq)
        new_pos=int(pos)-1        
        new_ref=''.join([extend_seq,new_ref])
        new_alt=''.join([extend_seq,new_alt])
        #print(str(i),new_ref,new_alt,ref,alt,chrom,new_pos,extend_seq)
        return right_trim(fasta,chrom,new_pos,new_ref,new_alt)      

def left_trim(chrom,pos,ref,alt):    
    if ref[0].upper() != alt[0].upper() or (len(ref)==1 or len(alt)==1):
        return chrom,pos,ref,alt

    #trim from left
    min_len=min(len(ref),len(alt))
    for i in range(min_len-1):
        flag=0
        if ref[i].upper() != alt[i].upper() : break
        flag=1
        
    if flag == 0 :
        new_ref=ref[i:]
        new_alt=alt[i:]  
        new_pos=int(pos)+i  
        return chrom,new_pos,new_ref,new_alt
    else: 
        new_ref=ref[i+1:]
        new_alt=alt[i+1:]         
        new_pos=int(pos)+i+1
        return chrom,new_pos,new_ref,new_alt

if __name__ == '__main__':
    outpre = args.outpre
    cpra = (args.cpra).split(',')
    #split = args.split
    fasta = args.reference
    inputf=args.input_file
    if '.gz' == inputf[-3:]:
        infile=gzip.open(inputf, 'rt')
    else:
        infile=open(inputf, 'r')

    outfile=open(outpre, 'w')
    for line in infile:
        if re.match("^#",line ) : outfile.write(line); continue 
        field_list=line.split("\t")
        chrom,pos,ref,alt=[field_list[int(i)-1] for i in cpra]
        alleles=alt.split(',')        
        for allele in alleles:
            code=check_SymOrRef(fasta,chrom,pos,ref,allele)
            if code:  # check failed  
                chrom,new_pos,new_ref,new_alt=[chrom,pos,ref,allele]                 
                #continue;
            else: 
                chrom,new_pos,new_ref,new_alt=right_trim(fasta,chrom,pos,ref,allele)
                chrom,new_pos,new_ref,new_alt=left_trim(chrom,new_pos,new_ref,new_alt)
                # code=check_SymOrRef(fasta,chrom,new_pos,new_ref,new_alt)
            for i,v in enumerate(cpra):
                field_list[int(v)-1]=[chrom,str(new_pos),new_ref,new_alt][i]
            outstr="\t".join(field_list)
            outfile.write("{}".format(outstr))   
