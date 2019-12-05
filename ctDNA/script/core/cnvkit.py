#Email:fanyucai1@126.com
#2019.7.11

import subprocess
import argparse
import sys
import glob
#http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refFlat.txt.gz
anno = "/software/CNVkit/cnvkit-master/data/refFlat.txt"
access = "/software/CNVkit/cnvkit-master/data/access-5k-mappable.hg19.bed"
python3 = "/software/python3/Python-v3.7.0/bin/python3"
cnvkit = "/software/python3/Python-v3.7.0/bin/cnvkit.py"
ref = "/data/Database/hg19/ucsc.hg19.fasta"
def run(tumor, normal, bed, outdir):
    cmd="%s %s batch %s --normal %s/*.bam --targets %s --fasta %s --output-reference my_reference.cnn --output-dir %s --annotate %s --access %s" %(python3,cnvkit,tumor,normal,bed,ref,outdir,anno,access)
    print(cmd)
    subprocess.check_call(cmd,shell=True)
    cns=glob.glob("%s/*.cns"%(outdir))
    infile=open(cns[0],"r")
    outfile=open("%s/cnv.final.tsv"%(outdir),"w")
    outfile.write("#Chr\tStart\tend\tgene\tlog2\ttype\tCopy\n")
    for line in infile:
        if not line.startswith("chromosome"):
            line=line.strip()
            array=line.split("\t")
            copy = 2 ** float(array[4]) * 2
            if float(array[4])>=0.585:#https://cnvkit.readthedocs.io/en/stable/calling.html
                type="gain"
                tmp = array[0] + "\t" + array[1] + "\t" + array[2] + "\t" + array[3] + "\t" + array[4] + "\t" + type + "\t" + str(copy)
                outfile.write("%s\n" % (tmp))
            if float(array[4]) <=-1:
                type="loss"
                tmp=array[0]+"\t"+array[1]+"\t"+array[2]+"\t"+array[3]+"\t"+array[4]+"\t"+type+"\t"+str(copy)
                outfile.write("%s\n"%(tmp))
    outfile.close()

if __name__=="__main__":
    if len(sys.argv)<5:
        print("python3 %s tumor_bam normal_bam bedfile outdir anno access python3 cnvkit ref\n"%(sys.argv[0]))
        print("Email:fanyucai1@126.com")
    else:
        tumor = sys.argv[1]
        normal = sys.argv[2]
        bed = sys.argv[3]
        outdir = sys.argv[4]
        run(tumor, normal, bed, outdir)