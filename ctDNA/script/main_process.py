import os
import sys
sub=os.path.abspath(__file__)
dir_name=os.path.dirname(sub)
sys.path.append(dir_name)
import core
import subprocess
outdir="/result/"
def run(pe1,pe2,bed,read_length,prefix,vaf):
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    ###################step1##################pre_UMI
    bam_out=outdir+"/"+"bam"
    if not os.path.exists(bam_out):
        os.mkdir(bam_out)
    if not os.path.exists("%s/bam/%s.consensus.filter.clipped.bam"%(outdir,prefix)):
        core.pre_UMI.run(pe1,pe2,bam_out,prefix,read_length)
    ###################step2##################call_SNV
    snv_out=outdir+"/"+"snv"
    if not os.path.exists(snv_out):
        os.mkdir(snv_out)
    if not os.path.exists("%s/%s.vardict.vcf"%(snv_out,prefix)):
        core.vardict_single.run(vaf,"%s/bam/%s.consensus.filter.clipped.bam"%(outdir,prefix),bed,prefix,snv_out)
    ###################step3###################annovar and filter
    anno_out = outdir + "/" + "anno"
    if not os.path.exists(anno_out):
        os.mkdir(anno_out)
    if not os.path.exists("%s/%s.annovar.filter.tsv"%(anno_out,prefix)):
        core.format_vcf.run(prefix,"%s/%s.vardict.vcf"%(snv_out,prefix),anno_out)
        core.anno_vcf.run("%s/%s.vcf"%(anno_out,prefix),anno_out,prefix)
        core.filter_vcf.run(0.01,"%s/%s.annovar.tsv"%(anno_out,prefix),anno_out,prefix)
    ##################step4#####################genefuse
    fuse_out = outdir + "/" + "genefuse"
    if not os.path.exists(fuse_out):
        os.mkdir(fuse_out)
    if not os.path.exists("%s/%s.tsv"%(fuse_out,prefix)):
        core.genefuse.run(pe1,pe2,prefix,fuse_out)
    #################step5######################call CNV
    cnv_out = outdir + "/" + "cnv"
    if not os.path.exists(cnv_out):
        os.mkdir(cnv_out)
    if not os.path.exists("%s/cnv.final.tsv"%(cnv_out)):
        core.cnvkit.run("%s/bam/%s.consensus.filter.clipped.bam"%(outdir,prefix),"/data/Database/CNV_baseline/27ctDNA/",bed,cnv_out)
    #####################step6########################
    stat_out = outdir + "/" + "qc"
    if not os.path.exists(stat_out):
        os.mkdir(stat_out)
    core.sample_qc.run(bed, bed, "%s/bam/%s.consensus.filter.clipped.bam" % (outdir, prefix), stat_out, prefix)
    #################step7#######################delete tmp file
    for (root,dirs,files) in os.walk(outdir):
        for file in files:
            tmp=os.path.join(root,file)
            if not tmp.endswith("filter.clipped.bam"):#####retain final bam
                if not tmp.endswith("filter.clipped.bai"):
                    if not tmp.endswith("tsv"):####retain snv and cnv final result
                        if not tmp.endswith("result") and not tmp.endswith("vcf") and not tmp.endswith("bed"):
                            subprocess.check_call("rm -rf %s"%(tmp),shell=True)

if __name__=="__main__":
    if len(sys.argv)<6:
        print("python3 %s pe1.fq pe2.fq bedfile readlength prefix vaf"%(sys.argv[0]))
        print("Email:fanyucai1@126.com")
    else:
        pe1 = sys.argv[1]
        pe2 = sys.argv[2]
        bed = sys.argv[3]
        read_length = sys.argv[4]
        prefix = sys.argv[5]
        vaf=0.001
        if len(sys.argv)==7:
            vaf=sys.argv[6]
        run(pe1, pe2, "%s/%s"%(outdir,bed), read_length, prefix, vaf)