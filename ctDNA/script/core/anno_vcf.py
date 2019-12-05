#Email:fanyucai1@126.com
#2019.7.5

import os
import subprocess
import re
import sys
Canonical_transcript_file="/data/Database/knownCanonical/clinvar_canonical_trans.txt"
out_name=['Chr','Start','End','Ref','Alt','Func.refGene','Gene.refGene','GeneDetail.refGene',
          'ExonicFunc.refGene',	'AAChange.refGene',	'cytoBand',	'1000g2015aug_all',	'avsnp150',	'snp138',
          'CLNALLELEID','CLNDN','CLNDISDB',	'CLNREVSTAT','CLNSIG','cosmic88_coding','SIFT_score','SIFT_pred',
          'Polyphen2_HDIV_score','Polyphen2_HDIV_pred','esp6500siv2_all','ExAC_ALL','ExAC_EAS','1000g2015aug_eas',
          '1000g2015aug_sas','1000g2015aug_afr','1000g2015aug_amr','1000g2015aug_eur','InterVar_automated','GT','AAChange.1',
          'Ref_Reads',	'Alt_Reads','Var']
annovar="/data/Database/Annovar/"
def run_hgvs(var_site):
    p1=re.search(r'[A-Z]fs\*\d+$',var_site)###匹配移码突变
    p2=re.search(r'del([ACGT]+)ins',var_site)###匹配del和ins
    if p1:
        new=re.sub(r'[A-Z]fs\*\d+$',"",var_site)
        new=new+"fs"
    else:
        new=var_site
    if var_site.endswith("X"):####终止密码子X替换*
        new1= re.sub(r'X$', "*", new)
    else:
        new1=new
    if p2:
        new2=re.sub(p2.group(1),"",new1,count=1)
    else:
        new2 = new1
    return new2
def run(vcf,outdir,prefix):
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    out=outdir+"/"+prefix
    #########################get Canonical transcript info
    transcript={}
    infile=open(Canonical_transcript_file,"r")
    for line in infile:
        line=line.strip()
        array=line.split("\t")
        transcript[array[0]]=[]
        for j in range(len(array)):
            if j!=0:
                tmp=array[j].split(".")
                transcript[array[0]].append(tmp[0])
    infile.close()
    ##########################run annovar
    par = " -protocol refGene,cytoBand,snp138,avsnp150,exac03,esp6500siv2_all,1000g2015aug_all,1000g2015aug_eas,gnomad211_exome,gnomad211_genome,cosmic88_coding,clinvar_20190305,ljb26_all,intervar_20180118"
    par += ",1000g2015aug_sas,1000g2015aug_afr,1000g2015aug_amr,1000g2015aug_eur "
    par += " -operation g,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f "
    par += " -nastring . -polish "
    subprocess.check_call("perl %s/table_annovar.pl %s %s/humandb -buildver hg19 -out %s -remove %s -vcfinput " % (annovar, vcf, annovar, out, par), shell=True)
    subprocess.check_call("rm -rf %s.hg19_multianno.vcf %s.avinput" % (out, out), shell=True)
    infile = open("%s.hg19_multianno.txt" % (out), "r")
    outfile = open("%s.annovar.tsv" % (out), "w")
    for i in range(len(out_name)):
        if i == 0:
            outfile.write("%s" % (out_name[i]))
        else:
            outfile.write("\t%s" % (out_name[i]))
    outfile.write("\n")
    dict = {}
    for line in infile:
        line = line.strip()
        array = line.split("\t")
        name = []
        if line.startswith("Chr"):
            for i in range(len(array)):
                name.append(array[i])
                dict[array[i]] = i
        else:
            p1 = re.compile(r'Alt_Reads=([0-9]+)')
            p2 = re.compile(r'Ref_Reads=([0-9]+)')
            p3 = re.compile(r'Var=([0-9.]+)')
            p4 = re.compile(r'GT=(\d+/\d+)')
            Alt_Reads = p1.findall(line)
            Ref_Reads = p2.findall(line)
            Var = p3.findall(line)
            GT = p4.findall(line)
            ##########################format output knownCanonical transcript
            tmp = array[dict['AAChange.refGene']].split(",")
            final_nm = ""
            if array[6] in transcript:
                for i in transcript[array[6]]:
                    if final_nm == "":
                        for k in tmp:
                            if final_nm=="" and re.search(i,k):
                                final_nm = k
                                continue
                            else:
                                pass
                    else:
                        continue
            if final_nm=="":
                final_nm=tmp[0]
            for l in range(len(out_name)):
                if l == 0:
                    outfile.write("%s" % (array[dict[out_name[l]]]))
                elif out_name[l] == "Var":
                    tmp_num = float(Var[0]) * 100
                    outfile.write("\t%.2f" % (tmp_num) + "%")
                elif out_name[l] == "Alt_Reads":
                    outfile.write("\t%s" % (Alt_Reads[0]))
                elif out_name[l] == "Ref_Reads":
                    outfile.write("\t%s" % (Ref_Reads[0]))
                elif out_name[l] == "AAChange.1":
                    outfile.write("\t%s" % (run_hgvs(final_nm)))
                elif out_name[l] == "GT":
                    outfile.write("\t%s" % (GT[0]))
                else:
                    outfile.write("\t%s" % (array[dict[out_name[l]]]))
            outfile.write("\n")
    infile.close()
    outfile.close()
    if os.path.exists("%s.hg19_multianno.txt" % (out)):
        subprocess.check_call("rm -rf %s.hg19_multianno.txt" % (out), shell=True)
    ###########################################################

if __name__=="__main__":
    if len(sys.argv)!=4:
        print("python3 %s format.vcf outdir prefix"%(sys.argv[0]))
        print("\nEmail:fanyucai1@126.com")
    else:
        vcf=sys.argv[1]
        outdir=sys.argv[2]
        prefix = sys.argv[3]
        run(vcf,outdir,prefix)