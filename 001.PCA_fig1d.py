import sys,os
sys.path.insert(1, "/home/zyang/software/MitEdit")
from basic import Basic

def getlabels(dir_path,datagroup):
    ##由于代码字符大小写出错了，所以对input的筛选失败了，在R中筛�?
    files = os.listdir(dir_path)
    bam_file = []
    labels = []
    for i in files:
        if "input" in i:
            continue
        if i.endswith("bam"):
            bam_filepath = dir_path+"/"+i
            if datagroup=="YAP":
                head_type = i.strip().split("_")[0]
                marker = i.strip().split("_")[2].split(".")[0]
                sample = i.strip().split("_")[1]
                label = head_type+"_"+marker+"_"+sample+"_"+datagroup
            else:
                head_type = "Cellline"
                marker=i.strip().split("_")[0]
                sample = i.strip().split("_",1)[1].split(".")[0]
                label = head_type + "_" + marker + "_" + sample + "_" + datagroup
            bam_file.append(bam_filepath)
            labels.append(label)
    return bam_file, labels




def multiBamSummary():


    CCLE_paired_path = "/home/zhluo/Project/TF_enrichment/bowtie_mapping_paired"
    CCLE_single_path = "/home/zhluo/Project/TF_enrichment/bowtie_mapping_single"
    YAP_bam_path = "/home/zhluo/Project/CRC/data_nazhang/step56_human_CRC/rename_bam"

    CCLE_paired_bam_file, CCLE_paired_labels = getlabels(CCLE_paired_path, "CCLE_paired")
    CCLE_single_bam_file, CCLE_single_labels = getlabels(CCLE_single_path, "CCLE_single")
    YAP_bam_file, YAP_labels = getlabels(YAP_bam_path, "YAP")
    ###multiBamSummary
    #datagroup = "CCLE_paired"
    #bam_file,labels=getlabels(datapath,datagroup)
    cmd = "multiBamSummary bins --binSize 100000 --bamfiles "+ " ".join(CCLE_paired_bam_file+CCLE_single_bam_file+YAP_bam_file) + " "
    cmd = cmd + "--minMappingQuality 30 -p 24 "
    cmd = cmd + "--labels " + " ".join(CCLE_paired_labels+CCLE_single_labels+YAP_labels) + " "
    cmd = cmd + "-out "+"readCounts.npz --outRawCounts " +"readCounts.tab"
    pbs_file = "/home/mwshi/project/CRC_enhancer/bincount/multiBamSummary.pbs"
    pbs = open(pbs_file,"w")
    pbs.write("cd /home/mwshi/project/CRC_enhancer/bincount/\n")
    pbs.write(cmd+"\n")
    pbs.close()
    cmd = "qsub -q batch -l nodes=1:ppn=24 %s" % pbs_file
    Basic.run(cmd, wkdir="/home/mwshi/project/CRC_enhancer/bincount/")


if __name__ == "__main__":
    ##step 1
    multiBamSummary("/home/zhluo/Project/TF_enrichment/bowtie_mapping_paired", "CCLE_paired")
    multiBamSummary("/home/zhluo/Project/TF_enrichment/bowtie_mapping_single", "CCLE_single")
    multiBamSummary("/home/zhluo/Project/CRC/data_nazhang/step56_human_CRC/rename_bam", "YAP")
