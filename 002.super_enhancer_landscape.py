###super enhancer绘制全局图
import sys,os
import pandas as pd
import linecache
sys.path.insert(1, "/home/zyang/software/MitEdit")
from basic import Basic


def mkdir(path):
    folder = os.path.exists(path)
    if not folder:                   #判断是否存在文件夹如果不存在则创建为文件夹
        os.makedirs(path)            #makedirs 创建文件时如果路径不存在会创建这个路径
    else:
        print "---  There is this folder!  ---"


def se_files(dir_paths):
    SuperEnhancer_files = []
    for dir_path in dir_paths:
        files = os.listdir(dir_path)
        for i in files:
            if i.endswith("Gateway_SuperEnhancers.bed") and "H3K27Ac".lower() in i.lower():
                SuperEnhancer_files.append(os.path.join(dir_path,i))
    return(SuperEnhancer_files)
    

def merge_bed():
    dir_paths = ["/home/zhluo/Project/CRC/data_nazhang/step56_human_CRC/super_enhancer/",
                 "/home/zhluo/Project/TF_enrichment/super_enhancer/"]
    ##sefile,用来筛选返回H3K27Acmarker 以及Gateway_SuperEnhancers文件名
    SuperEnhancer_files = se_files(dir_paths)
    ##直接写文件做merge
    with open("/home/mwshi/project/TF_enrichment/super_enhancer/merge.bed"，"w") as f:
        for f0 in SuperEnhancer_files:
            with open(f0,"r") as f1:
                f.write(f1.read()+"\n")


###创建合适的软连接，方便取名
def getlabels(dir_path,datagroup):
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

def lnfunc(scrpath,datagroup,newpath):
    bam_file, labels = getlabels(scrpath, datagroup)
    mkdir(newpath)
    for i,j in zip(bam_file,labels):
        os.symlink(os.path.join(scrpath,i), os.path.join(newpath,j+".bam"))
        os.symlink(os.path.join(scrpath,i+".bai"), os.path.join(newpath,j+".bam.bai"))
        os.symlink(os.path.join(scrpath,i+".flagstat"), os.path.join(newpath,j+".bam.flagstat"))

def lnfiles():
    
    CCLE_paired_path = "/home/zhluo/Project/TF_enrichment/bowtie_mapping_paired"
    CCLE_single_path = "/home/zhluo/Project/TF_enrichment/bowtie_mapping_single"
    YAP_bam_path = "/home/zhluo/Project/CRC/data_nazhang/step56_human_CRC/rename_bam"
    newCCLE_paired_path = "/home/mwshi/project/CRC_enhancer/chip_seq/CCLE_paired"
    newCCLE_single_path = "/home/mwshi/project/CRC_enhancer/chip_seq/CCLE_single"
    newYAP_bam_path = "/home/mwshi/project/CRC_enhancer/chip_seq/YAP"
    CCLE_paired_bam_file, CCLE_paired_labels = getlabels(CCLE_paired_path, "CCLE_paired")
    lnfunc(CCLE_paired_path,newCCLE_paired_path,"CCLEpaired")
    lnfunc(CCLE_single_path,newCCLE_single_path,"CCLEsingle")
    lnfunc(YAP_bam_path,newYAP_bam_path,"YAP")
    
    
    
    
    
def getlabepair(dirpath):
    filepath = []
    label = []
    for i in os.listdir(dirpath):
        if i.endswith("bam"):
        ####筛选出CO和H3K27ac相关的样本
            if i.startswith("Cellline") and "CO" not in i:
                continue
            if "H3K27ac".lower() in i.lower():
                filepath.append(os.path.join(dirpath,i))
                label.append(os.path.splitext(i)[0])
    return filepath, label
    
    
def multiBamSummary():
    newCCLE_paired_path = "/home/mwshi/project/CRC_enhancer/chip_seq/CCLE_paired"
    newCCLE_single_path = "/home/mwshi/project/CRC_enhancer/chip_seq/CCLE_single"
    newYAP_bam_path = "/home/mwshi/project/CRC_enhancer/chip_seq/YAP"
    ###获得文件名和对应的label
    a,b=getlabepair(newCCLE_paired_path)
    c,d=getlabepair(newCCLE_single_path)
    e,f=getlabepair(newYAP_bam_path) 
    filepaths = a+c+e
    labels = b+d+f
    cmd = "multiBamSummary BED-file --BED /home/mwshi/project/CRC_enhancer/super_enhancer/overlap_merge.bed --bamfiles " +  " ".join(filepaths)
    cmd = cmd + " --minMappingQuality 30 -p 24 "
    cmd = cmd + "--labels " + " ".join(labels) + " "
    cmd = cmd + "-out /home/mwshi/project/CRC_enhancer/super_enhancer/CO_H3K27ac_readCounts.npz "
    cmd = cmd + "--outRawCounts /home/mwshi/project/CRC_enhancer/super_enhancer/CO_H3K27ac_readCounts.tab"
    pbs_file = "/home/mwshi/project/CRC_enhancer/pbs/superenhancer_multiBamSummary.pbs"
    pbs = open(pbs_file,"w")
    pbs.write("cd /home/mwshi/project/CRC_enhancer/super_enhancer/\n")
    pbs.write(cmd+"\n")
    pbs.close()
    cmd = "qsub -q batch -l nodes=1:ppn=24 %s" % pbs_file
    Basic.run(cmd, wkdir="/home/mwshi/project/CRC_enhancer/pbs/")   

def getsumcount(columns):
    n = {}
    for i in columns:
        if i.endswith("CCLEpaired"):
            ff = r"/home/mwshi/project/CRC_enhancer/chip_seq/CCLE_paired/"+i+".bam.flagstat"
            with open(ff) as myfile:
                n[i]=int([next(myfile) for x in range(5)][4].split()[0])
        if i.endswith("CCLEsingle"):
            ff = r"/home/mwshi/project/CRC_enhancer/chip_seq/CCLE_single/"+i+".bam.flagstat"
            with open(ff) as myfile:
                n[i]=int([next(myfile) for x in range(5)][4].split()[0])               
        if i.endswith("YAP"):
            ff = r"/home/mwshi/project/CRC_enhancer/chip_seq/YAP/"+i+".bam.flagstat"
            with open(ff) as myfile:
                n[i]=int([next(myfile) for x in range(5)][4].split()[0])    
    return(n)

###计算rpm和tpm方便绘图
def calculate_rpm():
    df = pd.read_csv("/home/mwshi/project/CRC_enhancer/super_enhancer/CO_H3K27ac_readCounts.tab", sep="\t")
    df.columns = [x.strip("#").replace("'","") for x in df.columns] #改一下列名
    ##获得测序深度
    sums = getsumcount(list(df.columns))
    for column,value in df.iteritems():
        if column=="chr" or column=="start" or column=="end":
            continue
        df[column] = round(value / sums[column] * 1000000, 2)
    df.to_csv("/home/mwshi/project/CRC_enhancer/super_enhancer/super_enhancer_rpm.csv",index=0)

def calculate_tpm():
    df = pd.read_csv("/home/mwshi/project/CRC_enhancer/super_enhancer/CO_H3K27ac_readCounts.tab", sep="\t")
    df.columns = [x.strip("#").replace("'","") for x in df.columns] #改一下列名
    ##获得测序深度
    region_length = abs( df['end'] - df['start'])
    for column,value in df.iteritems():
        if column=="chr" or column=="start" or column=="end":
            continue
        sums = sum(df[column] / region_length)
        df[column] = round(value / region_length / sums * 1000000, 2)
    df.to_csv("/home/mwshi/project/CRC_enhancer/super_enhancer/super_enhancer_tpm.csv",index=0)  
    
    


if __name__ == "__main__":
    ##step 1
    merge_bed()
    lnfiles()
    ##step 2
    multiBamSummary()
    calculate_rpm()
    calculate_tpm()
    