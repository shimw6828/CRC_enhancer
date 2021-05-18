import sys,os
import pandas as pd
import csv
import re

sys.path.insert(1, "/home/zyang/software/MitEdit")
from basic import Basic

def transfer_file():
    #this function used to transfer the downloaded 
    sra_download_dir = "/home/mwshi/ncbi/public/sra"
    desti_dir = "/home/mwshi/project/CRC_enhancer/rawdata/Hotspots/sra/"
    
    handle = open("/home/mwshi/project/CRC_enhancer/rawdata/Hotspots/Hotspots.txt", "r")
    for line in handle:
        sample = line.strip("\n")
        cmd = "mv /home/mwshi/ncbi/public/sra/%s.sra %s;" % (sample, desti_dir)       
        Basic.run(cmd)

def extract_fastq():
    sra_dir = "/home/mwshi/project/CRC_enhancer/rawdata/Hotspots/sra"
    pbs_dir = "/home/mwshi/project/CRC_enhancer/pbs"
    for one_file in os.listdir(sra_dir):
        in_file = os.path.join(sra_dir, one_file)
        pbs_file = os.path.join(pbs_dir, one_file + ".pbs")
        cmd = """fastq-dump --gzip --outdir /home/mwshi/project/CRC_enhancer/rawdata/Hotspots/fastq --split-3 %s;""" %(in_file)
        handle = open(pbs_file, "w")
        handle.write(cmd)
        handle.close()
        #Basic.run(cmd, wkdir=pbs_dir)
        cmd = "qsub -q batch -l nodes=1:ppn=1 %s" % pbs_file
        Basic.run(cmd, wkdir=pbs_dir)

def parse_metadata():
    fastq_dir = "/home/mwshi/project/CRC_enhancer/rawdata/Hotspots/fastq"
    meta_data = "/home/mwshi/project/CRC_enhancer/rawdata/Hotspots/meta_hotspot.txt"
    new_fastq_dir = "/home/mwshi/project/CRC_enhancer/rawdata/Hotspots/rename_fastq/"
    Basic.mkdir(new_fastq_dir)
    meta_df = pd.read_csv(meta_data, sep=",")
    #print(meta_df.columns)
    meta_sub = meta_df[["Run", "source_name", "Tissue", "chip_antibody", "Sample Name"]]
    print(meta_sub[0:5])  
    #single
    for one_sample in meta_sub.iterrows():
        print(one_sample)
        sample_id = one_sample[1]['Run']
        if one_sample[1]['Tissue']=='colorectal cancer cell line':
            sample_type = "Cellline"
        elif one_sample[1]['Tissue']=='primary colorectal cancer tumor':
            sample_type = "tumor"
        elif one_sample[1]['Tissue']=='normal colon crypt':
            sample_type = "normal"
        if 'H3K27ac' in one_sample[1]['Tissue']:
            marker = "H3K27ac"
        elif 'H3K27me3' in one_sample[1]['Tissue']:
            marker = "H3K27me3"
        elif 'H3K4me1' in one_sample[1]['Tissue']:
            marker = "H3K4me1"
        else:
             marker = "ChIP_Input"
        new_id = "_".join([sample_type,marker,one_sample[1]['Sample Name'],"Hotspot"])
        old_path = os.path.join(fastq_dir, one_sample[1]['Run'] + ".fastq.gz")
        new_path = os.path.join(new_fastq_dir, new_id + ".fastq.gz")
        #print(new_path)
        if os.path.exists(new_path):
            print("error, soft link exits")
            exit(1)
        cmd = "ln -s %s %s" % (old_path, new_path)
        Basic.run(cmd)

 
 
def fastqc():
    merge_fastq_dir = "/home/mwshi/project/CRC_enhancer/rawdata/Hotspots/rename_fastq"
    fastqc_dir ="/home/mwshi/project/CRC_enhancer/fastqc/"
    pbs_dir = "/home/mwshi/project/CRC_enhancer/pbs/fastqc_job"
    Basic.mkdir(pbs_dir)
    Basic.mkdir(fastqc_dir)
    for one_sample in os.listdir(merge_fastq_dir):
        if ".log" in one_sample:
            continue
        if ".html" in one_sample:
            continue
        if ".json" in one_sample:
            continue
        
        file_path = os.path.join(merge_fastq_dir, one_sample)
        log_file = os.path.join(fastqc_dir, one_sample + ".log")
        cmd = "fastqc -t 7 -o %s -f fastq --noextract %s 2> %s" % (fastqc_dir, file_path, log_file)
        #Basic.run(cmd, wkdir= fastqc_dir)
         
        
        pbs_file = os.path.join(pbs_dir, one_sample + ".fastqc.pbs")
        handle = open(pbs_file, "w")
        handle.write(cmd)
        handle.close()
        #Basic.run(cmd, wkdir=pbs_dir)
        cmd = "qsub -l nodes=1:ppn=7 %s" % pbs_file
        print(cmd)
        Basic.run(cmd, wkdir=pbs_dir)
                   
    
def bowtie_index():
    #build bowtie 1 index
    #this function has not used in bowtie2, because I have build reference index long long ago
    #已经跑过，不用跑了
    reference_dir = "/home/zhluo/Project/CRC/data_nazhang/step56_human_CRC/reference_genome/"
    reference = "/home/zhluo/Project/CRC/data_nazhang/step56_human_CRC/reference_genome/GRCh37.p13.genome.fa"
    cmd = "bowtie-build --threads 20 %s %s" %(reference, reference)
    pbs_handle = open("/home/zhluo/Project/CRC/data_nazhang/step56_human_CRC/pbs/bowtie_idx.pbs", "w")
    pbs_handle.write(cmd)
    pbs_handle.close()
    cmd = "qsub -l nodes=1:ppn=20 %s" % ("/home/zhluo/Project/CRC/data_nazhang/step56_human_CRC/pbs/bowtie_idx.pbs")
    Basic.run(cmd, wkdir= reference_dir)

def cut_adapt():
    ##单端测序不跑这个
    fastq_dir = "/home/mwshi/project/CRC_enhancer/rawdata/Hotspots/rename_fastq"
    cut_adapt_fastq = "/home/mwshi/project/CRC_enhancer/rawdata/Hotspots/cut_fastq"
    #bowtie_dir = "/home/zhluo/Project/TF_enrichment/bowtie_mapping_paired"
    pbs_dir = "/home/mwshi/project/CRC_enhancer/pbs/cutadapt_pbs"
    Basic.mkdir(cut_adapt_fastq)
    Basic.mkdir(pbs_dir)
    #sample_list = []
    for one_sample in os.listdir(fastq_dir):
        if ".log" in one_sample:
            continue
        sample_id = one_sample.replace(".fastq.gz", "")
        #print(sample_id)
        fastq_1 = os.path.join(fastq_dir, one_sample)
        log_file = os.path.join(cut_adapt_fastq, sample_id + ".cut.log")
        cut_fastq_1 = os.path.join(cut_adapt_fastq, one_sample)
        cmd = "/home/zhluo/.local/bin/cutadapt -j 5 -u 10 -u -15 -U 10 -U -15 -m 40 -o %s %s 2> %s;" %(cut_fastq_1, fastq_1, log_file)
        #out_file = os.path.join(bowtie_dir, sample_id + ".sorted.bam")
        #tmp_prefix = os.path.join(bowtie_dir, sample_id + ".sorted.bam.tmp")
        #markdup_log = os.path.join(bowtie_dir, sample_id + ".mkdup.log")
        
        pbs_file = os.path.join(pbs_dir, "step3_" + one_sample + ".cutadapt.pbs")
        pbs_handle = open(pbs_file, "w")
        pbs_handle.write(cmd)
        pbs_handle.close()
        cmd = "qsub -l nodes=1:ppn=5 %s" % (pbs_file)
        #print(cmd)
        Basic.run(cmd, wkdir= pbs_dir)
    
def fastp():
    ###虽然短，但是还是跑一下
    fastq_dir = "/home/mwshi/project/CRC_enhancer/rawdata/Hotspots/rename_fastq/"
    fastp_dir = "/home/mwshi/project/CRC_enhancer/rawdata/Hotspots/fastp_fastq"
    outputDir = fastp_dir
    Basic.mkdir(fastp_dir)
    pbs_dir = "/home/mwshi/project/CRC_enhancer/pbs/cutadapt_pbs"
    for one_sample in os.listdir(fastq_dir):
        if ".log" in one_sample:
            continue
        if "_2.fastq.gz" in one_sample:
            continue
        sample_id = one_sample.replace(".fastq.gz", "")
        #print(sample_id)
        fastq1 = os.path.join(fastq_dir, one_sample)
        output1 = sample_id + ".fastq.gz"
        output1 = os.path.join(outputDir, output1)
        htmlPath = os.path.join(outputDir, sample_id + ".html")
        jsonPath = os.path.join(outputDir, sample_id + ".json")
        cmd = "fastp -z 7 -i %s -o %s -h %s -j %s" % (fastq1, output1, htmlPath, jsonPath)
        pbs_file = os.path.join(pbs_dir, "step2_" + sample_id + ".fastp.pbs")
        pbs_handle = open(pbs_file, "w")
        pbs_handle.write(cmd)
        pbs_handle.close()
        cmd = "qsub -l nodes=1:ppn=7 %s" % (pbs_file)
        #print(cmd)
        Basic.run(cmd, wkdir= pbs_dir)    


def bowtie_align_single():
    merge_fastq_dir = "/home/mwshi/project/CRC_enhancer/rawdata/Hotspots/fastp_fastq"
    bowtie_dir = "/home/mwshi/project/CRC_enhancer/chip_seq/Hotspots/"
    pbs_dir = "/home/mwshi/project/CRC_enhancer/pbs/mapping_single"
    Basic.mkdir(bowtie_dir)
    Basic.mkdir(pbs_dir)
    for one_sample in os.listdir(merge_fastq_dir):
        if ".log" in one_sample:
            continue
        if ".html" in one_sample:
            continue
        if ".json" in one_sample:
            continue
        sample_id = one_sample.replace(".fastq.gz", "")
        file_path = os.path.join(merge_fastq_dir, one_sample)
        log_file = os.path.join(bowtie_dir, sample_id + ".log")
        out_file = os.path.join(bowtie_dir, sample_id + ".sorted.bam")
        tmp_prefix = os.path.join(bowtie_dir, sample_id + ".sorted.bam.tmp")
        markdup_log = os.path.join(bowtie_dir, sample_id + ".mkdup.log")
        cmd = "bowtie --chunkmbs 320 -m 1 --best -p 7  /home/zhluo/Project/CRC/data_nazhang/step56_human_CRC/reference_genome/GRCh37.p13.genome.fa  -q %s -S 2> %s  | /home/zhluo/Project/CRC/data_nazhang/step56_human_CRC/samblaster/samblaster --removeDups | samtools view -Sb -F 4 - | samtools sort -m 2G -@ 7 -T %s -o %s 2> %s" % (file_path, log_file, tmp_prefix, out_file, markdup_log)
        pbs_file = os.path.join(pbs_dir, "step6_" + sample_id + ".mapping.pbs")
        pbs_handle = open(pbs_file, "w")
        pbs_handle.write(cmd)
        pbs_handle.close()
        cmd = "qsub -l nodes=1:ppn=7 %s" % (pbs_file)
        #print(cmd)
        Basic.run(cmd, wkdir= pbs_dir)
        
def mk_idx_flagstat():
    ##Input_GSM2058062文件质量太差，统计之后很小
    #bowtie_dir = "/home/zhluo/Project/TF_enrichment/bowtie_mapping_paired"
    bowtie_dir = "/home/mwshi/project/CRC_enhancer/chip_seq/Hotspots"
    pbs_dir = "/home/mwshi/project/CRC_enhancer/pbs/flagstat"
    Basic.mkdir(pbs_dir)
    for one_sample in os.listdir(bowtie_dir):
        print(one_sample)
        if ".flagstat" in one_sample:
            continue
        if ".log" in one_sample:
            continue
        if ".bai" in one_sample:
            continue
        if ".sorted.bam" in one_sample:
            file_path = os.path.join(bowtie_dir, one_sample)
            pbs_file = os.path.join(pbs_dir, "step6_" + one_sample + ".single.flagstat.pbs")
            #sample_id = one_sample.split(".")[0]
            flagstat_file = os.path.join(bowtie_dir, one_sample + ".flagstat")
            pbs_handle = open(pbs_file, "w")
            cmd = "samtools index %s;" % file_path
            pbs_handle.write(cmd)
            cmd = "samtools flagstat %s > %s;" % (file_path, flagstat_file)
            pbs_handle.write(cmd)
            pbs_handle.close()
            cmd = "qsub -l nodes=1:ppn=1 %s" % (pbs_file)
            Basic.run(cmd, wkdir= pbs_dir)

def makebw_single():
    downsample_dir = "/home/mwshi/project/CRC_enhancer/chip_seq/Hotspots"
    pbs_dir = "/home/mwshi/project/CRC_enhancer/pbs/makebw"
    bigwig_dir = "/home/mwshi/project/CRC_enhancer/makebw/Hotspots"
    Basic.mkdir(pbs_dir)
    Basic.mkdir(bigwig_dir)    
    
    #only sample, no substract input
    for one_file in os.listdir(downsample_dir):
        if ".log" in one_file:
            continue
        if ".bai" in one_file:
            continue
        file_path = os.path.join(downsample_dir, one_file)
        #sample_name = file_path.split(".")[0]
        #sample_name = sample_name.split("_")[2]
        #print(sample_name)
        #cmd = "bamCompare --bamfile1 %s --bamfile2 {input[0]} --normalizeUsingRPKM --ratio subtract --binSize 30 --smoothLength 300 -p 5  --extendReads 200 -o {output} 2> {log}" %(file_path)
        output_file = os.path.join(bigwig_dir, one_file + ".bw")
        log_file = os.path.join(bigwig_dir, one_file + ".bamCoverage.log")
        cmd = "bamCoverage -b %s --normalizeUsing RPKM --binSize 30 --smoothLength 300 -p 5 --extendReads 200 -o %s 2> %s;" %(file_path, output_file, log_file)
        pbs_file = os.path.join(pbs_dir, one_file + ".makebw.pbs")
        handle = open(pbs_file, "w")
        handle.write(cmd)
        handle.close()
        cmd = "qsub -l nodes=1:ppn=5 %s" % pbs_file
        Basic.run(cmd, wkdir=pbs_dir)
        
def rename_bam():
    bowtie_dir = "/home/zhluo/Project/CRC/data_nazhang/step56_human_CRC/bowtie_mapping"
    rename_dir = "/home/zhluo/Project/CRC/data_nazhang/step56_human_CRC/rename_bam"

if __name__ == "__main__":        
    ##step 1
    #transfer_file()
    ##step 2
    #extract_fastq()
    ##step3
    #parse_metadata()
    ##fastqc
    #fastqc()
    #cut_adapt()
    #fastp()
    #bowtie_index()
    #bowtie_align_paired()
    #bowtie_align_single()
    #mk_idx_flagstat()
    
    #mk_idx_flagstat()
    #create_chromhmm_matrix()
    #parse_flagstat()
    #makebw_paired()
    #makebw_single()
    rename_bam()
    #find_sample_input()
    #macs2_call_peaks()
    #rename_peak_file()
    #super_enhancer()