import sys,os
import pandas as pd
import csv
import re
import json
import sys

sys.path.insert(1, "/home/zyang/software/MitEdit")
from basic import Basic

def transfer_file():
    #this function used to transfer the downloaded 
    sra_download_dir = "/home/mwshi/ncbi/public/sra"
    desti_dir = "/home/mwshi/project/CEdb/sra/"
    handle = open("/home/mwshi/project/CRC_enhancer/rawdata/Hotspots/Hotspots.txt", "r")
    for line in handle:
        sample = line.strip("\n")
        cmd = "mv /home/mwshi/ncbi/public/sra/%s.sra %s;" % (sample, desti_dir)       
        Basic.run(cmd)
        
def extract_fastq():
    sra_dir = "/home/mwshi/ncbi/public/sra/"
    pbs_dir = "/home/mwshi/project/CEdb/pbs"
    Basic.mkdir(pbs_dir)
    for one_file in [f for f in os.listdir(sra_dir) if f.endswith(".sra")]:
        in_file = os.path.join(sra_dir, one_file)
        if os.path.exists(os.path.join('/home/mwshi/project/CEdb/tmp_fastq/', one_file[:-4]+"_1.fastq.gz")):
            print(one_file + " is ok")
            continue
        if os.path.exists(os.path.join('/home/mwshi/project/CEdb/tmp_fastq/', one_file[:-4]+".fastq.gz")):
            print(one_file + " is ok")
            continue            
        pbs_file = os.path.join(pbs_dir, one_file + ".pbs")
        cmd = """fastq-dump --gzip --outdir /home/mwshi/project/CEdb/tmp_fastq/ --split-3 %s;""" %(in_file)
        handle = open(pbs_file, "w")
        handle.write(cmd)
        handle.close()
        #Basic.run(cmd, wkdir=pbs_dir)
        cmd = "qsub -q batch -l nodes=1:ppn=1 %s" % pbs_file
        Basic.run(cmd, wkdir=pbs_dir)

def catfastq(GSM,srrlist,layout,ttype):
    tmp_dir="/home/mwshi/project/CEdb/tmp_fastq/"
    if layout=="PAIRED":
        fq1 = []
        fq2 = []
        if len(srrlist)==1:
            srr = srrlist[0]
            if os.path.exists(os.path.join(tmp_dir, srr+"_1.fastq.gz")) and os.path.exists(os.path.join(tmp_dir, srr+"_2.fastq.gz")):
                cmd = "cp %s %s" % (os.path.join(tmp_dir, srr+"_1.fastq.gz"), "/home/mwshi/project/CEdb/rawdata/"+ttype+"_"+GSM+"_1.fastq.gz")
                Basic.run(cmd)
                cmd = "cp %s %s" % (os.path.join(tmp_dir, srr+"_2.fastq.gz"), "/home/mwshi/project/CEdb/rawdata/"+ttype+"_"+GSM+"_2.fastq.gz")
                Basic.run(cmd)
        elif len(srrlist)>1:
            for srr in srrlist:
                if os.path.exists(os.path.join(tmp_dir, srr+"_1.fastq.gz")) and os.path.exists(os.path.join(tmp_dir, srr+"_2.fastq.gz")):
                    fq1.append(os.path.join(tmp_dir, srr+"_1.fastq.gz"))
                    fq2.append(os.path.join(tmp_dir, srr+"_2.fastq.gz"))
                else:
                    print(srr+"is error, run is not exist")
                    continue
            cmd = "cat %s > %s" % (" ".join(fq1), "/home/mwshi/project/CEdb/rawdata/"+ttype+"_"+GSM+"_1.fastq.gz")
            Basic.run(cmd)
            cmd = "cat %s > %s" % (" ".join(fq2), "/home/mwshi/project/CEdb/rawdata/"+ttype+"_"+GSM+"_2.fastq.gz")
            Basic.run(cmd)
    if layout=="SINGLE":
        fq1 = []
        if len(srrlist)==1:
            srr=srrlist[0]
            if os.path.exists(os.path.join(tmp_dir, srr+".fastq.gz")):
                cmd = "cp %s %s" % (os.path.join(tmp_dir, srr+".fastq.gz"), "/home/mwshi/project/CEdb/rawdata/"+ttype+"_"+GSM+".fastq.gz")
                Basic.run(cmd)
        elif len(srrlist)>1:                
            for srr in srrlist:
                if os.path.exists(os.path.join(tmp_dir, srr+".fastq.gz")):
                    fq1.append(os.path.join(tmp_dir, srr+".fastq.gz"))
                else:
                    print(srr+"is error, run is not exist")
                    continue
            cmd = "cat %s > %s" % (" ".join(fq1), "/home/mwshi/project/CEdb/rawdata/"+ttype+"_"+GSM+".fastq.gz")            
            Basic.run(cmd)

def merge_fastq(jsfile):
    sample = json.load(open(jsfile,"r"))
    srrs=sample["srr"]
    if sample["layout"]=="PAIRED":
        if os.path.exists("/home/mwshi/project/CEdb/rawdata/cancer"+"_"+sample["GSM"]+"_1.fastq.gz") and os.path.exists("/home/mwshi/project/CEdb/rawdata/cancer"+"_"+sample["GSM"]+"_2.fastq.gz"):
            print(sample["GSM"]+" is already merge")
            continue
        record = {"GEO":sample["GEO"],"GSM":sample["GSM"],"TUSSUE":sample["TUSSUE"],"CANCER_TYPE":sample["CANCER_TYPE"],
        "fq":{"fq1":"/home/mwshi/project/CEdb/rawdata/cancer"+"_"+sample["GSM"]+"_1.fastq.gz",
              "fq2":"/home/mwshi/project/CEdb/rawdata/cancer"+"_"+sample["GSM"]+"_2.fastq.gz"},
        "INPUT":None,"normal":None,"norinp":None,"layout":"PAIRED"}
        catfastq(sample["GSM"],srrs,"PAIRED","cancer")
        if len(sample["input_srr"])>0:
            catfastq(sample["GSM"],sample["input_srr"],"PAIRED","input")
            record['INPUT']={"fq1":"/home/mwshi/project/CEdb/rawdata/input"+"_"+sample["GSM"]+"_1.fastq.gz",
                             "fq2":"/home/mwshi/project/CEdb/rawdata/input"+"_"+sample["GSM"]+"_2.fastq.gz"}
        if len(sample["normal"])!="":
            catfastq(sample["normal"]+"_"+sample["GSM"],sample["nor_srr"],"PAIRED","normal")
            record['normal']={"fq1":"/home/mwshi/project/CEdb/rawdata/normal"+"_"+sample["normal"]+"_"+sample["GSM"]+"_1.fastq.gz","fq2":"/home/mwshi/project/CEdb/rawdata/normal"+"_"+sample["normal"]+"_"+sample["GSM"]+"_2.fastq.gz"}
        if len(sample["normal_input"])!="":
            catfastq(sample["normal"]+"_"+sample["GSM"],sample["nor_input_srr"],"PAIRED","norinp")
            record['norinp']={"fq1":"/home/mwshi/project/CEdb/rawdata/norinp"+"_"+sample["normal"]+"_"+sample["GSM"]+"_1.fastq.gz",
                             "fq2":"/home/mwshi/project/CEdb/rawdata/norinp"+"_"+sample["normal"]+"_"+sample["GSM"]+"_2.fastq.gz"}
    elif sample["layout"]=="SINGLE":
        if os.path.exists("/home/mwshi/project/CEdb/rawdata/cancer"+"_"+sample["GSM"]+".fastq.gz"):
            print(sample["GSM"]+" is already merge")
            continue
        record = {"GEO":sample["GEO"],"GSM":sample["GSM"],"TUSSUE":sample["TUSSUE"],"CANCER_TYPE":sample["CANCER_TYPE"],
        "fq":{"fq1":"/home/mwshi/project/CEdb/rawdata/cancer"+"_"+sample["GSM"]+".fastq.gz"},
        "INPUT":None,"normal":None,"norinp":None,"layout":"SINGLE"}
        catfastq(sample["GSM"],srrs,"SINGLE","cancer")
        if len(sample["input_srr"])>0:
            catfastq(sample["GSM"],sample["input_srr"],"SINGLE","input")
            record['INPUT']={"fq1":"/home/mwshi/project/CEdb/rawdata/input"+"_"+sample["GSM"]+".fastq.gz"}
        if len(sample["normal"])!="":
            catfastq(sample["normal"]+"_"+sample["GSM"],sample["nor_srr"],"SINGLE","normal")
            record['normal']={"fq1":"/home/mwshi/project/CEdb/rawdata/normal"+"_"+sample["normal"]+"_"+sample["GSM"]+".fastq.gz"}
        if len(sample["normal_input"])!="":
            catfastq(sample["normal"]+"_"+sample["GSM"],sample["nor_input_srr"],"SINGLE","norinp")
            record['norinp']={"fq1":"/home/mwshi/project/CEdb/rawdata/norinp"+"_"+sample["normal"]+"_"+sample["GSM"]+".fastq.gz"}
    else:
        print(sample["GSM"]+" has no layout") 
        continue
    with open("/home/mwshi/project/CEdb/metajson/"+sample["GSM"]+".json", 'w') as w:
        json.dump(record,w)


        
if __name__ == "__main__":
    merge_fastq(sys.argv[1])






    
    