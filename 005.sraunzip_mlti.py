import sys,os
import pandas as pd
import csv
import re
import json

sys.path.insert(1, "/home/zyang/software/MitEdit")
from basic import Basic

f = open("/home/mwshi/project/CEdb/srr2gsm.json","r")
tmp_dir = "/home/mwshi/project/CEdb/tmp_fastq/"
ok_srr = set([oo.split(".")[0].split("_")[0] for oo in os.listdir(tmp_dir)])
records = []
pbs_dir = "/home/mwshi/project/CEdb/pbs/"
for line in f:
    sample = json.loads(line)
    srrs=sample["srr"]
    all_srr = sample["srr"]+sample["input_srr"]+sample["nor_srr"]+sample["nor_input_srr"]
    if not all([word in ok_srr for word in all_srr]):
        print(sample["GSM"]+"the srr list is not complete")
        continue
    with open("/home/mwshi/project/CEdb/prejson/"+sample["GSM"]+".json","w") as w:
        w.write(line)
    cmd = "python /home/mwshi/github/CRC_enhancer/004.sraunzip.py /home/mwshi/project/CEdb/prejson/%s.json" % (sample["GSM"])
    pbs_file = os.path.join(pbs_dir, sample["GSM"] + ".mergesrr.pbs")
    handle = open(pbs_file, "w")
    handle.write(cmd)
    handle.close()
    #Basic.run(cmd, wkdir=pbs_dir)
    cmd = "qsub -l nodes=1:ppn=2 %s" % pbs_file
    print(cmd)
    #Basic.run(cmd, wkdir=pbs_dir)

        