#!/usr/bin/env python3

# Args:
# 1: file with list of samples
# 2: suffix without a dash

import sys

with open(sys.argv[1], "r") as f:
    x = f.readlines()
    f.close()

assemblies = ""
readsL = ""
readsR = ""
for i in x:
    assemblies += f"{i.strip()}-{sys.argv[2]}.fasta,"
    readsL += f"{i.strip()}_R1_001.fastq,"
    readsR += f"{i.strip()}_R2_001.fastq,"
assemblies = assemblies[0:-1]
readsL = readsL[0:-1]
readsR = readsR[0:-1]

with open("assemblies.txt", "w") as f:
    f.write(assemblies)
    f.close()
with open("readsR.txt", "w") as f:
    f.write(readsR)
    f.close()
with open("readsL.txt", "w") as f:
    f.write(readsL)
    f.close()
