#!/usr/bin/env python3

"""
Retrieve annotations closest to the query position from a GTF file.
"""

import sys

def indexGTF(file):
    index = dict()
    with open(file, "r") as f:
        buffer = f.readlines()
    f.close()

    n = 0
    for line in buffer:
        if len(line.split("\t")) <= 1:
            n += 1
        else:
            break

    # Remove header
    # print(f"Removing {n} lines from GTF file")
    for i in range(n + 1):
        buffer.pop(0)

    # Remove bad footers from annotation file
    if len(buffer[len(buffer) - 1].split("\t")) <= 1:
        buffer.pop(len(buffer) - 1)

    for line in buffer:
        x = line.split("\t")
        try:
            name = x[0]
            source = x[1]
            feature = x[2]
            start = x[3]
            stop = x[4]
        except IndexError:
            print(f"Bad line: {line}")

        index[name] = [feature, start, stop]

    return index

def showNearestStart(query, index):
    print(f"> Showing nearest items for starting position query: {query}")

    x = list()
    for i in index:
        pos = index[i][1]
        dist = int(query) - int(pos)
        index[i] = index[i].append(dist)
        x.append([i, dist])

    y = list()
    for i in x:
        y.append(abs(i[1]))
    print(min(y))

query = sys.argv[1]
index = indexGTF(sys.argv[2])

showNearestStart(query, index)

# 1. Index GTF file (get start and end positions, with annotation name)
# 2. Get query position, and get distance from start position
# 3. Print annotations with lowest distance
