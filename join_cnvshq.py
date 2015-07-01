#! /usr/bin/python

import sys

def usage():
    print("python join_cnvshq.py primary1_cnvs.hq primary2_cnvs.hq relapse1_cnvs.hq relapse2_cnvs.hq")
    sys.exit(1)

def read_into_hash(window_hash, filename, id):
    print id;
    fh = open(filename, 'r')
    window_hash[id] = {}
    for l in fh:
        l = l.rstrip("\n")
        chr, pos, tumor, normal, diff = l.split()
        locus = chr + ":" + pos
        window_hash[id][locus] = tumor

def write_hash(window_hash):
    p1_hash = {}
    p1_hash = window_hash["p1"]
    for key in p1_hash:
        if key in window_hash["r1"] and key in window_hash["r2"] and key in window_hash["p2"]:
            print "\t".join(key.split(":")) + "\t" + window_hash["p1"][key] + "\t" + window_hash["p2"][key] + \
                    "\t" + window_hash["r1"][key] + "\t" + window_hash["r2"][key]

def main():
    p1 = sys.argv[1]
    p2 = sys.argv[2]
    r1 = sys.argv[3]
    r2 = sys.argv[4]
    print(p1, p2, r1, r2)
    window_hash = {}
    read_into_hash(window_hash, p1, "p1")
    read_into_hash(window_hash, p2, "p2")
    read_into_hash(window_hash, r1, "r1")
    read_into_hash(window_hash, r2, "r2")
    write_hash(window_hash)

if __name__ == "__main__":
    if(len(sys.argv) != 5):
        usage()
    main()

