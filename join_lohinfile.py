#! /usr/bin/python

import sys

def usage():
    print("python join_lohinfile.py primary1_loh.infile primary2_loh.infile relapse1_loh.infile relapse2_loh.infile")
    sys.exit(1)

def read_into_hash(window_hash, filename, id):
    fh = open(filename, 'r')
    window_hash[id] = {}
    for l in fh:
        l = l.rstrip("\n")
        chr, pos, nfreq, tfreq, status = l.split()
        locus = chr + ":" + pos
        window_hash[id][locus] = tfreq

def write_hash(window_hash):
    p1_hash = {}
    p2_hash = {}
    r1_hash = {}
    r2_hash = {}
    p1_hash = window_hash["p1"]
    p2_hash = window_hash["p2"]
    r1_hash = window_hash["r1"]
    r2_hash = window_hash["r2"]
    print "chr" + "\t" + "pos" + "\t" + "Primary1" + "\t" + "Primary2" + \
        "\t" + "Relapse1" + "\t" + "Relapse2"
    #Go through all of p1
    for key in p1_hash:
        p1_val = p1_hash[key]
        if key not in p2_hash:
            p2_val = "NA"
        else:
            p2_val = p2_hash[key]
        if key not in r1_hash:
            r1_val = "NA"
        else:
            r1_val = r1_hash[key]
        if key not in r2_hash:
            r2_val = "NA"
        else:
            r2_val = r2_hash[key]
        print "\t".join(key.split(":")) + "\t" + p1_val + "\t" + p2_val + \
                "\t" + r1_val + "\t" + r2_val
    #Go through P2 keys not in P1
    for key in p2_hash:
        if key in p1_hash:
            continue
        p1_val = "NA"
        p2_val = p2_hash[key]
        if key not in r1_hash:
            r1_val = "NA"
        else:
            r1_val = r1_hash[key]
        if key not in r2_hash:
            r2_val = "NA"
        else:
            r2_val = r2_hash[key]
        print "\t".join(key.split(":")) + "\t" + p1_val + "\t" + p2_val + \
                "\t" + r1_val + "\t" + r2_val
    #Go through R1 keys not in P1 and P2
    for key in r1_hash:
        if key in p1_hash or key in p2_hash:
            continue
        p1_val = "NA"
        p2_val = "NA"
        r1_val = r1_hash[key]
        if key not in r2_hash:
            r2_val = "NA"
        else:
            r2_val = r2_hash[key]
        print "\t".join(key.split(":")) + "\t" + p1_val + "\t" + p2_val + \
                "\t" + r1_val + "\t" + r2_val
    #Go through R2 keys not in P1, P2 and R1
    for key in r2_hash:
        if key in p1_hash or key in p2_hash or key in r1_hash:
            continue
        p1_val = "NA"
        p2_val = "NA"
        r1_val = "NA"
        r2_val = r2_hash[key]
        print "\t".join(key.split(":")) + "\t" + p1_val + "\t" + p2_val + \
                "\t" + r1_val + "\t" + r2_val

def main():
    p1 = sys.argv[1]
    p2 = sys.argv[2]
    r1 = sys.argv[3]
    r2 = sys.argv[4]
    sys.stderr.write(p1 + "\n" +  p2 + "\n" + r1 + "\n" + r2)
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

