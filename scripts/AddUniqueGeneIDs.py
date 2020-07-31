#!/usr/bin/env python
import argparse
import os 
import sys
import re

# global var for inputs
args=None 

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("infile", help="positional input")
    parser.add_argument("-s", "--string", help="string option")
    parser.add_argument("-n", "--number", help="numeric option", type=int, default=5)
    parser.add_argument("-l", "--list", nargs="*", help="list with zero or more entries")
    parser.add_argument("-l2", "--list2", nargs="+", help="list one or more entries")
    parser.add_argument('-d', help="store args.d as true if -d",  action="store_true", default=False)
    args = parser.parse_args()

    seen = set()

    transcript = "fake"
    for line in open(args.infile):
        t = line.strip().split()
        line = line.replace("CDStopAdjusted","cDStopAdjusted")
        if(line[0]=="#"): 
            sys.stdout.write(line)
            continue
        atts = t[8].split(";")
        match = re.match("(Parent|ID)=([^;]+)", atts[0] )
        tag, ID = match.groups()
        typ = t[2]
        
        if(typ == "locus"):
            seen.add(transcript)
            sys.stdout.write(line)
            continue
        elif(typ == "transcript"):
            seen.add(transcript)
            transcript = ID
            if(transcript in seen):
                transcript = ID+"_dup"
        
        assert ID in transcript , ID + transcript

        sys.stdout.write(line.replace(ID, transcript))


