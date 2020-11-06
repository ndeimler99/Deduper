#!/usr/bin/env python

import argparse
import re

def get_args():
    parser = argparse.ArgumentParser(description = "Deduper - To aid in the removal of PCR Duplicates")
    parser.add_argument("-f", "--sam_file", help = "what is your input SAM file?", required = True)
    parser.add_argument("-u", "--umi_file", help="What are your UMI's?", required=True)
    parser.add_argument("-o", "--output", help = "what is the name you wish on your output file?", default = "output", required=False)
    parser.add_argument("-p", "--paired", help="True for paired end", required = False, default = False)
    return parser.parse_args()

args = get_args()

input_sam = args.sam_file
umi_file = args.umi_file
output_name = args.output
paired_end = args.paired

#sort file by strandedness, chrom, UMI, left most, quality?

def create_UMI_dict(umi_file, number): 
    """This function takes the input of a file containing a list of UMI's that should appear in the SAM file and adds them to a dictionary as the keys and assigns each a value of 0"""
    #initialize empty dictionary
    umi_dict = {}
    #open the file
    with open(umi_file, "r") as file:
        #parse through the file line by line
        for line in file:
            #remove the newline character
            line = line.strip()
            #add to dictionary as key with value of 0
            umi_dict[line] = number
    #return finalized dictionary
    return(umi_dict)

def get_Factors(read):
    """This Function takes a strippined line of the sam file and dissects it to retrieve UMI, chromosome #, starting position, and strandedness"""
    attribute = read.split("\t")
    chrom = attribute[2]
    UMI = attribute[0].split(":")[-1]
    if(int(attribute[1]) & 16 == 16):
        stranded = "reverse"
    else:
        stranded="forward"
    if(int(attribute[1])&4 == 4):
        mapped=False
    else:
        mapped=True
    
    start_pos = calc_start(stranded, attribute[3], attribute[5])
    #print(start_pos)
    return UMI, chrom, stranded, start_pos, mapped


def calc_start(stranded, start_pos, cigar):
    """function that aids in the calculation of starting position based on strandedness"""
    #to calculate starting position if forward read add soft clipping if reverse read do not add soft clipping
    if stranded == "forward":
        if "S" in cigar:
            soft = cigar.split("S")[0]
            actual_pos = int(start_pos) - int(soft)
        else:
            actual_pos = start_pos
    else:
        #I only care about what happens after soft clipping
        if "S" in cigar:
            cigar = cigar.split("S")[1]
            num = re.findall(r"\d+", cigar)
            num = [int(i) for i in num]
            actual_pos = int(start_pos) + sum(num)
        else:
            num = re.findall(r"\d+", cigar)
            num = [int(i) for i in num]
            actual_pos = int(start_pos) + sum(num)
    
    return actual_pos

def error_correct(umi, umi_dict):
    for item in umi_dict:
        mistake = 0
        for i in range(len(umi)):
            if item[i] != umi[i]:
                mistake += 1
        if mistake == 1:
            umi = item
    return umi


def main_function(umi_file, sam_file, output_name):
    forward = create_UMI_dict(umi_file, [])
    reverse = create_UMI_dict(umi_file, [])
    count_dup = create_UMI_dict(umi_file, 0)
    current_chrom = -1
    #print(forward)
    #create reverse umi dictionary
    #create forward UMI dictionary
    #create UMI count dictionary
    #do it seperately since dictionaries point to eachother
    with open(sam_file, "r") as sam, open("{}.deduplicated.out".format(output_name), "w") as dedup_out, open("{}.duplicates.out".format(output_name), "w") as duplicates, open("{}.unmapped.out".format(output_name),"w") as unmapped, open("{}.unknownUMI.out".format(output_name), "w") as unknown:
        for line in sam:
            if line.startswith("@"):
                dedup_out.write(line)
                duplicates.write(line)
                unmapped.write(line)
                unknown.write(line)
            else:
                line.strip()
                factors = get_Factors(line)

                umi = factors[0]
                chrom = factors[1]
                strand = factors[2]
                start_pos = int(factors[3])
                #print(factors[3])
                #print(start_pos)
                if factors[4] == False:
                #check for mapping first
                    unmapped.write(line + "\n")
                else:
                    #it mapped
                    if current_chrom != factors[1]:
                        #meaning current chrom progressed
                        forward = {item:[] for item in forward}
                        reverse = {item:[] for item in reverse}
                        current_chrom = factors[1]
                    if current_chrom == factors[1]:
                        #we are on same chromosome
                        if umi not in count_dup:
                            #if UMI does not appear in UMI count dictionary (error correct UMI)
                            umi = error_correct(factors[0], count_dup)
                        if umi in count_dup:
                            #if UMI appears in UMI count dictionary
                            if strand == 'forward':
                                if start_pos in forward[umi]:
                                #check forward dictionary[umi] for starting position
                                    #increment count[umi]
                                    count_dup[umi] += 1
                                    #write out read to duplicated
                                    duplicates.write(line + "\n")
                                else:
                                    #write to good_reads
                                    dedup_out.write(line + "\n")
                                    #add to forward dictionary
                                    forward[umi].append(start_pos)
                                    #print(forward[umi])
                            else:
                                #read is reverse
                                if start_pos in reverse[umi]:
                                    #increment count[umi]
                                    count_dup[umi] += 1
                                    duplicates.write(line + "\n")
                                    #write out to read duplicated
                                    pass
                                else:
                                    #write to good_reads
                                    dedup_out.write(line + "\n")
                                    #add to reverse dictionary
                                    reverse[umi].append(start_pos)
                                    pass
                        else:
                            #write out to UMI_not_found file
                            unknown.write(line+"\n")
    return "Finished"

main_function(umi_file, input_sam, output_name)

# def sort_SAM(file_name):
# 	```takes SAM input file specified by argparse and sorts it by chromsome and then starting position```
# 	Call bash command that sorts SAM file?
# 	If this doesnt work just create a bash script that sorts SAM file and then calls this python script

