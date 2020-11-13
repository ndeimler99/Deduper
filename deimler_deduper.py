#!/usr/bin/env python

import argparse
import re
import os

def get_args():
    """argparse function that takes user input"""
    parser = argparse.ArgumentParser(description = "Deduper - To aid in the removal of PCR Duplicates from a SAM file based on chromosome number, UMI, and read start position; A conda environment with samtools version 1.7 must be activated")
    parser.add_argument("-f", "--file", help = "what is your input SAM file?", required = True)
    parser.add_argument("-u", "--umi", help="What are your UMI's? should be in a text file with an umi on each row of the file", required=True)
    parser.add_argument("-o", "--output", help = "what is the name you wish on your output file?", default = "output", required=False)
    parser.add_argument("-p", "--paired", help="True for paired end. If this is specified as true the program will end", required = False, default = False)
    #argparse automatically has a -h --help function which returns the help messages specified above
    #parser.add_argument("-h", "--help_function", help="This program allows for the import of an unsorted sam file and a file containing a list of all used UMIS. (-f and -u respectively).  It will return a deduplicated sam file, a sam file containing all duplicates, a sam file containing unknown umis, and a sam file containing unmapped reads. The prefix of these files can be set with the -o option", required=False)
    return parser.parse_args()

args = get_args()

input_sam = args.file
umi_file = args.umi
output_name = args.output
paired_end = args.paired

if paired_end == 'True':
    print("This program cannot deduplicate a paired end sam file at this point in time")
    exit()


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
        #s must be first letter
        first_cigar_letter = re.findall(r"\D", cigar)[0]
        if first_cigar_letter == "S":
            soft = cigar.split("S")[0]
            actual_pos = int(start_pos) - int(soft)
        else:
            actual_pos = start_pos
    else:
        #I only care about what happens after soft clipping
        first_cigar_letter = re.findall(r"\D", cigar)[0]
        if first_cigar_letter == "S":
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
    """this function takes in an incorrect umi and compares it to those in the dictionary trying to correct it, it will only work given there is a single nucleotide differene"""

    for item in umi_dict:
        mistake = 0
        for i in range(len(umi)):
            if item[i] != umi[i]:
                mistake += 1
        if mistake == 1:
            umi = item
    return umi

def print_headers(sam):
    """goes through original sam files header lines and prints them to output files"""
    with open(sam, "r") as sam:
        for line in sam:
            if line.startswith("@"):
                #print("@")
                dedup_out.write(line)
                duplicates.write(line)
                unmapped.write(line)
                unknown.write(line)
            else:
                #print("else")
                return True

def main_function(umi_file, sam_file, output_name):
    reads = create_UMI_dict(umi_file, [])
    count_dup = create_UMI_dict(umi_file, 0)
    current_chrom = -1
    current_strand = ""
    count=0

    with open(sam_file, "r") as sam:
        #open original input sam file
        for line in sam:
            #if its a header line 
            if line.startswith("@"):
                pass
            #if it isnt a header line
            else:
                line = line.strip()
                factors = get_Factors(line)
                #print(factors)
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
                        #if current_strand != strand:
                            #print("strandd_switched")
                        reads = {item:[] for item in reads}
                        #print("dict cleared")
                        current_chrom = chrom
                    if current_chrom == chrom:
                        #we are on same chromosome and same direction
                        if umi not in count_dup:
                            #if UMI does not appear in UMI count dictionary (error correct UMI)
                            umi = error_correct(factors[0], count_dup)
                            #print("error correct")
                        if umi in count_dup:
                            #if UMI appears in UMI count dictionary
                            
                            if start_pos in reads[umi]:
                                #check forward dictionary[umi] for starting position
                                    #increment count[umi]
                                count_dup[umi] += 1
                                    #write out read to duplicated
                                duplicates.write(line+ "\n")
                            else:
                                    #write to good_reads
                                dedup_out.write(line+ "\n")
                                    #add to forward dictionary
                                reads[umi].append(start_pos)
                                #count+=1
                                #print(count)

                                    #print(forward[umi])
                        else:
                            #write out to UMI_not_found file
                            unknown.write(line+ "\n")
    return "Finished"




#samtools commands to create files go here, should create a positive and negative.sam temp file

sorted_sam = "samtools sort -O sam {} > ./{}.sorted.sam".format(input_sam, output_name)
os.system(sorted_sam)
positive_command = "samtools view -F 16 ./{}.sorted.sam > ./positive.sam".format(output_name)
os.system(positive_command)
negative_command = "samtools view -f 16 ./{}.sorted.sam > ./negative.sam".format(output_name)
os.system(negative_command)

positive = "./positive.sam"
negative = "./negative.sam"

with open("{}.deduplicated.out".format(output_name), "w") as dedup_out, open("{}.duplicates.out".format(output_name), "w") as duplicates, open("{}.unmapped.out".format(output_name),"w") as unmapped, open("{}.unknownUMI.out".format(output_name), "w") as unknown:
    #open all files
    #print header lines to all files
    print_headers(input_sam)

    #parses through original samfile (only through lines that start with @)
    main_function(umi_file, positive, output_name)
    main_function(umi_file, negative, output_name)

#remove temporary files
os.system("rm ./positive.sam")
os.system("rm ./negative.sam")
