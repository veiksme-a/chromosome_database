#!/usr/bin/env python3

"""
Program:    codon_usage_calc
File:       codon_usage_calc.py

Version:    V1.0
Date:       27.04.18
Function:   Calculate codon usage for the entire chromosome (coding regions
            only), including frequency per 100 codons and amino acid ratio,
            and write INSERT statements for each codon.

Author:     AO
            
--------------------------------------------------------------------------

This program was produced as a part of Bioinformatics MSc programme assessment.

--------------------------------------------------------------------------
Description:
============
A set of functions for calculating codon frequency (per 100 codons), linking
codons to corresponding amino acids and finding amino acid ratios. The result
is a set of INSERT statements for each of 64 codons suitable for uploading
into the chromosome database (codon_usage table).

--------------------------------------------------------------------------
Usage:
======
The codon_parser() function combines all other functions to produce the output,
hence the user will obtain the result by executing this function. User needs to
provide a file that contains codons with their amino acid translations and a
file where the INSERT statements will be written.

--------------------------------------------------------------------------
Revision History:
=================
V1.0   27.04.18   Original   By: AO
"""

#*************************************************************************
# Import libraries

import pymysql.cursors
import re
from collections import Counter
from collections import defaultdict



#*************************************************************************

def codon_counts():

    """Extract all CDS sequences from the database, split in codons and count
    each codon occurrence. The function only includes 64 existing codons and
    any codons that might contain extra characters like "n".

    Return: output              --- A dictionary, where the keys are the
                                codons and the counts as values.

    27.04.18  Original   By: AO
    """

    # extract all CDS sequences
    query = 'select whole_seq from chr_db.cds'
    db = pymysql.connect(db='chr_db', user='root', passwd='', host='localhost', port=3306)

    cursor = db.cursor()

    q = cursor.execute(query)

    p_seq_string = re.compile(r"\('(.*)',\)")
    # thimine is not present in RNA, so it will have to be removed
    p_thimine = re.compile(r't')

    all_cds_joined = ''


    for sequence in cursor:
        match = p_seq_string.search(str(sequence))
        seq_string = match.group(1)
        # check if the string length is correct length, since codons composed
        # of 3 nucleotides
        if len(seq_string)%3==0:
            rna_seq = p_thimine.sub('u', seq_string)
            # join all CDS sequences into one
            all_cds_joined += rna_seq
        else:
            pass


    n = 3
    # split the joined CDS sequence into 3-nt long codons
    codons = [all_cds_joined[i:i+n] for i in range(0, len(all_cds_joined), n)]

    codon_content = re.compile(r'[augc]{3}')
    #to avoid any extra characters like 'n' for unidentified nucleotides

    codon_list=[] # correct length of this list would be 64
    for i in codons:
        match = codon_content.search(i)
        if match:
            # if a codon contains only A, T, C, G, add them to codon list
            codon_list.append(i)
        else:
            pass
    # count codon occurrence
    counted_codons = dict(Counter(codon_list))

    return(counted_codons)



#*************************************************************************


def codon_names(file):

    """Create a codon map, where an amino acid is matched to all codons
    encoding it. The information is obtained from a file, which avoids using
    Biopython or other specific libraries, since they are not allowed for
    the assignment.

    Input:  file                --- An input file with one codon and a matching
                                amino acid per line.
    Return: output              --- A dictionary, where the keys are 
                                amino acids and the values contain
                                all codons encoding particular amino acid.

    27.04.18  Original   By: AO
    """

    # create dictionary where a value can be a list with miltiple elements
    codon_dict = defaultdict(list)

    p = re.compile(r'([AUCG]{3}).*\s(\w{3})\s')
    with open(file, 'r') as f:
        lines = f.read().splitlines()

        for line in lines:
            match = p.search(line)
            codon = str(match.group(1))
            amino_acid = str(match.group(2))
            # forms a dictionary with amino acids as keys and all corresponding
            # codons as values
            codon_dict[amino_acid.lower()].append(codon.lower())

    return(codon_dict)



#*************************************************************************


def frequencies(codon_count_dict):

    """Count codon frequencies from counts.

    Input:  codon_count_dict    --- A dictionary with codons as keys and
                                their counts as values.
    Return: codon_percentages   --- A dictionary, where the keys are 
                                codons and the values show the occurrence
                                percentages for eacg codon (per 100).

    27.04.18  Original   By: AO
    """
    
    triplet_sum = 0
    for k,v in codon_count_dict.items():
        triplet_sum += v # calculates a total sum of all codon counts

    codon_percentages={}
    for k,v in codon_count_dict.items():
        percentage = v/triplet_sum * 100 #calculates an occurrence percentage for each codon
        codon_percentages[k] = round(percentage,2)
        # percentages were rounded to two significant figures, as it wasd in the
        # example provided

    return(codon_percentages)



#*************************************************************************


def aa_ratios(codon_count_dict, codon_name_dict):

    """Count amino acid ratios: the proportion of occurrence of each codon
    in context of an amino acid it encodes.

    Input:  codon_count_dict    --- A dictionary with codons as keys and
                                their counts as values.
            codon_name_dict     --- A dictionary with amino acids as keys and
                                their codons as values.
    Return: aa_ratio_dict       --- A dictionary, where the keys are 
                                codons and the values are ratios of 
                                occurrence regarding a particular amino acid.

    27.04.18  Original   By: AO
    """

    aa_ratio_dict = {}

    for k,v in codon_name_dict.items():
        aa_count = 0
        # count total number of codons encoding a particular amino acid
        for i in v:
            codon_count = int(codon_count_dict[i])
            aa_count += codon_count
        # calculate the ratio for each codon e+ncoding the amino acid 
        for i in v:
            if i in codon_count_dict.keys():
                aa_ratio = round(codon_count_dict[i]/aa_count, 2)
                aa_ratio_dict[i] = aa_ratio
        
    return(aa_ratio_dict)



#*************************************************************************


class Codon:
    """ Class for handling the codon data. Attributes contain a codon sequence
    its frequence, amino acid it is encoding and corresponding amino acid ratio

    27.04.18   Original   By: AO
    """


    def __init__(self, codon, frequency, codon_to_aa_dict, aa_ratio_dict):
        self.codon = codon
        self.frequency = frequency

        for k,v in codon_to_aa_dict.items():
            # asign corresponding amino acid
            if codon in v:
                self.amino_acid = k
            else:
                pass


        for k,v in aa_ratio_dict.items():
            if codon == k:
                # get a relevant amino acid ratio from the precalculated values
                self.aa_ratio = v
            


#*************************************************************************


def codon_parser(codon_file, output_file):

    """Parser creates INSERT statements for the chromosomal codon usage table
    in the assignment database. Output contains 64 statements - for each codon.

    Input:  codon_file          --- File for codon-amino acid mapping.
            output_file         --- A file where the user wants INSERT statements
                                to be written to.
    Return: Writes 64 INSERT statements into a file of choice.

    27.04.18  Original   By: AO
    """
    
    count_data = codon_counts() 
    frequency_data = frequencies(count_data)
    codon_map = codon_names(codon_file)
    amino_acid_ratios = aa_ratios(count_data, codon_map)
    # generate INSERT statements for each codon
    insert_strings = ''
    for codon, frequency in frequency_data.items():
        entry = Codon(codon, frequency, codon_map, amino_acid_ratios)

        codon_insert_string = "INSERT INTO codon_usage(codon, frequency_100, ratio_aa, amino_acid, chr_name) \nVALUES ('" + str(entry.codon) + "', " + str(entry.frequency) + ", " + str(entry.aa_ratio) + ", '" + str(entry.amino_acid) + "', 'chromosome 16');\n"
        insert_strings += codon_insert_string

    with open(output_file, "a") as db_entries:
        db_entries.write(insert_strings)

        
        

    
