#!/usr/bin/env python3

"""
Program:    parser
File:       parser.py

Version:    V1.0
Date:       25.04.18
Function:   Extract information from the GenBank chromosome file and produce
            INSERT statements that can be loaded into the chromosome database -
            tables Chromosome, Locus, Accession and CDS.

Author:     AO
            
--------------------------------------------------------------------------

This program was produced as a part of Bioinformatics programme assessment.

--------------------------------------------------------------------------
Description:
============
A set of functions designed to retrieve required information from the
original GenBank file. GenBank entries contain information for each locus
separately. The parser extracts all information accosiated with a locus,
including its accession numbers and all CDS, and writes it into INSERT
statements suitable for uploading into a database.

--------------------------------------------------------------------------
Usage:
======
Since the parser() function combines all other functions together, the user
need to execute only this function, providing a GenBank file, species name,
chromosome and the name of the file where the INSERT statements will be
written into.

--------------------------------------------------------------------------
Revision History:
=================
V1.0   25.04.18   Original   By: AO
"""

#*************************************************************************
# Import libraries

import re
import time



#*************************************************************************


def entry_format(file):

    """Format the GenBank file by splitting it into separate entries, removing
    unnecessary characters, and changing it to lower case.

    Input:  file                --- A file containing multiple entries in
                                standard GenBank format.
    Return: formatted_entries   --- A list of formatted GenBank entries suitable
                                for data extraction.

    25.04.18  Original   By: AO
    """
    
    formatted_entries = []
    with open(file, 'rt') as f:
        data = f.read() # read the whole document into memory
        entries = data.split("\n//") # split at the end of an entry (each finish with //)
        # format each entry, e.g. all information was stored in the database in
        # lower case, to simplify querying due to a case-sensitivity in Linux
        for d in entries:
            clean = d.replace('http://', '').replace('\n', ' ').lstrip().lower() + "//" #last // added in order to show the end of the line for regex
            if clean.startswith("locus"):
                formatted_entries.append(clean)
    
    return(formatted_entries)



#*************************************************************************


class Locus:

    """Class for handling locus data of the entry. Attributes hold information
    on locus name, cytogenetic location, chromosome name, source, and whole
    locus sequence.

    Input:  entry                 --- a formatted GenBank entry; it is used to
                                  extract information into class attributes
            id_number             --- locus ID number that will be used to load
                                  it into a database
    
    25.04.18 Original By: AO
    """

    def __init__(self, entry, id_number):
        self.raw_data=entry
        self.id = id_number

        p_locus_name = re.compile(r'locus\b\s*(\w*)')
        match = p_locus_name.search(entry)
        if match:
            self.locus_name = match.group(1)
        else:
            self.locus_name = "" # in most of following regexes no match will
            # produce to an empty string, rather than stop the program from
            # working, since some information might be missing for eaxh entry.


        p_chr_location = re.compile(r'/map="(.*?)"')
        match = p_chr_location.search(entry)
        if match:
            self.chr_location = match.group(1)
        else:
            self.chr_location = ""


        p_chr_name = re.compile(r'/chromosome="(\w*)"')
        match = p_chr_name.search(entry)
        if match:
            self.chr_name = "chromosome " + match.group(1)
        else:
            self.chr_name = "chromosome 16"  # since the group assessment is about chromosome 16
            # in reality, no match would indicate "unknown chromosome"


        p_source = re.compile(r'/organism="([a-zA-Z0-9_\s]*)"')
        match = p_source.search(entry)
        if match:
            unformatted_string = match.group(1)
            # in every case when the extracted information might be composed of
            # several words, there is a chance that the words are split by a tab
            # due to the nature of the file. To avoid having large gaps between
            # the words, split the frase and re-join with a single space between them
            string_list = unformatted_string.split()
            formatted_string = ' '.join(string_list)
            self.source = formatted_string
        else:
            self.source = ""


        p_whole_seq = re.compile(r'origin\b\s*(.*?)//')
        match = p_whole_seq.search(entry)
        # store a whole sequence as a single string without numbering and spaces
        if match:
            seq = match.group(1)
            # remove numbers
            seq_no_digits = ''.join(i for i in seq if not i.isdigit())
            # remove all spaces
            self.whole_seq = seq_no_digits.replace(" ", "")
        else:
            self.whole_seq = ""



#*************************************************************************


def accession_extractor(entry):

    """Extract locus accession numbers for a particular locus. If the 
    locus has multiple accession numbers, distinguishes the primary 
    accession number (latest update version) from older accessions.

    Input:  entry               --- An entries in standard GenBank format.
    Return: output              --- A dictionary, where the keys are the
                                accession numbers and the values indicate
                                if the number is the latest version.

    25.04.18  Original   By: AO
    """
    
    p_accession = re.compile(r'accession\b\s*(.*?)\s*version\b')
    match = p_accession.search(entry)
    output = {}
    if match:
        string = match.group(1)
        # an entry might have multiple accession numbers, so need to extract
        # separate words.
        sep_ac =re.compile(r'\w+')
        accession_list = sep_ac.findall(string)
        counter = 0
        # The first number on the list is the most recent, so for convenience
        # and query simplicity we treated as a primary number and designated
        # as the latest version.
        latest_version = 'T'
        for i in accession_list:
            output[i] = latest_version
            counter +=1
            if counter != 0:
                latest_version = 'F'
            # The output is a dictionary with the accession number as a key and
            # T/F as a value.
            
    else:
        # if an entry is missing an accession number, take locus number instead,
        # because they often coincide.
        p_locus_name = re.compile(r'locus\b\s*(\w*)')
        match = p_locus_name.search(entry)
        if match:
            accession = match.group(1)
            latest_version = 'T'
            output[accession]=latest_version
        else:
            # if even locus number is missing
            accession = "unknown " + str(self.locus_id)
            latest_version = ''
            output[accession]=latest_version
    return(output)



#*************************************************************************


class Accession:
    """Class for handling accession data of the entry. Attributes hold
    information on the accession number, corresponding locus ID number and
    whether the accession number is the latest version.

    25.04.18   Original   By: AO
    """

    def __init__(self, accession, latest_version, locus_id):
        self.accession = accession
        self.latest_version = latest_version
        self.locus_id = locus_id
        

        
#*************************************************************************


def cds_extractor(entry, output):

    """Extract information on all CDS in a particular locus, as well as a
    locus sequence. 

    Input:  entry               --- An entries in standard GenBank format.
            output              --- A desired output - either CDS information
                                or locus sequence.
    Return: depending on user's desired output, it is either a list with data
    on all CDS in a locus, or locus sequence.

    25.04.18  Original   By: AO
    """

    if output == 'cds':
        p_cds = re.compile(r'cds\b.*?/translation="[a-zA-Z0-9_\s]*"')
        cds = p_cds.findall(entry)
        if cds:
            return(cds)
        else:
            return('')
    if output == 'locus_seq':
        p_whole_seq = re.compile(r'origin\b\s*(.*?)//')
        match = p_whole_seq.search(entry)
        if match:
            whole_seq = match.group(1)
            seq_no_digits = ''.join(i for i in whole_seq if not i.isdigit())
            format_seq = seq_no_digits.replace(" ", "")
            return(format_seq)
        else:
            return('')

       

#*************************************************************************    


class Cds:
    
    """Class for handling CDS data of the entry. The attributes and methods
    extract information on gene name, product name, protein ID, DNA sequence 
    encoding the gene, translated amino acid sequence, if the sequence is on a
    reverse strand (complement), if it is alternative splicing version, if the
    product is unknown and if the DNA sequence is correct length.

    25.04.18  Original   By: AO
    """

    def __init__(self, entry_cds, locus_seq):
        self.raw_data=entry_cds
        self.locus_seq=locus_seq


        p_gene_name= re.compile(r'/gene="(.*?)"')
        match = p_gene_name.search(entry_cds)
        if match:
            unformatted_string = match.group(1)
            string_list = unformatted_string.split()
            formatted_string = ' '.join(string_list)
            self.gene_name = formatted_string
        else:
            self.gene_name = ""

        p_product_name= re.compile(r'/product="(.*?)"')
        match = p_product_name.search(entry_cds)
        if match:
            unformatted_string = match.group(1)
            string_list = unformatted_string.split()
            formatted_string = ' '.join(string_list)
            self.product_name = formatted_string
        else:
            self.product_name = ""

        p_protein_id= re.compile(r'/protein_id="(.*?)"')
        match = p_protein_id.search(entry_cds)
        if match:
            unformatted_string = match.group(1)
            string_list = unformatted_string.split()
            formatted_string = ' '.join(string_list)
            self.protein_id = formatted_string
        else:
            self.protein_id = ""

        p_seq_location= re.compile(r'cds\b\s*.*join\((.*?)\)')
        match = p_seq_location.search(entry_cds)
        if match:
            seq_loc = match.group(1)
            self.seq_location = seq_loc.replace(" ", "")
        else:
            self.seq_location = ""

        p_translation = re.compile(r'/translation="([a-zA-Z\s]*)"')
        match = p_translation.search(entry_cds)
        if match:
            transl = match.group(1)
            self.translation = transl.replace(" ", "")
        else:
            self.translation = ""

        p_complement = re.compile(r'cds\b\s*complement\(')
        # if the CDS sequence is a reverse complement, the CDS sequence will
        # have to be constructed accordingly.
        match=p_complement.search(entry_cds)
        if match:
            self.complement='T'
        else:
            self.complement='F'

        p_alt_splicing = re.compile(r'alternative')
        # Detects if the CDS has alternative splicing. However, in some cases
        # it is stated for all variants, but in others only for actual alternative
        # variants, so it will be dealt with when processing the whole entry.
        match=p_alt_splicing.search(entry_cds)
        if match:
            self.alternative_splicing='T'
        else:
            self.alternative_splicing='F'


        p_unknown_product= re.compile(r'/product="unknown')
        match = p_unknown_product.search(entry_cds)
        if match:
            self.unknown_product = 'T'
        else:
            self.unknown_product = 'F'


    def cds_seq(self):
        """Method generates the cds sequence using the sequence of the locus
        and cds location"""

        sequence = ""
        #extract coordinate pairs
        p_coordinate_pairs = re.compile(r'(\d+)\.\.(\d+)')
        coordinate_pairs = p_coordinate_pairs.finditer(self.seq_location)
        for pair in coordinate_pairs:
            start = int(pair.group(1))-1
            end = int(pair.group(2))
            # extract the sequence for every coordinate pair and join to the
            # final sequence
            fragment = self.locus_seq[start:end]
            sequence = sequence + fragment

        if self.complement=='T':
            # in case if the sequence is a complement, reverse it, split it to
            # make a list and build a new list from complementary bases and
            # re-join the sequence.
            reverse_seq = sequence[::-1]
            seq_list = list(reverse_seq)
            complement_list = []
            for i in seq_list:
                if i=='a':
                    complement_list.append('t')
                elif i=='t':
                    complement_list.append('a')
                elif i=='g':
                    complement_list.append('c')
                elif i=='c':
                    complement_list.append('g')
            reverse_complement = "".join(complement_list)
            return(reverse_complement)
        else:
            return(sequence)
    


    def cds_length_correct(self):
        """Many entries contain sequences that need to be joined from different
        loci to form a product. The method checks it if the nucleotide
        sequence length corresponds to the length of protein sequence."""
        protein_length = (len(self.cds_seq())/3)-1 #subtract 1 because translation
        # sequence does not contain stop codon
        if protein_length == len(self.translation):
            return('T')
        else:
            return('F')
        
            

#*************************************************************************


def parser(input_file, output_file, species, chromosome):

    """The function uses the GenBank file to build INSERT statements
    suitable for uploading the information into the database created for
    the chromosome browser project. The complete set of statements fills in
    Chromosome, Locus, Accession and CDS tables.

    Input:  input_file          --- The original GenBank file with all loci
                                in the chromosome.
            output_file         --- The file where the INSERT statements
                                will be written.
            species             --- The source of genetic information (species).
            chromosome          --- Chromosome contained in GenBank file.
    Return: INSERT statements for loci, corresponding CDS and accession numbers
    are written into user's file of choice.

    25.04.18  Original   By: AO
    """
    
    start = time.time()
    # First, format the document with raw data to produce a list of entries.
    entry_list = entry_format(input_file)
    print("File (%s) has been loaded in %s" % (input_file, (time.time() - start)))
    print('Total records loaded: ', len(entry_list))

    
    # Tha database design requires chromosome table to be filled in first
    with open(output_file, "a") as db_entries:
        db_entries.write('INSERT INTO chromosome(name, source)\nVALUES ("' + chromosome + '", "' + species +'");\n')

    # The counter will correspond to locus ID.
    counter = 1
    start = time.time()
    all_entry_sql = []

    # Process every entry in the list.
    for entry in entry_list:
        entry_sql = process_entry(entry, counter)
        if entry_sql:
            all_entry_sql.append(entry_sql)
            counter += 1

    print("\nEntry processing finished in %s" % (time.time() - start))
    entry_count = int(counter)-1
    print("Number of valid entries: %s" % entry_count)

    # Write each processed entry into an output document.
    with open(output_file, "a") as db_entries:
        for unit in all_entry_sql:
            db_entries.write(unit)



#*************************************************************************


def process_entry(entry, counter):

    """The function extracts all relevant locus, CDS and accession information
    from a single GenBank entry and builds all corresponding INSERT statements. 

    Input:  entry               --- The entry in a standard GenBank format.
            counter             --- The number that corresponds to locus ID 
    Return: a string with all INSERT statements for a locus, its CDS and
    accession numbers.

    25.04.18  Original   By: AO
    """
    
    locus_seq = cds_extractor(entry, 'locus_seq')
    cds_list = cds_extractor(entry, 'cds')
    cds_insert_strings = ""
    gene_list = []
    cds_counter = 0
    for i in cds_list:
        cds = Cds(i,locus_seq)
        if cds.cds_length_correct()=='F':
            pass # the assessment description suggests to remove all CDS, which
        # products have to be combined from multiple loci.
        else:
            if cds.gene_name not in gene_list:
                # hence it is the first time the gene appears,
                # even it has alternative splicing appears
                cds_counter += 1
                cds_insert_strings = cds_insert_strings + 'INSERT INTO cds(gene_name, product_name, product_id, seq_location, whole_seq, translation, complement, locus_id)\nVALUES ("' + str(cds.gene_name) + '", "' + str(cds.product_name) + '", "' + str(cds.protein_id) + '", "' + str(cds.seq_location) + '", "' + str(cds.cds_seq()) + '", "' + str(cds.translation) + '", "' + str(cds.complement) + '", '+ str(counter) + ');\n'
                gene_list.append(cds.gene_name)
            elif cds.gene_name in gene_list:
                if cds.alternative_splicing=='T':
                    pass # this way only the first version of a transcript is
                # captured, but all subsequent alternative splicing is removed.
                if cds.alternative_splicing=='F' and cds.unknown_product=='T':
                    cds_counter += 1
                    cds_insert_strings = cds_insert_strings + 'INSERT INTO cds(gene_name, product_name, product_id, seq_location, whole_seq, translation, complement, locus_id)\nVALUES ("' + str(cds.gene_name) + '", "' + str(cds.product_name) + '", "' + str(cds.protein_id) + '", "' + str(cds.seq_location) + '", "' + str(cds.cds_seq()) + '", "' + str(cds.translation) + '", "' + str(cds.complement) + '", '+ str(counter) + ');\n'
                    # this accounts for cases when no genes in the locus are
                    # known/identified (since many entries lack gene names), and
                    # yet the transcripts are not alternative splicings of the
                    # same gene.
                else:
                    pass

                    
    if cds_counter == 0:
        return
    # Assessment description mentioned to keep only loci with (valid) CDS.
    else:
        accession_dict = accession_extractor(entry)
        accession_insert_strings = ""
        for k,v in accession_dict.items():
            accession_num = k
            version = v
            accession = Accession(accession_num,version,counter)
            accession_insert_strings = accession_insert_strings + 'INSERT INTO accession(accession_num, locus_id,latest_version)\nVALUES ("' + str(accession.accession) + '", ' + str(accession.locus_id) + ', "' + str(accession.latest_version) +'");\n'


        locus = Locus(entry,counter)
        locus_insert_string = 'INSERT INTO locus(id, whole_seq, chr_location, locus_name, chr_name)\nVALUES ('+str(counter) + ', "' +  str(locus.whole_seq) + '", "' + str(locus.chr_location) + '", "' + str(locus.locus_name)+ '", "' + str(locus.chr_name) + '");\n'
        # Return all INSERT statements for a locus, with all its accession number
        # and all sorresponding CDS information in the INSERT statement format.
        return(locus_insert_string + cds_insert_strings + accession_insert_strings)
    

            
