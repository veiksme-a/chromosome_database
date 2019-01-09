#!/usr/bin/env python3

"""
Program:    restriction_site_calc
File:       restriction_site_calc.py

Version:    V1.0
Date:       28.04.18
Function:   Calculate the restriction sites for each of the preset restriction
            enzymes EcoRI, BamHI and BsuMI in each of the locus sequences in
            the database, and produce INSERT statements for bulk loading into
            the database.

Author:     AO
            
--------------------------------------------------------------------------

This program was produced as a part of Bioinformatics MSc programme assessment.

--------------------------------------------------------------------------
Description:
============
A set of functions that extract locus sequences and finds restriction sites
for each of the preset restriction enzymes. The output is a set of INSERT
statements to fill in the restriction_enzyme and restriction_sites tables
of the chromosome database.

--------------------------------------------------------------------------
Usage:
======
To obtain INSERT statements, the user needs to execute restr_site_parser()
function, providing an output file where the statements will be written.

--------------------------------------------------------------------------
Revision History:
=================
V1.0   28.04.18   Original   By: AO
"""

#*************************************************************************

#Import libraries
import pymysql.cursors
import re



#*************************************************************************


def restr_enzymes(output_file):

    """Write INSERT statements for the preset restriction enzymes and return a
    dictionary with preset enzymes, their IDs and recognition sequences for
    further use.

    Input:  output_file          --- An file where the user wants the statements
                                 to be written into.
    Return: preset_restr_enzymes --- A dictionary, where the keys are 
                                enzyme IDs and values are enzyme names and
                                their recognition sequences.

    28.04.18  Original   By: AO
    """
    
    # a dictionary of preset restriction enzymes with their recognition sequences
    preset_restr_enzymes = {1:['EcoRI', 'gaattc'], 2:['BamHI', 'ggatcc'], 3:['BsuMI','ctcgag']}

    for k,v in preset_restr_enzymes.items():
        # first, write the INSERT statements for the enzymes (so call this
        # first and use it once to avoid duplication
        with open(output_file, "a") as db_entries:
            db_entries.write('INSERT INTO restriction_enzyme(id,name, recogn_seq)\nVALUES (' + str(k) + ', "' + str(v[0]) + '", "' + str(v[1]) + '");\n')
    # return the restriction enzyme dictionary for further use
    return(preset_restr_enzymes)



#*************************************************************************


def locus_seq():

    """Retrieve all loci sequences and their IDs.

    Return: output              --- A list of dictionaries, where the keys are 
                                the column names from the database and values
                                contain loci IDs and sequences.

    28.04.18  Original   By: AO
    """
    
    query = 'select id, whole_seq from chr_db.locus;'
    db = pymysql.connect(db='chr_db', user='root', passwd='', host='localhost', port=3306, cursorclass = pymysql.cursors.DictCursor)

    cursor = db.cursor()

    q = cursor.execute(query)

    output = cursor.fetchall()

    return(output)



#*************************************************************************


class RestrictionSites:
    """ A class for restriction site information. Attributes contain locus ID,
    restriction enzyme ID and a list of start end end positions for each enzyme
    restriction site.
    28.04.18    Original    By: AO
    """

    def __init__(self, locus_id, re_id, locus_seq, recognition_seq):
        self.locus_id = locus_id
        self.re_id = re_id

        # produce a position list containing all restriction sites
        position_list = []
        # construct an iterator by searching for enzyme recognition
        # sequence within locus sequence
        all_matches = re.finditer(re.escape(recognition_seq), locus_seq)

        if all_matches:
            # if a particular enzyme has restriction sites present,
            # calculate the site span in the sequence
            for match in all_matches:
                start_end_pos = match.span()
                # produces a tuple with start and end positions of the site 
                position_list.append(start_end_pos)
        else:
            pass

        self.position_list = position_list



#*************************************************************************


def restr_site_parser(output_file):
    
    """Parser creates INSERT statements for each preset restriction enzyme and
    for each restriction site in each locus.

    Input:  output_file         --- A file where the user wants INSERT statements
                                to be written to.
    Return: Writes INSERT statements into a file of choice.

    27.04.18  Original   By: AO
    """

    enzymes = restr_enzymes(output_file)
    locus_data = locus_seq()
    # Write an INSERT statements for each locus [for every enzyme in preset list)
    for locus in locus_data:
        for enz_id,enzyme in enzymes.items():
            # first, calculate all restriction sites using RestrictionSites class
            sites = RestrictionSites(locus['id'], enz_id, locus['whole_seq'], enzyme[1])
            
            if not sites.position_list:
                pass
            else:
                # if a particular enzyme has at least one restriction site,
                # compile an INSERT statement for each restriction site
                for site in sites.position_list:
                    insert = 'INSERT INTO restriction_sites(locus_id, re_id, start_position, end_position)\nVALUES (' + str(sites.locus_id) + ', ' + str(sites.re_id) + ', ' + str(site[0]) + ', ' + str(site[1]) + ');\n'
                    with open(output_file, "a") as db_entries:
                        db_entries.write(insert)
