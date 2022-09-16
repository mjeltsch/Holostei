#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# This is a first implementation of an algorhythm that takes three lists
# of proteins (VEGF-A, VEGF-C, VEGF-D) and finds the corresponding mRNA
# sequences. The purpose of this is to do a codon analysis (synonymous
# versus non-synomynous mutations) in order to identify the degree of
# diversifying versus purifying evolution (i.e. conservation). There
# are many problems with writing such a script, most notably that the
# protein and gene (prediction) databases are not kept in sync. Therefore,
# we needed to implement the possibility of manual overrides by providing
# local uniprot.xml (= protein sequence) and mRNA.fasta (= mRNA sequence)
# files. The script is still a mess and needs refactoring, but it works.
#   
# Another problem was BioPython, which seems to have a problem with
# parsing uniprot.xml files. We could not extract the links to the Ensembl
# sequences and therefore switched to BeautifulSoup to accomplish this
# task.

import os, re, sys, glob, shutil, gzip, argparse, time, datetime, requests, json
#from ssl import TLSv1
from pickle import FALSE
from Bio import SeqIO
from Bio import Entrez
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastpCommandline
from io import StringIO
#from time import sleep, time
from bs4 import BeautifulSoup
import lxml
from phylolib import execution_time_str
import xml.etree.ElementTree as ET

def extract_organism(seq_record, VERBOSE = True):
    # get organism from square brackets (if there are such)
    re_search_result = re.search('\[(.*)\]', str(seq_record))
    if re_search_result:
        #if VERBOSE: print('Protein sequence detected.')
        organism = re_search_result.group(0)[1:-1]
        # get rid of three-word organisms
        organism = ' '.join(organism.split(' ')[:2])
        if VERBOSE: print('Organism (from protein record {0}): {1}.'.format(seq_record.id, organism))
    # Check whether it is a protein sequence or mRNA sequence
    elif seq_record.id[:3] in ['sp|', 'tr|']:
        #if VERBOSE: print('Protein sequence detected.')
        organism = re.search('OS=(.*) OX=', str(seq_record)).group(1)
        if VERBOSE: print('Organism (from protein record {0}): {1}.'.format(seq_record.id, organism))
    # The rest should be mRNA sequences
    else:
        #if VERBOSE: print('mRNA sequence detected.')
        if 'PREDICTED' in str(seq_record):
            new_str = str(seq_record).split('PREDICTED: ')[1]
        else:
            new_str = seq_record.description.split(' ', 1)[1]
            #print('new_str: {0}'.format(new_str))
        organism = ' '.join(new_str.split(' ', 2)[:2])
        if VERBOSE: print('Organism (from mRNA record {0}): {1}.'.format(seq_record.id, organism), end = '')
    return organism

def tblastn(PROTEINSEQUENCE, ORGANISM):
    BLAST_DATABASE = 'nr'
    EVALUE = 0.01
    BLAST_XMLFILE = 'results/'+ORGANISM+'_mRNA.xml'
    start_time = time.time()
    print('\nRunning remote tblastn requesting one result with the following query sequence: {0}.'.format(ORGANISM))
    print('Blast job started at {0}'.format(str(datetime.datetime.now())[:-7]))
    try:
        result_handle = NCBIWWW.qblast('tblastn', BLAST_DATABASE, PROTEINSEQUENCE, hitlist_size = 1, expect = EVALUE, entrez_query=''+ORGANISM+'[organism]', format_type = 'XML')
        # Write blast result to xml file
        with open(BLAST_XMLFILE, "w") as out_handle:
            out_handle.write(result_handle.read())
            print('{0} written to disk.'.format(BLAST_XMLFILE))
        result_handle.close()
        with open(BLAST_XMLFILE) as blastresult:
            blast_record = NCBIXML.read(blastresult)
            print('Blast records: {0}'.format(len(blast_record.descriptions)))
            if len(blast_record.descriptions) != 1:
                nt_accession_number = False
            else:
                for description in blast_record.descriptions:
                    print(str(description)+"\n")
                    print("\ndescription: " + description.title)
                    for align in blast_record.alignments:
                        for hsp in align.hsps:
                            nt_accession_number = align.hit_id.split('|')[-2]
                            print('Checking nt_accession_number "{0}"'.format(nt_accession_number))
                            # Check whether a valit accession number has been extracted
                            check = nt_accession_number.split('.')
                            if check[0][0].isupper() and check[-1].isdigit():
                                print(check[0][0])
                                print(check[-1])   
                            else:
                                nt_accession_number = False
    except Exception as ex:
        how_long = execution_time_str(time.time()-start_time)
        print('Blasting failed. It took {0}. The error was: '.format(how_long, ex))
        return False
    else:
        how_long = execution_time_str(time.time()-start_time)
        if nt_accession_number == False:
            print('Blasting completed successfully in {0} but not returning any sequences.\n'.format(how_long))
        else:
            print('Blasting completed successfully in {0}, returning sequence {1}\n'.format(how_long, nt_accession_number))
        return nt_accession_number


def write_fasta_to_file(outputfile, sequence):
    print('Writing to file {0}:'.format(outputfile))
    print('Content:\n'+content)
    with open(outputfile, 'w') as output_file:
        output_file.write(str(sequence))

def e_fetch(outputfile, nt_accession_number):
    try:
        print('Trying to e_fetch {0}'.format(nt_accession_number))
        with Entrez.efetch(db="nucleotide", rettype="fasta", retmode="text", id=nt_accession_number) as handle:
            seq_record = SeqIO.read(handle, "fasta")
            content = seq_record.format("fasta")
            write_fasta_to_file(outputfile, content)
    except Exception as ex:
        print("Problem contacting Entrez server. Skipping " + nt_accession_number + ". Error: " + str(ex))

def get_from_ensemble(accession_number):
    print('Getting Ensemble transcript id from Uniprot using accession number {0}'.format(accession_number))
    URL = 'https://rest.uniprot.org/uniprotkb/{0}.xml'.format(accession_number)
    response = requests.get(URL)
    soup_content = BeautifulSoup(response.content, 'lxml')
    find_result = soup_content.findAll('property')
    if len(find_result) > 0:
        for item in find_result:
            if item.get('type') == 'gene ID':
                return item.get('value').split('.')[0]
    else:
        return False

def fetch_transcript(id):
    def fetch_endpoint(server, request, content_type):
        r = requests.get(server+request, headers={ "Accept" : content_type})
        if not r.ok:
            r.raise_for_status()
            sys.exit()
        if content_type == 'application/json':
            return r.json()
        else:
            return r.text
    def fetch_endpoint_POST(server, request, data, content_type):
        r = requests.post(server+request, headers={ "Accept" : content_type}, data=data)
        if not r.ok:
            r.raise_for_status()
            sys.exit()
        if content_type == 'application/json':
            return r.json()
        else:
            return r.text
    transcripts = []
    server = "http://rest.ensembl.org/"
    con = "application/json"
    ext_get_gene = "lookup/id/" + id + "?expand=1;"
    get_gene = fetch_endpoint(server, ext_get_gene, con)
    for transcript in get_gene['Transcript']:
        transcripts.append(transcript['id'])
    data = json.dumps({ "ids" : transcripts })
    ext_sequence = '/sequence/id?type=cdna'
    sequences = fetch_endpoint_POST(server, ext_sequence, data, "text/x-fasta")
    return sequences

def run():
    parser = argparse.ArgumentParser()
    parser.add_argument('inputfiles', nargs='+')
    args = parser.parse_args()
    mRNA_list = []
    for inputfile in args.inputfiles:
        records = list(SeqIO.parse(inputfile, "fasta"))
        if len(records) > 0:
            unsuccessful_list = []
            Entrez.email = "michael@jeltsch.org"
            if not os.path.exists('results'):
                os.makedirs('results')
            for prot_record in records:
                organism = extract_organism(prot_record)
                outputfile = 'results/' + organism.replace(" ", "_")  +'_mRNA.fasta'
                # LOCAL COPY
                if not os.path.isfile(outputfile):
                    accession_number = prot_record.id.split('|')[-2]
                    # E_FETCH
                    nt_accession_number = tblastn(prot_record.seq, organism)
                    if nt_accession_number != False:
                        e_fetch(outputfile, nt_accession_number)
                        print('mRNA sequence for {0} retrieved from Entrez'.format(accession_number))
                    else:
                        # ENSEMBLE
                        transcript_id = get_from_ensemble(accession_number)
                        if transcript_id != None:
                            sequence = fetch_transcript(transcript_id)
                            write_fasta_to_file(outputfile, sequence)
                        else:
                            print('All methods to retrieve mRNA sequence for {0} were unsuccessful'.format(accession_number))
                            unsuccessful_list.append(organism)
                else:
                    print('Using local copy {}'.format(outputfile))
        else:
            print('{} is not a fasta files or is malformed. Ignoring...'.format(inputfile))            
    print('Retrieval attempts for {0} sequences completed.\nFor the following {1} sequences, no mRNA sequences were found:\n{2}'.format(len(records), len(unsuccessful_list), unsuccessful_list))

if __name__ == '__main__':
    run()
