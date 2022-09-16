#!/usr/bin/python3
# -*- coding: utf-8 -*-
#

import os, sqlite3, subprocess, re, time, requests
from Bio import Entrez
from Bio import SeqIO
from ete3 import NCBITaxa

# This is the css style and javascript section for html documents to ensure a uniform display
def get_script_style_section():
    script_style_string = '''<script>
function myFunction(id) {
  var x = document.getElementById(id);
  if (x.style.display === "none") {
    x.style.display = "block";
  } else {
    x.style.display = "none";
  }
}
</script>
<style>
body {
    background-color: white;
}
h1 {
    color: ;
}
p {
    color: DarkGray;
}
thead, tfoot {
    background-color: #3f87a6;
    color: #fff;
}
tbody {
    background-color: #e4f0f5;
}
caption {
    padding: 10px;
    caption-side: bottom;
}
table {
    border-collapse: collapse;
    border: 2px solid rgb(200, 200, 200);
    letter-spacing: 1px;
    font-family: sans-serif;
    font-size: .8rem;
}
td,
th {
    border: 1px solid rgb(190, 190, 190);
    padding: 5px 10px;
}
td {
    text-align: left;
    vertical-align: top
}
.a1 {
    display: inline-block;
    width: 240px;
}
</style>'''
    return script_style_string

# Sanitize malformed sceintific species names (some people use single quotes, parenthesis, etc. in the scientific name!)
def sanitize_species_name(SPECIES_NAME):
    # Replace single quotes, commas, and round parenthesis by spaces 
    translation_table = str.maketrans('\',() ', '_____')
    return SPECIES_NAME.translate(translation_table) 

# Download a single fasta file
def download_fasta_file(ID_OR_ACCNO, FASTA_FILE):
    Entrez.email = "michael@jeltsch.org"
    Entrez.tool = "local_script_under_development"
    with open(FASTA_FILE, 'w') as filehandle:
        print('Retrieving sequence {0} from Entrez...'.format(ID_OR_ACCNO))
        time.sleep(1)
        try:
            with Entrez.efetch(db="protein", rettype="fasta", retmode="text", id=ID_OR_ACCNO) as seqhandle:
                seq_record = SeqIO.read(seqhandle, "fasta")
                filehandle.write(seq_record.format("fasta"))
        except Exception as err:
            print('Problem contacting Blast server. Skipping {0}. Error: {1}.'.format(ID_OR_ACCNO, err))
            return False
        else:
            return True


# Get all fasta sequences from the master_dictionary via Entrez
def download_proteins(target_dir, master_dict, rename_fasta_description_after_key=False, overwrite=False, filename='after_key'):
    Entrez.email = "michael@jeltsch.org"
    Entrez.tool = "local_script_under_development"
    for key, value in master_dict.items():
        if filename == 'after_key':
            FASTA_FILE = '{0}/{1}.fasta'.format(target_dir, key)
        elif filename == 'after_value[0]':
            FASTA_FILE = '{0}/{1}.fasta'.format(target_dir, value[0])
        if os.path.isfile(FASTA_FILE) and os.path.getsize(FASTA_FILE) > 0:
            print('Not overwriting existing non-zero-size file {0}.'.format(FASTA_FILE))
        else:
            with open(FASTA_FILE, 'w') as filehandle:
                print('Retrieving sequence {0} ({1}) from Entrez...'.format(value[0], key))
                try:
                    with Entrez.efetch(db="protein", rettype="fasta", retmode="text", id=value[0]) as seqhandle:
                        seq_record = SeqIO.read(seqhandle, "fasta")
                        if rename_fasta_description_after_key == True:
                            seq_record.id = key
                        filehandle.write(seq_record.format("fasta"))
                except Exception as err:
                    print('Problem contacting Blast server. Skipping {0} - {1}. Error: {2}.'.format(key, value[0], err))
            time.sleep(2)

# This takes complex taxon data and returns a tuple (a string of the main taxon and a list of all taxons to be subtracted)
# E.g. '8504-1329911' -> 8504, [1329911]
def expand_complex_taxa(taxon_data):
    print('taxon_data to expand: {0}'.format(taxon_data))
    TAXA_LIST = []
    try:
        print('Expanding complex taxon data to a list: {0} -> '.format(taxon_data), end = '')
        try:
            taxon_list = re.split('([-\+])',taxon_data)
            print('taxon_list: {0}'.format(taxon_list))
        except Exception as err:
            taxon_list = [taxon_data]
        # if taxon list starts with a taxon name without prepended +/-, handle the first element separately
        if taxon_list[0] not in '+-':
            TAXA_LIST.append('+'+taxon_list[0])
            taxon_list.pop(0)
        # Now iterate through the list
        for item in taxon_list:
            if item in '+-':
                sign = item
            else:
                if sign == '+':
                    TAXA_LIST.append('+'+item)
                elif sign == '-':
                    TAXA_LIST.append('-'+item)
                else:
                    print('Formating error in taxon_id definition.')
    except Exception as err:
        print('Nothing to expand from {0} or error expanding. Error: {1}'.format(taxon_data, err))
    print('taxa_list: {0}'.format(TAXA_LIST))
    return TAXA_LIST

def load_blacklist():
    # All the phyla in the list below have been found in the first screen not to have any VEGF-like molecules
    #blacklist = ['ctenophora', 'porifera', 'placozoa', 'xenacoelomorpha', 'cyclostomata', 'onychophora', 'pycnogonida', 'myriapoda', 'nematomorpha', 'loricifera', 'kinorhyncha', 'chaetognatha', 'bryozoa', 'entoprocta', 'cycliophora', 'nemertea', 'phoroniformea', 'gastrotricha', 'platyhelminthes', 'gnathostomulida', 'micrognathozoa', 'orthonectida', 'dicyemida']
    # All hits for the phyla below were manually checked and changed in taxon_data.py => therefore, we do not want to chenge them anymore
    #blacklist = ['porifera', 'xenacoelomorpha', 'myriapoda', 'bryozoa', 'entoprocta', 'platyhelminthes', 'orthonectida']
    #blacklist = ['kinorhyncha']
    blacklist = []
    return blacklist

def blast_formatter(RID, FORMAT, FILE, VERBOSE = True):
    # The blast formatter requires access to the same database for the reformatting since some
    # data needed for the xml file is not contained in the html file, but is retrieved during formatting.
    # In this case a local database must either be available to fetch the data or, alternatively,
    # if the RID is known, you can request the blast_reformater to fetch you a refomatted
    # file from the blast server. The RID can be found from the html file.
    comment = 'Converting html into xml:'
    #bash_command = 'blast_formatter -rid {0} -outfmt 5 > {1}'.format(RID, BLAST_XMLFILE)
    #print('comment: {0}\nbash_command: {1}'.format(comment, bash_command))
    #output, error = execute_subprocess(comment, bash_command)
    bash_command_as_list = ['blast_formatter', '-rid', RID, '-outfmt', '5']
    print('{0} {1}'.format(comment, bash_command_as_list))
    working_directory = '.'
    try:
        result = subprocess.run(bash_command_as_list, stdout = subprocess.PIPE, stderr = subprocess.PIPE, cwd = working_directory)
    except Exception as ex:
        if VERBOSE: print('Subprocess error: {0}'.format(result.stderr))
        return result.stderr.decode('utf-8')
    else:
        if VERBOSE: print('Subprocess output: {0}'.format(result.stdout))
        #if result.stdout == 'Network error: Error fetching sequence data from BLAST databases at NCBI, please try again later'
        # Write blast results to an html file
        with open(FILE, "w") as out_handle:
            out_handle.write(result.stdout.decode('utf-8'))
            print('New format {0} written to file {1}.'.format(FORMAT, FILE))
        out_handle.close()
        return result.stdout.decode('utf-8')

def execute_subprocess(comment, bash_command, working_directory='.'):
    print('\nSubprocess: {0} {1} {2}'.format(comment, bash_command, time.ctime()))
    process = subprocess.Popen(bash_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, cwd=working_directory)
    output, error = process.communicate()
    process_status = process.wait()
    output = output.decode('utf-8')
    error = error.decode('utf-8')
    if output != '':
        print('Subprocess output: {0}'.format(output))
    if error != '':
        print('Subprocess error: {0}'.format(error))
    print('Subprocess output/error: {0}/{1}\n{2}'.format(output, error, time.ctime()))
    return output, error

def execute_subprocess_new(comment, bash_command_as_list, working_directory='.'):
    print('\nSubprocess function entered.')
    print('{0} {1}'.format(comment, bash_command_as_list))
    print('start time: {0}'.format(time.ctime()))
    result = subprocess.run(bash_command_as_list, stdout = subprocess.PIPE, stderr = subprocess.PIPE, cwd = working_directory)
    print('Subprocess output: {0}'.format(result.stdout))
    print('Subprocess error: {0}'.format(result.stderr))
    print('completion time: {0}'.format(time.ctime()))
    return result

def make_synonym_dictionary(master_dictionary):
    #print('master_dictionary:\n{0}'.format(master_dictionary))
    synonym_dictionary = {}
    for key, value in master_dictionary.items():
        #print('key:\n{0}'.format(key))
        #print('value:\n{0}'.format(value))
        #print('value[3][0]:\n{0}'.format(value[3][0]))
        canonical_ortholog_group_name = value[3][0]
        print('canonical_ortholog_group_name: {0}'.format(canonical_ortholog_group_name))
        for item in value[3]:
            synonym_dictionary[item] = canonical_ortholog_group_name
    return synonym_dictionary

def load_dictionary(FILENAME, VERBOSE=True):
    #print("Enter subroutine")
    # Load a sequence_dictionary if it exists
    if os.path.isfile(FILENAME):
        try:
            preamble, dictionary = read_file_to_dict(FILENAME)
            print("\nReading in the dictionary " + FILENAME + ":\n")
            if VERBOSE == True:
                for key, value in dictionary.items():
                    print(key, value)
            print("Done reading dictionary.")
        except Exception as e:
            preamble = '#'
            dictionary = {}
            print('Could not read taxon dictionary {0}'.format(FILENAME))
    else:
        preamble = '#'
        dictionary = {}
    return preamble, dictionary

def insert_line_breaks(file_name):
    try:
        with open(file_name, "r") as file:
            content = file.read()
            file.close()
        content = content.replace("}], '","}],\n     '")
        content = content.replace("], '","],\n '")
        with open(file_name, "w") as file:
            file.write(content)
            file.close()
            return True
    except Exception as ex:
        print('Could not insert line breaks into file {0} Error: {1}'.format(filename, str(ex)))
        return False

def write_dict_to_file(preamble, dictionary, file_name):
    try:
        with open(file_name, "w") as file:
            file.write(preamble)
            file.write(str(dictionary))
            file.close()
            return True
    except Exception as ex:
        print("Could not write dictionary to file. Error: " + str(ex))
        return False

def read_file_to_dict(file_name):
    try:
        with open(file_name, "r") as file:
            string = file.read()
            file.close()
            #print('File content:\n' + string)
            dictionary = eval(string)
            # Store the first three lines of the dictionary
            buf = string.split('\n')
            preamble = ''
            for i in range(0,len(buf)):
                if buf[i][:1] == '#':
                    preamble += buf[i] + '\n'
                else:
                    break
            return preamble, dictionary
    except Exception as ex:
        print('Could not read dictionary from file {0}. Error: '.format(file_name) + str(ex))
        return '', {}

def create_sqlite_file(FILE_NAME):
    try:
        conn = sqlite3.connect(FILE_NAME)
        print(sqlite3.version)
        commands = ['''CREATE TABLE IF NOT EXISTS `species` (`scientific_name` TEXT NOT NULL, `taxon_id` INTEGER NOT NULL, `phylum` TEXT NOT NULL, PRIMARY KEY(`taxon_id`))''',
        '''CREATE TABLE IF NOT EXISTS `protein` (`id` VARCHAR NOT NULL, `accession_no` TEXT, `fasta_description` TEXT NOT NULL, `species` TEXT NOT NULL, `ortholog_group` TEXT, `curated_manually_by` TEXT, PRIMARY KEY(`id`))''',
        '''CREATE TABLE IF NOT EXISTS `ortholog_groups` (`ortholog_group` TEXT NOT NULL)''',
        '''CREATE TABLE IF NOT EXISTS `curator` (`curator` TEXT NOT NULL)''']
        cur = conn.cursor()
        for command in commands:
            print('Executing command: {0}'.format(command))
            cur.execute(command)
            conn.commit()
    except Exception as err:
        print(err)
        return False
    else:
        return True
    finally:
        conn.close()

# Converts the number of seconds into a days/hours/minutes/seconds string
def execution_time_str(elapsed_time_seconds):
    min, sec = divmod(elapsed_time_seconds, 60)
    hours, min = divmod(min, 60)
    days, hours = divmod(hours, 24)
    return (str(days) + " days " if days != 0 else '') + (str(hours) + " hours " if hours != 0 else '') + str(min)[:-2] + " min " + str(round(sec, 1)) + " sec"

# This should be replaced by the corresponding function from ETE3
def get_taxon_id_from_NCBI(species_name, VERBOSE=True):
    ncbi = NCBITaxa()
    taxon_id_list = ncbi.get_name_translator([species_name])
    try:
        taxon_id = taxon_id_list[species_name][0]
        # Check that something sensinble was returned
        if isinstance(taxon_id, int) and taxon_id > 0:
            print('taxon_id for {0}: {1}'.format(species_name, taxon_id))
        else:
            taxon_id = None
            print('No taxon_id for {0}. NCBI returned: {1}'.format(species_name, taxon_id_list))
    except:
        taxon_id = None
        print('No taxon_id for {0}. NCBI returned: {1}'.format(species_name, taxon_id_list))
    #return TAXON_ID
    return taxon_id

# This sends a species id to NCBI and returns to which of the phyla within taxon_dictionary.py
# this species belongs to.
def get_phylum_from_NCBI(TAXON_ID_OR_NAME, VERBOSE=True):
    preamble1, species_dictionary = load_dictionary('prog_data/species_data.py', VERBOSE)
    # Check first the locally cached data
    if TAXON_ID_OR_NAME in species_dictionary and species_dictionary[TAXON_ID_OR_NAME] != 'unknown':
        phylum = species_dictionary[TAXON_ID_OR_NAME]
        print('Species {0} assigned from local cache to phylum {1}.'.format(TAXON_ID_OR_NAME, phylum))
    else:
        print('Species {0} not found in local cache. Asking NCBI...'.format(TAXON_ID_OR_NAME))
        preamble2, taxon_dictionary = load_dictionary('prog_data/taxon_data.py', VERBOSE)
        # If a latin species name is given, translate it into a taxonomic id
        if isinstance(TAXON_ID_OR_NAME, int):
            TAXON_ID = str(TAXON_ID_OR_NAME)
        elif TAXON_ID_OR_NAME.isdigit():
            TAXON_ID = TAXON_ID_OR_NAME
        else:
            TAXON_ID = get_taxon_id_from_NCBI(TAXON_ID_OR_NAME)
        # Counter to increase waiting period upon incresing number of unsuccesful trials
        i = 0
        while True:
            try:
                URL = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id={0}&retmode=xml&rettype=full'.format(TAXON_ID)
                r = requests.get(URL)
                text = r.text
                if VERBOSE: print('request: {0}\nresponse: {1}'.format(URL, text))
                i = 0
                for taxon, taxon_data in taxon_dictionary.items():
                    # Determine, whether a taxon is complex (with some branches to exclude)
                    try:
                        TAXON = taxon_data[0].split('-')[0]
                        SUBTRACT_TAXON = taxon_data[0].split('-')[1]
                    except Exception as err:
                        SUBTRACT_TAXON = '-'
                    if VERBOSE: print('Looking for TAXON_ID {0}, exluding TAXON_ID {1}'.format(TAXON, SUBTRACT_TAXON))
                    if VERBOSE:
                        if TAXON == '8782':
                            print('TAXON = 8782 = aves')
                    if '<TaxId>{0}</TaxId>'.format(TAXON) in text and '<TaxId>{0}</TaxId>'.format(SUBTRACT_TAXON) not in text:
                        if VERBOSE: print('Adding phylum info {0} to species {1}'.format(taxon, TAXON_ID))
                        phylum = taxon
                        print(phylum)
                        break
                    i += 1
                # The following line will be only executed if calling 'phylum' does not give a NameError
                # When i = 51, none of the 51 phyla was in the text that was received from NCBI.
                if i == 51:
                    if VERBOSE: print('Species {0} does not belong to any phylum in the phylum list'.format(TAXON_ID))
                    phylum = 'unknown'
                elif i < 51:
                    i += 1
                    if VERBOSE: print('{0} phyla (of 51) checked before hit was found.'.format(i))
                else:
                    if VERBOSE: print('Unknown error when trying to identify phylum for species {0}.'.format(TAXON_ID))
                    phylum = 'unknown'
            except NameError as err:
                if VERBOSE: print('Sleeping due to server error. Will retry in 60 seconds. Error: {0}'.format(err))
                time.sleep(60)
                continue
            # If there were no exceptions (= we got the phylum), break out of the while loop
            else:
                break
        # Add result from NCBI to dictionary, write new dictionary file, and replace the old dictioanry file with the new one
        try:
            species_dictionary[TAXON_ID_OR_NAME] = phylum
            write_dict_to_file(preamble1, species_dictionary, 'prog_data/.species_data.py')
            os.rename('prog_data/.species_data.py', 'prog_data/species_data.py')
        except:
            print('Could not write phylum to prog_data/species_data.py')            
    return phylum
