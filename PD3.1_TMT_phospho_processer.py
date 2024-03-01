""" program "PD1.4_TMT_phospho_processer_v2.py"
Process TMT export files from Proteome Discoverer. For proteins
with a minimum number of unique peptides, total reporter ion
intensities are computed.

The MIT License (MIT)

Copyright (c) 2017 Phillip A. Wilmarth and OHSU

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

Direct questions to:
Technology & Research Collaborations, Oregon Health & Science University,
Ph: 503-494-8200, FAX: 503-494-4729, Email: techmgmt@ohsu.edu.

added median reporter ion intensity cutoff, spring 2015 -PW
added deltamass filter, spring 2015 -PW
double checked things and added some additional comments, Nov. 2015, -PW
added reformatted PSM data dump, Nov. 2015, -PW
parses extra data from phospho searches, Feb. 2016, -PW
site positions were off by one (too large), Aug. 2016, -PW

written by Phil Wilmarth 2015, 2024.

3/20/2017 (PW): added flag to skip or not skip unmodified peptides
6/15/2017 (PW): moved missing data inputation to the combined PSMs only; added DB to results files
1/8/2020 (PW): added log file and brief output file

TODO:
rewrite intensity functions to return test value and do the
testing in the calling function/method.
    
"""
# import common modules here
import os
import sys
import time
import operator
import copy
import PAW_lib

# column separator character
SEPARATOR = '\t'

# missing value replacement (for combined PSMs)
MISSING = 150.0

# minimum trimmed average intensity
INTENSITY = 500.0

# q-Value cutoff
QVALUE = 0.01

# plus/minus PPM deltamass window
PPM = 20.0

# skip unmodified peptides?
SKIP_UNMOD = True

def make_header_map(line, separator='\t'):
    """Returns a dictionary of "column header":index
    """
    header_map = {}
    if not line:
        return header_map

    # split line at separator and populate dictionary (remove any double quotes around headers)
    headers = [x.replace('"', '') for x in line.split(separator)]
    for i, header in enumerate(headers):
        header_map[header] = i
        
    return header_map
    # end

def test_median_intensity(intensities, intensity, zero_input):
    """Find median intensity of channels, excluding "empty" channels.
    Return True if median above intensity, False otherwise.
    """
    int_vector = [x for x in intensities if x > zero_input]

    # exclude PSM if too few channels with minimal intensity
    if len(int_vector) < 3:
        return False

    # compute the median
    int_vector.sort()
    length = len(int_vector)
    if int_vector:
        if not length % 2:
            median = (int_vector[length/2] + int_vector[length/2 - 1]) / 2.0
        else:
            median = int_vector[length/2]
    else:
        median = 0.0
        
    # test threshold
    if median >= intensity:
        return True
    else:
        return False
    
def test_average_intensity(intensities, intensity, zero_input):
    """Find average intensity of channels, excluding "empty" channels.
    Return True if average above "intensity", False otherwise.
    """
    int_vector = [x for x in intensities if x > zero_input]

    # exclude PSM if too few channels with minimal intensity
    if len(int_vector) < 3:
        return False

    # compute the average
    average = sum(int_vector)/float(len(int_vector))
        
    # test threshold
    if average >= intensity:
        return True
    else:
        return False
    
def test_trimmed_average_intensity(intensities, intensity):
    """Parse line and get intensities. Find average intensity of channels,
    excluding the top and bottom values. It does not skip any low values.
    Return True if average above "intensity", False otherwise.
    """
    int_vector = intensities

    # compute the trimmed average
    average = sum(int_vector[1:-1])/float(len(int_vector[1:-1]))
        
    # test threshold
    if average >= intensity:
        return True
    else:
        return False
    
class PSM():
    """Container for PSM information
    """
    def __init__(self, line, header_map, proteins, prot_index, missing=0.0, separator='\t'):
        """Parse line and populate fields
        """
        items = [x.replace('"', '') for x in line.split(separator)]  # PD exports fields may be enclosed in double quotes
        
    
        
        self.channels = ['126C', '127N', '127C', '128N', '128C', '129N',
                         '129C', '130N', '130C', '131N', '131C', '132N',
                         '132C', '133N', '133C', '134N', '134C', '135N']

        # parse row items into object attributes
        # 
        self.checked = items[header_map['Checked']]
        self.confidence = items[header_map['Confidence']]
        self.identifying_node = items[header_map['Identifying Node']]
        self.psm_ambiguity = items[header_map['PSM Ambiguity']]
        self.annotated_sequence = items[header_map['Annotated Sequence']]
        self.sequence = self.annotated_sequence.split('.')[1]
        self.peptide_length = len(self.annotated_sequence.split('.')[1])    # compute peptide length 
        self.modifications = items[header_map['Modifications']] 
        self.number_phospho_groups = self.modifications.count('Phospho')   # add count of phospho sites
        self.number_proteins = self._int(items[header_map['# Proteins']]) # convert to integer
        self.master_protein_accessions = items[header_map['Master Protein Accessions']]
        self.protein_accessions = items[header_map['Protein Accessions']]
        self.protein_descriptions = lookup_headers(self.protein_accessions, proteins, prot_index)
        self.missed_cleavages = self._int(items[header_map['# Missed Cleavages']])  # convert to integer
        self.charge = self._int(items[header_map['Charge']])
        self.delta_score = self._float(items[header_map['DeltaScore']])
##        # not sure if this is needed
##        if self.delta_score == 1.0:
##            self.delta_score = 0.0
        self.delta_cn = self._float(items[header_map['DeltaCn']])
        self.rank = self._int(items[header_map['Rank']])
        self.se_rank = self._int(items[header_map['Search Engine Rank']])
        self.m_over_z = self._float(items[header_map['m/z [Da]']])
        self.mhplus = self._float(items[header_map['MH+ [Da]']])
        self.theo_mhplus = self._float(items[header_map['Theo. MH+ [Da]']])
        self.deltamass_ppm = self._float(items[header_map['DeltaM [ppm]']])
        self.deltamass_da = self._float(items[header_map['Deltam/z [Da]']])
        self.activation_type = items[header_map['Activation Type']]
        self.ms_order = items[header_map['MS Order']]
        self.isolation_interference = self._float(items[header_map['Isolation Interference [%]']])
        self.sps_mass_matches = self._float(items[header_map['SPS Mass Matches [%]']])
        self.ave_reporter_sn = self._float(items[header_map['Average Reporter S/N']])
        self.ion_inject_time = self._float(items[header_map['Ion Inject Time [ms]']])
        self.rt_min = self._float(items[header_map['RT [min]']])
        self.first_scan = items[header_map['First Scan']]
        self.file_id = items[header_map['File ID']]
        # get the array of reporter ion intensities
        self.plex = len([x for x in list(header_map.keys()) if x.startswith('Abundance:')])
        self.intensities = [0.0 for x in self.channels[:self.plex]]
        for i, channel in enumerate(self.channels[:self.plex]):
            try:
                key = 'Abundance: ' + channel
                self.intensities[i]  =  self._float(items[header_map[key]], missing)
            except KeyError:
                key = 'Abundance: ' + channel[:-1]
                self.intensities[i]  =  self._float(items[header_map[key]], missing)
        self.total = sum(self.intensities)
        self.quant_info = items[header_map['Quan Info']]
        self.xcorr = self._float(items[header_map['XCorr']])
        self.protein_groups = self._int(items[header_map['# Protein Groups']])
        self.q_value = self._float(items[header_map['q-Value']])
        self.pep = self._float(items[header_map['PEP']])
        self.svm_score = self._float(items[header_map['SVM Score']])
        try:
            index = header_map['ptmRS: Best Site Probabilities']
        except KeyError:
            index = header_map['PhosphoRS: Best Site Probabilities']
        try:
            self.best_site_probabilities = items[index]
        except:
            self.best_site_probabilities = ''
                
        # some computed attributes and attributes that are set later
##        self.grouper = self.sequence.upper() + '_' + str(int(round(self.theo_mhplus, 0)))  # key for grouping
        self.grouper = self.sequence.upper() + '_' + str(self.number_phospho_groups)  # key for grouping
        self.meets_all_criteria = False
        self.psm_number = None
        self.grouped_psms = ''
        self.psm_count = None
        self.localization = 'NA'
        return

    def _float(self, string, default=0.0):
        """Converts string to a float, set to "default" if ValueError (missing)
        """
        try:
            val = float(string)
        except ValueError:
            val = default
        return val

    def _int(self, string):
        """Converts string to an integer, zero if ValueError
        """
        try:
            val = int(string)
        except ValueError:
            val = 0
        return val

    def _snoop(self):
        """Diagnostic dump to console
        """
        print(f'Accessions: {self.protein_accessions}')
        print(f'Annotated Sequence: {self.annotated_sequence}')
        print(f'Sequence: {self.sequence}')
        print(f'Ambiguity: {self.psm_ambiguity}')
        print(f'Modifications: {self.modifications}')
        print(f'Number of phospho groups: {self.number_phospho_groups}')
        print(f'Plex: {self.plex}')
        for i, channel in enumerate(self.channels[:self.plex]):
            print(f'{channel}: {self.intensities[i]:,.0f}')
        print(f'Total Intensity: {self.total:,.0f}')
        print(f'Quan Info: {self.quant_info}')
        print(f'q-Value: {self.q_value}')
        print(f'XCorr: {self.xcorr}')
        print(f'DeltaCn: {self.delta_score}')
        print(f'# Missed Cleavages: {self.missed_cleavages}')
        print(f'First Scan: {self.first_scan}')
        print(f'Charge: {self.charge}')
        print(f'RT [min]: {self.rt_min}')
        print(f'MH+ [Da]: {self.mhplus}')
        print(f'Delta Mass [Da]: {self.deltamass_da}')
        print(f'Delta Mass [PPM]: {self.deltamass_ppm}')
        print(f'MeetsAllCriteria: {self.meets_all_criteria}')
        print(f'Group Key: {self.grouper}')
        return

    def make_header(self):
        """Makes a header line for PSM data
        """
        header_line = (['Counter', 'Checked', 'Confidence', 'Identifying Node', 'PSM Ambiguity', 'Annotated Sequence',
                        'Modifications', '# Proteins', 'Master Protein Accessions', 'Protein Accessions',
                        '# Missed Cleavages', 'Charge', 'DeltaScore', 'DeltaCn', 'Rank', 'Search Engine Rank',
                        'm/z [Da]', 'MH+ [Da]', 'Theo. MH+ [Da]', 'DeltaM [ppm]', 'Deltam/z [Da]',
                        'Activation Type', 'MS Order', 'Isolation Interference [%]', 'SPS Mass Matches [%]',
                        'Average Reporter S/N', 'Ion Inject Time [ms]', 'RT [min]', 'First Scan', 'File ID'] +
                        self.channels[:self.plex] +
                        ['Intensity Total', 'Quan Info', 'XCorr', '# Protein Groups', 'q-Value', 'PEP',
                         'SVM Score', 'ptmRS: Best Site Probabilities', 'PSM Number', 'Meets All Criteria',
                         'New Sequence', 'New Modifications', 'Number Phospho Groups', 'Peptide Length',
                         'New Site Prob. Peptide', 'New Site Probabilities', 'Maximum Prob.', 'Minimum Prob.',
                         'Sites', 'Localization', 'Grouping Key', 'Grouped PSM Numbers', 'Number Grouped PSMs',
                         'Match String'])
        return '\t'.join(header_line)

    def make_data(self):
        """Makes a data line for PSM data
        """
        data_list = ([1, self.checked, self.confidence, self.identifying_node, self.psm_ambiguity,
                      self.annotated_sequence, self.modifications, self.number_proteins,
                      self.master_protein_accessions, self.protein_accessions, self.missed_cleavages,
                      self.charge, self.delta_score, self.delta_cn, self.rank, self.se_rank, self.m_over_z,
                      self.mhplus, self.theo_mhplus, self.deltamass_ppm, self.deltamass_da,
                      self.activation_type, self.ms_order, self.isolation_interference, self.sps_mass_matches,
                      self.ave_reporter_sn, self.ion_inject_time, self.rt_min, self.first_scan, self.file_id] +                    
                     self.intensities +                     
                     [self.total, self.quant_info, self.xcorr, self.protein_groups, self.q_value, self.pep,
                      self.svm_score, self.best_site_probabilities, self.psm_number, self.meets_all_criteria,
                      self.new_sequence, self.new_modifications, self.number_phospho_groups, self.peptide_length,
                      self.new_site_prob_peptide, self.new_site_prob, self.max_prob, self.min_prob, self.sites,
                      self.localization, self.grouper, self.grouped_psms, self.psm_count, self.match])
        return '\t'.join([str(x) for x in data_list])

    def make_header_brief(self):
        """Makes a brief header line.
        """
        header_line = ['Counter', 'Group Key', 'Protein Accessions', 'Protein Descriptions',
                       'Quan Info', 'Isolation Inteference [%]', 'SPS Mass Matches',
                       'New Sequence', 'New Modifications', 'Number Phospho Sites',
                       'Peptide Length', 'New Site Prob Peptide', 'New Site Prob Protein',
                       'Site List', 'Localization Status', 'Used PSM Count'] + self.channels[:self.plex]
        return '\t'.join(header_line)

    def make_data_brief(self):
        """Makes a brief data line for PSM data
        """
        data_list = ([1, self.grouper, self.protein_accessions, self.protein_descriptions,
                      self.quant_info, self.isolation_interference, self.sps_mass_matches,
                      self.new_sequence, self.new_modifications, self.number_phospho_groups,
                      self.peptide_length, self.new_site_prob_peptide, self.new_site_prob, self.sites,
                      self.localization, self.psm_count] + self.intensities[:self.plex])
        return '\t'.join([str(x) for x in data_list])
    
    def make_key(self):
        """Makes a column header key
        """
        return 
    # end class

def parse_psm_lines(psm_file, proteins, prot_index, max_qvalue=0.05, max_ppm=20.0, min_intensity=500, missing=0.0, separator='\t'):
    """Parses PSM exports and adds information to PSM objects.
    Returns a list of all PSMs passing q-value, deltamass, and intensity cutoffs.
    """
    # define some fixed ranges of valid PSM attributes
    min_charge = 2
    max_charge = 4
    min_length = 7
    max_length = 40

    # initialize counters and psm list
    total = 0       #0
    top = 0         #1
    unmod = 0       #2
    mod = 0         #3
    valid = 0       #4
    reject = 0      #5
    qval_good = 0   #6
    qval_bad = 0    #7
    in_ppm = 0      #8
    out_ppm = 0     #9
    above_int = 0   #10
    below_int = 0   #11    
    psm_list = []

    # start parsing PSM file
    start = False    
    for line in open(psm_file, 'r'):
        line = line.strip()
        if not line:
            continue    # skip blank lines
        
        if not start:   # skip lines until header
            if line.startswith('Checked\tConfidence'):   # look for header line
                psm_map = make_header_map(line)                
                start = True
                total = 1
                continue
        else:
            # parse table line
            psm = PSM(line, psm_map, proteins, prot_index, missing, separator)
            total += 1

##            # this may filter out low q-value matches (?)
##            if not psm.protein_accessions: 
##                total += -1
##                reject += -1
##                continue

            # skip non-top-ranked matches (Percolator considers top 5 matches)
            if psm.rank == 1:
                top += 1
            else:
                continue

            # add PSM to list
            psm_list.append(psm)

            # test various criteria, count passing/failing PSMs, set "meets_all_criteria" flag (False by default)
            # test q-value cutoff
            if psm.q_value <= max_qvalue:
                qval_good += 1
            else:
                qval_bad += 1
                reject += 1
                continue
                
            if abs(psm.deltamass_ppm) <= max_ppm:
                in_ppm += 1
            else:
                out_ppm += 1
                reject += 1
                continue

            if test_trimmed_average_intensity(psm.intensities, min_intensity):
                above_int += 1
            else:
                below_int += 1
                reject += 1
                continue

            if min_charge <= psm.charge <= max_charge:
                pass
            else:
                reject += 1
                continue

            if min_length <= psm.peptide_length <= max_length:
                pass
            else:
                reject += 1
                continue

            if psm.number_phospho_groups != 0:
                mod += 1                
            else:
                unmod += 1
                reject += 1
                continue

            psm.meets_all_criteria = True
            valid += 1
            
    counters = [total, top, unmod, mod, valid, reject, qval_good, qval_bad, in_ppm, out_ppm, above_int, below_int]                
    return counters, psm_list

def parse_mods(modstring):
    """Parses PD modification descriptions to get
    modification types and count of affected residues.
    
    modifications returned as dictionary of modtypes and dictionary of residues and counts.
    empty mod descriptions should return empty structures.
    """
    mods = {}

    # split modification description string
    modlist = modstring.split(';')
    for mod in modlist:
        mod = mod.strip()   # get rid of whitespace
        temp = mod[:-1].split('(') # separate part inside ()
        residue = ''.join([c for c in temp[0] if not c.isdigit()]) # ignore positions
        modtype = temp[1]

        # mods is a dictionary of modtype where each modtype value is a count dictionary {residue: count}
        if modtype in mods:
            if residue in mods[modtype]:
                mods[modtype][residue] += 1
            else:
                mods[modtype][residue] = 1
        else:
            mods[modtype] = {residue: 1}

    return mods

def amino_acid_frequency(pepstring, aa_freq):
    """Counts amino acids frequencies of peptide sequences from PD,
    sequence and frequency dictionary are passed as arguments.
    """
    pepstring = pepstring.upper()
    for aa in pepstring:
        if aa in aa_freq:
            aa_freq[aa] +=1
        else:
            aa_freq[aa] = 1
    return

def update_dictionary(big, little):
    """Big and little are dictionaries with dictionaries of count values.
    Big probably has more keys than little.
    """
    for key in little:
        if key in big:
            for secondkey in little[key]:
                if secondkey in big[key]:
                    big[key][secondkey] += little[key][secondkey]
                else:
                    big[key][secondkey] = little[key][secondkey]
        else:
            big[key] = little[key]
    return

def fixed_or_variable(all_mods, aa_freq):
    """Assigns special symbols to variable mods. None assigned to fixed mods.
    """
    symbols = ['*', '#', '@', '^', '~', '$', '[', ']']

    mod_type = {}
    global_mod_freq = {}
    for mod in all_mods:
        for residue in all_mods[mod]:
            if all_mods[mod][residue] == aa_freq[residue]:
                mod_type[mod] = None
            else:
                mod_type[mod] = True
    for mod in all_mods:
        if mod_type[mod]:
            global_mod_freq[mod] = sum(all_mods[mod].values())
    variable_freq = sorted(global_mod_freq.items(), key=lambda x: x[1], reverse=True)
    for i, (mod, count) in enumerate(variable_freq):
        mod_type[mod] = symbols[i]
    return mod_type

def get_variable_positions(psm, mod_type):
    """Parses PD modification descriptions to get modification positions and symbols
    """
    modmask = {}

    # split modification description string
    modlist = psm.modifications.split(';')
    for mod in modlist:
        mod = mod.strip()   # get rid of whitespace
        temp = mod[:-1].split('(') # separate part inside ()
        position = ''.join([c for c in temp[0] if c.isdigit()]) # get positions
        try:
            position = int(position)
        except ValueError:
            if temp[0] == 'N-Term':
                position = 0
            elif temp[0] == 'C-Term':
                position = len(psm.sequence)
            else:
                position = -1
        modtype = temp[1]

        if mod_type[modtype]:   # just variable mods
            modmask[position-1] = mod_type[modtype]

    return modmask

def fix_PTM_info(psm, mod_type):
    """Makes SEQUEST-style sequence stings and removes static mods from modifications strings
    """
    new_seq = list(psm.sequence.upper())
    new_symbols = ['' for x in new_seq]
    modmask = get_variable_positions(psm, mod_type)
    for index in modmask:
        new_symbols[index] = modmask[index]
    psm.new_sequence = ''.join([j for i in zip(new_seq, new_symbols) for j in i])

    new_mod_list = [x.strip() for x in psm.modifications.split(';') if mod_type[x[:-1].split('(')[1]]]
    psm.new_modifications = '; '.join(new_mod_list)

def sequence_length_analysis(psm_list):
    length_count = [0 for x in range(60)]
    print()
    for i in range(60):
        length_count[i] = len([x for x in psm_list if x.peptide_length == i])

    total = float(sum(length_count))
    running_tot = 0.0
    for i in range(60):
        running_tot += length_count[i]
        print(i, length_count[i], round(100.0*running_tot/total, 2))
    print()
    return

class RSProb:
    """Container for parsed PhosphoRS data
    """
    def __init__(self, one_site):
        if '(Phospho)' not in one_site:
            print('...WARNING: non-phospho site:', one_site)
        try:
            residue_position = one_site.split('(Phospho)')[0]
            residue = residue_position[0]
            position = int(residue_position[1:])
            self.residue = [residue]   # amino acid (S, T, or Y)
            self.position = [position]   # relative position in peptide string
            self.probability = float(one_site.split(':')[1])  # assigned probability in percent
        except:
            print('...WARNING: site parsing failed:', one_site)

def parse_RSProb(psm):
    """Parse phosphoRS probability string and returns N sites of the highest probability.
    """
    RS_full_list = []
    if ((not psm.best_site_probabilities) or
        (psm.best_site_probabilities == 'Too many NL-allowing PTMs') or
        (psm.best_site_probabilities == 'Too many isoforms') or
        (psm.best_site_probabilities == 'Inconclusive data')):
        psm.new_site_prob = ''
        psm.new_site_prob_peptide = ''
        psm.max_prob = 0.0
        psm.min_prob = 0.0
        psm.sites = ''
        return
    for prob in psm.best_site_probabilities.split(';'):
        prob = prob.strip()
        prob = RSProb(prob)     # this can return empty object for odd PhosphoRS site descriptions
        prob.position[0] += (psm.start-1)
        RS_full_list.append(prob)        

    # sort descending and get the top N sites
    RS_full_list.sort(reverse=True, key=lambda x: x.probability)
    RS_list = RS_full_list[:psm.number_phospho_groups]
    psm.max_prob = RS_list[0].probability
    psm.min_prob = RS_list[-1].probability

    # check for prob ties
    for rs in RS_full_list[psm.number_phospho_groups:]:
        if rs.probability == RS_list[-1].probability:
            RS_list[-1].residue += rs.residue
            RS_list[-1].position += rs.position
    
    # put together output string           
    site_list_prob = []
    for rs in RS_list:
        residue_position = ', '.join(rs.residue) + ' (' + ', '.join([str(x) for x in rs.position]) + ')'
        # maybe combine probabliities for ties?
        site_list_prob.append('%s: %0.1f' % (residue_position, len(rs.position)*rs.probability))
    psm.new_site_prob = '; '.join(site_list_prob)
#    psm.sites = '; '.join(site_list)
    
    # put together output string with relative positions          
    site_list_prob = []
    for rs in RS_list:
        positions = [x-psm.start+1 for x in rs.position]
        residue_position = ', '.join(rs.residue) + ' (' + ', '.join([str(x) for x in positions]) + ')'
        site_list_prob.append('%s: %0.1f' % (residue_position, rs.probability))
    psm.new_site_prob_peptide = '; '.join(site_list_prob)

    # make site list (no ties), sorted by residue position
    site_set = set()
    for rs in RS_list:
        [site_set.add(x) for x in zip(rs.residue, rs.position)]
    site_list = sorted(list(site_set), key=lambda x: x[1])
    site_list = [str(x[0])+str(x[1]) for x in site_list]
    psm.sites = '; '.join(site_list)

    return

def make_protein_index(proteins):
    """Indexes proteins
    """
    prot_index = {}
    skip = set(['sp', 'tr', 'gi', 'ref', ''])
    for i, p in enumerate(proteins):
        accs = p.accession.split('|')
        for acc in accs:
            if acc in skip:
                continue
            prot_index[acc] = i
    return prot_index

def lookup_headers(acc_string, proteins, prot_index):
    """Fetches header lines for proteins.
    """
    header_list = []
    for acc in [x.strip() for x in acc_string.split(';')]:
        prot = proteins[prot_index[acc]]
        header_list.append(prot.description)
    return '; '.join(header_list)


def lookup_peptides(psm_list, proteins, prot_index):
    """Finds starting residue number for peptide sequence in protein sequence.
    """
    for psm in psm_list:
        try:
            # eventually add lookup of all accessions, just first to test
            acc = psm.protein_accessions.split(';')[0].strip()
            prot = proteins[prot_index[acc]]
            psm.match = prot.findPeptide(psm.new_sequence, pad_count=3)
            psm.start = psm.match[0][1]
        except IndexError:
            print()
            print('...peptide lookup issue:')
            print('pre-acc:', psm.accessions)
            print('acc:', acc)
            print('index:', prot_index[acc])
            print('full acc:', prot.accession)
            print('peptide:', psm.new_sequence)
            print('base_peptide:', prot.base_peptide_sequence(psm.new_sequence))
            print('peptide in sequence?', psm.new_sequence in prot.sequence)
            print()
            psm.start = 0
        
def make_group_index(psm_list):
    """Makes a dictionary of PSM index lists (values) that have the same grouper string (key)
    """
    group_index = {}
    for i, psm in enumerate(psm_list):
        if psm.meets_all_criteria:
            if psm.grouper in group_index:
                group_index[psm.grouper].append(i)
            else:
                group_index[psm.grouper] = [i]
    return group_index

def analyze_modifications(psm_list):
    """gets freqeuncies of amino acids in sequences, gets frequencies of modifications
    """
    aa_freq = {}
    all_mods = {}
    
    aa_freq['N-Term'] = aa_freq['C-Term'] = len(psm_list)
    for psm in psm_list:
        amino_acid_frequency(psm.sequence, aa_freq)
        mods = parse_mods(psm.modifications)
        update_dictionary(all_mods, mods)

    return aa_freq, all_mods

def print_modification_report(all_mods, mod_type, log_obj):
    """Prints a summary of variable and static modifications
    """
    variable = []
    for obj in [None, log_obj]:
        print('\nVariable modifications:', file=obj)   # print mod type, symbol, affected residues
    for mod in mod_type:
        if mod_type[mod]:
            variable.append([mod, mod_type[mod], sorted(all_mods[mod].keys()), sum(all_mods[mod].values())])
    variable = sorted(variable, reverse=True, key=lambda x: x[3])   # order mods by decreasing frequency
    for mod in variable:
        for obj in [None, log_obj]:
            print('  ', mod[0], mod[1], mod[2], file=obj)

    static = []
    for obj in [None, log_obj]:
        print('Static modifications:', file=obj)   # print mod type and afected residues
    for mod in mod_type:
        if not mod_type[mod]:
            static.append([mod, sorted(all_mods[mod].keys()), sum(all_mods[mod].values())])
    static = sorted(static, reverse=True, key=lambda x: x[2])   # order mods by decreasing freqiuency
    for mod in static:
        for obj in [None, log_obj]:
            print('  ', mod[0], mod[1], file=obj)
    for obj in [None, log_obj]:
        print(file=obj)
    return

def test_phosphoRS_localization(psm_list):
    """Need some localization probability cutoff testing here.
    Should depend on number of sites in peptide.
    """
    for psm in psm_list:
        if psm.number_phospho_groups == 0:
            if SKIP_UNMOD:
                psm.meets_all_criteria = False  # this rejects unmodified peptides
        elif psm.number_phospho_groups == 1 and psm.min_prob < 45.0:
            psm.meets_all_criteria = False
        elif psm.number_phospho_groups == 2 and psm.min_prob < 0.30:
            psm.meets_all_criteria = False
        elif psm.number_phospho_groups == 3 and psm.min_prob < 22.5:
            psm.meets_all_criteria = False

def make_grouped_PSMs(psm_list, group_index, missing=150.0):
    """Groups PSMs by grouper key, keeps psm with best q-value, sums intensities
    """
    new_psm_list = []
    for key in group_index:
        if len(group_index[key]) == 1:
            # nothing to group
            best_psm = copy.deepcopy(psm_list[group_index[key][0]])
            best_psm.localization = 'consistent'
            best_psm.grouped_psms = str(best_psm.psm_number)
            best_psm.psm_count = 1            
            # replace any zeros
            best_psm.intensities = [(x if x > 0.0 else missing) for x in best_psm.intensities]
            new_psm_list.append(best_psm)
        else:
            # pick group member with best q-value to represent group
            psm_group_list = sorted([psm_list[i] for i in group_index[key]], key=lambda x: x.q_value)
            best_psm = copy.deepcopy(psm_group_list[0])
            best_psm.psm_count = len(psm_group_list)
            best_psm.grouped_psms = '; '.join([str(x.psm_number) for x in psm_group_list])
            # sum the intensities for the group
            summed = [0.0 for i in best_psm.channels]
            for i in range(len(best_psm.channels[:best_psm.plex])):
                for j in range(len(psm_group_list)):
                    summed[i] += psm_group_list[j].intensities[i]
            best_psm.intensities = copy.deepcopy(summed)
            # do zero replacement
            best_psm.intensities = [(x if x > 0.0 else missing) for x in best_psm.intensities]
            # sum reporter ions
            best_psm.total = sum(best_psm.intensities)
            # see if site localization is the same for all PSMs
            site_set = set([x.sites for x in psm_group_list])
            if len(site_set) == 1:
                best_psm.localization = 'consistent'
            else:
                best_psm.localization = 'varies'
            new_psm_list.append(best_psm)
    return new_psm_list
            
    

###########################################
######## main program starts here #########
###########################################

print('====================================================')
print(' program "PD_TMT_phospho_processer.py", version 3.1 ')
print('====================================================')
print('Ran on:', time.ctime())
print('INTENSITY = %s, QVALUE = %s, MISSING = %s, PPM = %s' % (INTENSITY, QVALUE, MISSING, PPM))
    
# get the PSM results PD export file information
default_location = r'F:\PSR_Core_Analysis'
if not os.path.exists(default_location):
    default_location = os.getcwd()
print('\nSelect the PSM export file (tab-delimited text)')
psm_filename = PAW_lib.get_file(default_location,
                                [('Text files', '*.txt'), ('All files', '*.*')],
                                'Select a PD 3.x phospho PSM export file')
if not psm_filename: sys.exit()     # exit if no file selected
print('...', psm_filename, 'was selected')

# get PSM path and basename, open log file
path, base = os.path.split(psm_filename)
base_log = base.replace('.txt', '_log.txt')
log_obj = open(os.path.join(path, base_log), 'a')

print('\n====================================================', file=log_obj)
print(' program "PD_TMT_phospho_processer.py", version 3.1 ', file=log_obj)
print('====================================================', file=log_obj)
print('Ran on:', time.ctime(), file=log_obj)
print('INTENSITY = %s, QVALUE = %s, MISSING = %s, PPM = %s' % (INTENSITY, QVALUE, MISSING, PPM), file=log_obj)

# get the FASTA database file
print('Select the FASTA protein database file')
db_filename = PAW_lib.get_file(default_location,
                               [('Fasta files', '*.fasta'), ('All files', '*.*')],
                               'Select a FASTA database')
if not db_filename: sys.exit()
print('...', db_filename, 'was selected')

# read in the protein sequences
print('Reading proteins...')
f = PAW_lib.FastaReader(db_filename)
p = PAW_lib.Protein()
proteins = []
while f.readNextProtein(p, False):
    prot = copy.deepcopy(p)
    proteins.append(prot)
prot_index = make_protein_index(proteins)
print('Number of proteins in FASTA file:', len(proteins))

# get the psm information from the PSM export file
counts, psm_list = parse_psm_lines(psm_filename, proteins, prot_index, QVALUE, PPM, INTENSITY, 0.0, SEPARATOR)

# analyze PTMs: fixed or variable? what residues? what names?
aa_freq, all_mods = analyze_modifications(psm_list)

# analyze mods and print summary
mod_type = fixed_or_variable(all_mods, aa_freq)
print_modification_report(all_mods, mod_type, log_obj)    

# add alternatively formatted sequence string (SEQUEST style), reformat modifications string (remove static mods)
for psm in psm_list:
    fix_PTM_info(psm, mod_type)

# find peptide starting residue numbers (do before RS probabilities)
lookup_peptides(psm_list, proteins, prot_index)

# add a simplified RS probability string (top N sites, check for ties)
for psm in psm_list:
    parse_RSProb(psm)

# test for some minimal phosphoRS localization and update "meets_all_criteria"
# need to set "meets_all_criteria" since it controls combining of similar PSMs
test_phosphoRS_localization(psm_list)

# add PSM number for indexing purposes (DO NOT sort psms until after grouping is done)
for i, psm in enumerate(psm_list):
    psm.psm_number = i

# figure out which PSMs to group together (sum intensities)
group_index = make_group_index(psm_list)    # this is a filtered list ("meets_all_criteria")
new_psm_list = make_grouped_PSMs(psm_list, group_index, MISSING)    # do the zero replacement here

######################################################
# when everything is done, sort and write to new files

# get original file name/path and ext (for making otput filenames)
path, extension = os.path.splitext(psm_filename)

# open full results file for all filtered PSMs
new_path = path + '_psm_filtered_all'
psmout = open(new_path + extension, 'w')

# print the header line
print(psm_list[0].make_header(), file=psmout)

# print PSM data to file sorted by decreasing total intensity
psm_list = sorted(psm_list, key=lambda x: x.total, reverse=True)
for psm in psm_list:
    print(psm.make_data(), file=psmout)

# print keys at end of table (psm should still be last psm in psm_list)
print(psm.make_key(), file=psmout)

# print the parsing statistics
for out in [None, psmout, log_obj]:
    print('The PSM export file was:', base, file=out)
    print('The FASTA database was:', db_filename, file=out)
    print('\nThere were', counts[0], 'Total rows in PSM table export', file=out)    
    print('There were %i top-ranked PSMs (%i non-top-ranked)' % (counts[1], counts[0]-counts[1]), file=out)    
    print('There were %i PSMs passing q-value cutoff (%i failed)' % (counts[6], counts[7]), file=out)
    print('There were %i PSMs passing deltamass ppm cutoff (%i failed)' %(counts[8], counts[9]), file=out)
    print('There were %i PSMs with reporter ions above intensity cutoff (%i below)' % (counts[10], counts[11]), file=out)     
    print('   (%0.2f%% of reporter ion scans were too low intensity)' % (100.0 * counts[11] / (counts[10]+counts[11]),), file=out)
    print('There were %i phosphorylated PSMs (%i unmodified)' % (counts[3], counts[2]), file=out)
    print('\nThere were %i PSMs that passed all criteria (%i failed tests)' % (counts[4], counts[5]), file=out)
    print('   (%0.0f PSMs are estimated to be incorrect matches)' % (QVALUE * counts[4],), file=out)
    
# close files
psmout.close()
log_obj.close()

# open results file for combined PSMs, print header lines
new_path = path + '_psm_filtered_combined'
psmout = open(new_path + extension, 'w')
new_path = path + '_psm_filtered_combined_brief'
psmout_brief = open(new_path + extension, 'w')
print(new_psm_list[0].make_header(), file=psmout)
print(new_psm_list[0].make_header_brief(), file=psmout_brief)

# print PSM data to file sorted by decreasing total intensity
new_psm_list = sorted(new_psm_list, key=lambda x: x.total, reverse=True)
for psm in new_psm_list:
    print(psm.make_data(), file=psmout)
    print(psm.make_data_brief(), file=psmout_brief)

# print keys at end of tables
print('\nThe FASTA database was:', db_filename, file=psmout)
print(psm.make_key(), file=psmout)
    
# close files
psmout.close()
psmout_brief.close()

# end
