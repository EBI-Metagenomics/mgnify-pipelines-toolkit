
from collections import defaultdict
import os
import subprocess

from mgnify_pipelines_toolkit.constants.regex_ambiguous_bases import _AMBIGUOUS_BASES_DICT, _AMBIGUOUS_BASES_DICT_REV

def split_dir_into_sample_paths(_DIR):

    file_list = os.listdir(_DIR)
    file_list = [ file for file in file_list if '.fastq' in file and ('_1' in file or '_2' in file) ]
    sample_set = set()
    [ sample_set.add(f"{_DIR}/{file.split('_')[0]}") for file in file_list ]
    sample_list = sorted(list(sample_set))

    return sample_list

def get_read_count(read_path, type='fastq'):

    cmd = []
    
    if type == 'fastq':
        cmd = [
            'zgrep',
            '-c',
            '^@',
            read_path
        ]

    elif type == 'fasta':
        cmd = [
            'grep',
            '-c',
            '^>',
            read_path
        ]

    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()

    read_count = int(stdout.strip())

    return read_count

def build_cons_seq(cons_list, read_count, cons_threshold=0.80, do_not_include=[], counter=1):

    cons_seq = ''
    cons_confs = []

    for count_dict in cons_list:
        max_count = 0
        cons_dict = defaultdict(float)

        if counter in do_not_include:
            counter += 1
            cons_seq += 'N'
            continue 
        
        for base, count in count_dict.items():
            if base not in ('A', 'T', 'C', 'G'):
                continue

            cons_dict[base] = count/read_count
            
            if count > max_count:
                max_count = count

        counter += 1
        
        try:
            max_prop = max_count/read_count

            cons_bases = []
            curr_prop = 0.0
            sorted_cons_dict = dict(sorted(cons_dict.items(), key=lambda x:x[1], reverse=True))

            for base, prop in sorted_cons_dict.items():
                cons_bases.append(base)
                curr_prop += prop
                if curr_prop >= cons_threshold:
                    break

            cons_bases = sorted(cons_bases)

            if len(cons_bases) == 1:
                cons_seq += cons_bases[0]
            else:
                amb_string = ','.join(cons_bases)
                amb_base = _AMBIGUOUS_BASES_DICT_REV[amb_string]
                cons_seq += amb_base
                
        except ZeroDivisionError:
            max_prop = 0.0

        cons_confs.append(max_prop)


    return cons_seq, cons_confs

def primer_regex_query_builder(primer):

    query = ''

    for char in primer:
        if char in ('A', 'C', 'T', 'G'):
            query += char
        else:
            query += str(_AMBIGUOUS_BASES_DICT[char])

    return query

def build_mcp_cons_dict_list(mcp_count_dict, mcp_len):
    """
    Generate list of dictionaries of base conservation for mcp output (mcp_cons_list)
    e.g. [{'A':0.9, 'C':0.1}, {'T':1.0}, ....] for every base position
    """

    mcp_cons_list = []

    for i in range(mcp_len):
        index_base_dict = defaultdict(int)
        for mcp in mcp_count_dict.keys():
            if len(mcp) < mcp_len:
                continue
            base = mcp[i]
            index_base_dict[base] += mcp_count_dict[mcp]
        mcp_cons_list.append(index_base_dict)
    
    return mcp_cons_list

def fetch_mcp(fastq, prefix_len, start=1, rev=False):
    """
    Runs a the "find_most_common_prefixes.sh" script to generate the mcps from a fastq file.

    Outputs dictionary containing counts for each generated MCP in the fastq.

    """

    start = str(start)
    prefix_len = str(prefix_len)

    # Check strand requested
    if not rev:
        rev = '0'
    else:
        rev = '1'

    cmd = [
        'bash',
        './mgnify_pipelines_toolkit/analysis/amplicon/find_most_common_prefixes.sh',
        '-i',
        fastq,
        '-l',
        prefix_len,
        '-c',
        '10',
        '-b',
        start,
        '-r',
        rev
    ]

    # Run command
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    
    # Capture stdout
    output = stdout.decode('ascii')
    output = output.strip()
    output = output.split('\n')

    mcp_count_dict = defaultdict(int)

    # Process output into mcp count dictionary
    for line in output:
        line = line.strip()
        temp_lst = line.split(' ')
        if len(temp_lst) == 2:
            mcp_count_dict[temp_lst[1]] = int(temp_lst[0])

    return mcp_count_dict