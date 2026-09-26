# TN93 Cluster script

# Imports -------------------------------------------------------------
import os
import sys
import argparse
import json
import shutil
import random

# Declares
# Argparse here
arguments = argparse.ArgumentParser(description='Cluster an MSA with genetic distance (TN93)')

arguments.add_argument('-i', '--input',            help = 'MSA file to process',                                  required = True, type = str )
arguments.add_argument('-o', '--output_fasta',           help = 'Output json file',                                     required = True, type = str)
arguments.add_argument('-j', '--output_json',           help = 'Output json file',                                     required = True, type = str)
arguments.add_argument('--threshold',              help = 'Distance threshold for clustering query sequences',    required = True, type = float)
#arguments.add_argument('--step',                   help = 'Distance threshold for clustering query sequences',    required = True, type = float)
arguments.add_argument('-m', '--max_retain',    help = 'The maximum number of sequences to retain',               required = True, type = int)
arguments.add_argument('-r', '--reference_seq',    help = 'The maximum number of sequences to retain',               required = False, type = str)
arguments.add_argument('--preserve_list',       help = 'File containing sequence names that must be retained',     required = False, type = str)
arguments.add_argument('--preserve',            help = 'Comma-separated sequence names that must be retained',     required = False, type = str)

settings = arguments.parse_args()

input_stamp = os.path.getmtime(settings.input)

# Output is actuall the compressed.fas, thats what we want.
# Input is the MSA.
# also pass in {GENE}.{query/reference}.json filename
_ref_seq_name = ""

if settings.reference_seq:
    with open (settings.reference_seq) as fh:
        for l in fh:
            if l[0] == '>':
                _ref_seq_name = l[1:].split (' ')[0].strip()
                break
            #end if
        #end for
    #end with
#end if
            
#print ("Reference seq_name %s" % _ref_seq_name)

# files
cluster_json = settings.output_json          #output file # This is actually a json file
compressed_fasta = settings.output_fasta           # .compressed.fas
msa_strike_ambigs = settings.input

# values
threshold = settings.threshold
#step = settings.step
max_toRetain = settings.max_retain
preserve_names = set()

if settings.preserve_list:
    with open(settings.preserve_list) as fh:
        preserve_names.update(line.strip() for line in fh if line.strip())

if settings.preserve:
    preserve_names.update(x.strip() for x in settings.preserve.split(",") if x.strip())

# Need to load this from a config.json file.
task_runners = {}
task_runners['tn93-cluster'] = "tn93-cluster"

# Helper functions -----------------------------------------------------

def run_command (exec, arguments, filename, tag):
    global input_stamp
    cmd = " ".join ([exec] + arguments)
    print(cmd)
    result = os.system (cmd)
    if result != 0:
        print ('Command exection failed code %s' % result)
        sys.exit(result)
        return input_stamp
    return os.path.getmtime(filename)
#end method

def read_fasta_records(in_file):
    records = {}
    current_name = None
    current_seq = []

    with open(in_file, "r") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if current_name is not None:
                    records[current_name] = "".join(current_seq).replace(" ", "")
                current_name = line[1:].strip().split(" ")[0]
                current_seq = []
            else:
                current_seq.append(line.strip())

    if current_name is not None:
        records[current_name] = "".join(current_seq).replace(" ", "")

    return records
#end method

def write_fasta_record(out_handle, seq_name, seq, check_uniq):
    seq_id = ">" + seq_name
    while seq_id in check_uniq:
        seq_id = ">" + seq_name + '_' + ''.join(random.choices ('0123456789abcdef', k = 10))
    check_uniq.add(seq_id)
    print(seq_id + "\n" + seq.replace(" ", "") + "\n", file = out_handle)
#end method

def cluster_to_fasta (in_file, out_file, source_fasta, ref_seq = None, preserve = None):
    preserve = preserve or set()
    input_records = read_fasta_records(source_fasta)
    # assert that in_file exists, otherwise this will crash if tn93-cluster did not run.
    with open (in_file, "r") as fh:
        cluster_data = json.load (fh)
        #print (colored('Running ... converting representative clusters to .FASTA\n', 'cyan'))
        check_uniq = set ()
        retained_count = 0
        print("# Saving to fasta:", out_file)
        with open (out_file, "w") as fh2:
            for c in cluster_data:
                cc = c['centroid'].split ('\n')
                members = c.get('members', [])
                emitted_preserved = False

                for member in members:
                    member_name = member[1:].strip().split(" ")[0] if member.startswith(">") else member.strip().split(" ")[0]
                    if member_name in preserve:
                        seq = input_records.get(member_name)
                        if seq is None:
                            print("# WARNING: preserved sequence not found in current FASTA:", member_name, file=sys.stderr)
                            continue
                        write_fasta_record(fh2, member_name, seq, check_uniq)
                        retained_count += 1
                        emitted_preserved = True
                    #end if
                #end for

                if emitted_preserved:
                    continue
                #end if

                if ref_seq:
                    if ref_seq in members:
                        seq = input_records.get(ref_seq, cc[1])
                        write_fasta_record(fh2, ref_seq, seq, check_uniq)
                        retained_count += 1
                        continue
                    #end if
                #end if

                centroid_name = cc[0][1:].strip().split(" ")[0] if cc[0].startswith(">") else cc[0].strip().split(" ")[0]
                centroid_seq = cc[1] if len(cc) > 1 else input_records.get(centroid_name, "")
                write_fasta_record(fh2, centroid_name, centroid_seq, check_uniq)
                retained_count += 1
            #end for
         #end with
    #end with
    return (os.path.getmtime(out_file), len(cluster_data), retained_count)
#end method

# Main subroutine -----------------------------------------------------

input_file = msa_strike_ambigs

while True:
    input_stamp = run_command (task_runners['tn93-cluster'], ['-f', '-o', cluster_json, '-t', "%g" % threshold, input_file], cluster_json, "extract representative clusters at threshold %g" % threshold)    
    if _ref_seq_name != "":
        input_stamp, cluster_count, retained_count = cluster_to_fasta (cluster_json, compressed_fasta, input_file, _ref_seq_name, preserve_names)
    else:
        input_stamp, cluster_count, retained_count = cluster_to_fasta (cluster_json, compressed_fasta, input_file, preserve = preserve_names) # changes the json to fasta, also returns the count len().
    #end if
    
    print("# Current number of sequences", cluster_count)
    print("# Current retained FASTA records", retained_count)
    if preserve_names:
        print("# Preserving %d requested sequences during TN93 clustering" % len(preserve_names))
    if cluster_count <= max_toRetain:
        #shutil.copy (msa_SA, msa)
        break
    else:
       step = threshold * 0.25 
       threshold += step
       input_file = compressed_fasta
    #end if
#end while

sys.exit(0)
# End of file
