#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  2 13:04:15 2020
@author: Alexander G. Lucaci
The idea for this script is that:
    Given a protein sequence and a transcript sequence
    I find the codons by stepping over the transcript sequence until the translated sequence matches the protein sequence
    that way, I have only the codons and not the additional sequences from the transcript
    (Which may be useful later)
    I will also create two output files
        One with the STOP codon stripped (this makes it hyphy compatible.)
        One with the STOP codons (may be useful later, codon bias?)
"""

# =============================================================================
# Imports
# =============================================================================
from Bio import SeqIO
import sys
import argparse

# =============================================================================
# Declares
# =============================================================================
PROTEIN_FASTA = snakemake.params.Prot
TRANSCRIPTS_FASTA = snakemake.params.Nuc
OUTPUT = snakemake.params.Out

results = []
no_match = []

logfile = OUTPUT + ".log"

with open(logfile, "w") as fh2:
    print("", file=fh2)
#end with

fh2.close()

# =============================================================================
# Helper functions
# =============================================================================
def already_in_results(transcript_desc):
    global results
    Found = False
    for record in results: # results stores transcript records that passed
        if transcript_desc == record.description: # already exists?
            #print("# Already found this, move on.")
            #start += 1
            #continue
            Found = True
            break
        #end if
    #end for
    return Found
#end method

def log(msg, logfile):
    with open(logfile, "a") as fh2:
        print(msg, file=fh2)
    #end with
    fh2.close()
#end method

def Process(protein_desc, protein_seq, TRANSCRIPTS_FASTA, species, seq_threshold = None): #protein species
    global results, no_match, logfile
    start = 0
    NT_SEQ_LENGTH = len(protein_seq) * 3
    # loop over all of the TRANSCRIPTS_Fasta seq

    #if seq_threshold != None:
    #    seq_threshold = seq_threshold / 2
    #    log("# Sequence length treshold is: " + str(seq_threshold), logfile)
    #end if
  
    with open(TRANSCRIPTS_FASTA, "r") as transcript_handle:
        for m, transcript_record in enumerate(SeqIO.parse(transcript_handle, "fasta")):
            DONE = False
            # Grab Transcript Data
            transcript_id = transcript_record.id
            transcript_desc = transcript_record.description
            transcript_seq = transcript_record.seq

            #if seq_threshold != None:
            #    #seq_threshold = seq_threshold / 2
            #    #log("# Sequence length treshold is: " + seq_threshold, logfile)
            #    if len(transcript_seq) < seq_threshold:
            #        with open(logfile, "a") as fh2:
            #            print("# Sequence length minimum is:", seq_threshold, file=fh2)
            #            print("# Skipping this sequence due to it failing the sequence length minimum:", transcript_record.description, file = fh2)
            #        #end with
            #        fh2.close()
            #        continue
            #    #end if
            ##end if
 
            if species not in transcript_desc: 
                #print("# Mismatch between species") # move on to the next one
                continue # only look at sequences from your species, not something similar.
            #end if

            #print("TX DESC:", transcript_desc)
            # can be a separate subroutine.
            start = 0
            NT_SEQ_LENGTH = len(protein_seq) * 3
            #print(NT_SEQ_LENGTH)
            while start < len(str(transcript_seq)):
                coding_dna = ""
                try:
                    coding_seq = transcript_seq[start: start + NT_SEQ_LENGTH]
                    #print(coding_seq)
                    coding_dna = coding_seq.translate() #translated, universal code
                except:
                    pass
                    #print("#ERROR", coding_dna)
                    #report ERROR
                #end try
                
                exists = already_in_results(transcript_desc)
                if coding_dna == str(protein_seq) and exists == False: # Exit upon first match, may be useful to see how many matches.
                  DONE = True
                  #print(codon_dna)
                  break
                else:
                  start += 1
                #end if
            #end while
            if DONE == True:  break
        #end for
    #end with

    if DONE == True:
        #return transcript_id, transcript_desc, coding_seq
        transcript_record.seq = coding_seq        
        # If it fails, return a known
        return transcript_record
    else:
        return "NO_MATCH"
    #end if
#end method


# Average length of sequences?
# Skip those that are half the average.

def average_seq_length():
    pass

# =============================================================================
# Main subroutine.
# =============================================================================
def main(PROTEIN, TRANSCRIPTS): # Really to verify things.
    print("# TRANSCRIPT INPUT FILE:", TRANSCRIPTS)
    print("# PROTEIN INPUT FILE:", PROTEIN)
    protein_list = []
    transcript_list = []

    with open(TRANSCRIPTS, "r") as handle:
        trans_count = 0 
        for record in SeqIO.parse(handle, "fasta"):
            trans_count +=1
            transcript_list.append(record.description)
        print("# Transcripts:", trans_count)    
    handle.close()
    
    with open(PROTEIN, "r") as handle:
        #x = SeqIO.parse(handle, "fasta")
        #print(len(x))
        prot_count = 0 
        for record in SeqIO.parse(handle, "fasta"):
            prot_count +=1
            protein_list.append(record.description)
        print("# Proteins:", prot_count)
    handle.close()
    
    return trans_count, prot_count
#end method

# =============================================================================
# Main
# =============================================================================
#Verify files exist
print("# =============================================================================")
print("# Processing... ")
print("# =============================================================================")
trans_count, prot_count = main(PROTEIN_FASTA, TRANSCRIPTS_FASTA)
#Looks like species all match up in transcript and protein fasta.
#This is exceptional, will need to look for species name (from protein desc.) in transcript desc.

# DEBUG
#sys.exit(1)
      
# Create empty output file.
with open(OUTPUT, "w") as fh:
    fh.write("")
fh.close()

# Main program
successful_count = 0 
num_errors = 0
errors_IDs = []

# Get average sequence length.

prot_seq_lengths = []

with open(PROTEIN_FASTA, "r") as prot_handle:
    for n, record in enumerate(SeqIO.parse(prot_handle, "fasta")):
        prot_seq_lengths.append(len(record.seq))
    #end for
#end with
prot_handle.close()

avg_sequence_length = sum(prot_seq_lengths) / len(prot_seq_lengths) 
avg_sequence_length_nt = avg_sequence_length * 3

     
with open(logfile, "a") as fh2:
    #print()
    #print("# Processing:", protein_desc, file=fh2)
    print("# Average sequence length is (PROTEIN AA ungapped):", avg_sequence_length, file=fh2)
    print("# Average sequence length is (NUCLEOTIDE NUC):", avg_sequence_length_nt, file=fh2)
#end with
fh2.close()
#print("# Average sequence length is:", avg_sequence_length)


# Grab the mrna transcript
# iterate over the possible proteins from it.
# Does any of them match a protein within the protein file? if so, have a dict store the transcript, protein pair.
# output the fasta
#{} = {"NT ID -- Short form i.e. the id": {"NT ID -- Full version i.e. the description": NT_ID, 
  #             "Codon Sequence": SEQ, 
  #             "Protein ID": PROTID, 
  #             "PROTEIN_SEQ": Protein_seq}}


# Grab the protein, can I find a match in any of the mRNA transcript? if so how many?
with open(PROTEIN_FASTA, "r") as prot_handle:
    for n, record in enumerate(SeqIO.parse(prot_handle, "fasta")):
        # Grab protein data.
        protein_id = record.id 
        protein_desc = record.description
        protein_seq = record.seq
       
        #print() 
        # Send to log 
        with open(logfile, "a") as fh2:
            print("\n", file=fh2)
            print("# Processing:", protein_desc, file=fh2)
        #end with
        fh2.close()

        if "LOW QUALITY PROTEIN" in str(protein_desc) or "partial" in str(protein_desc):
            with open(logfile, "a") as fh2:
                print("# Skipping this sequence due to quality issues", record.description, file = fh2) 
            #end with
            fh2.close()
            continue
        #end if

        species = ""

        if "[" in protein_desc:
            species = record.description.split("[")[1].replace("]", "")
        #end if


        with open(logfile, "a") as fh2:
            print("# Species:", species, file=fh2)
        #end with
        fh2.close()

        # Threshold here based on protein length.
        #if int(avg_sequence_length) > 0:
        #    prot_threshold = avg_sequence_length / 2
        #    if len(protein_seq.ungap("-")) < prot_threshold:
        #        log("# Skipping protein sequence, failed to meet minimum length requirement", logfile)
        #        continue 
        #    #end if
        #    #continue
        ##end if 

        # Heavy lifting here.
        tx_record = Process (protein_desc, protein_seq, TRANSCRIPTS_FASTA, species, int(avg_sequence_length_nt))

        if type(tx_record) != str:
            with open(logfile, "a") as fh2:
                print("# Match:", tx_record.description, "\n", file=fh2)
            #end with
            fh2.close()
            results.append(tx_record)
        else:
            # No match...
            with open(logfile, "a") as fh2:
                print("# -- NO Match -- \n", file=fh2)
            #end with
            fh2.close()
            no_match.append(protein_desc)
            #pass
    #end for
#end with

# Bad binary operator fix here?

#for item in results:
#    a = str(item.description).split(" ")
#    #item.description = "_".join([a[0], a[1], a[2]])
#    item.id = "_".join([a[0], a[1], a[2]])
#    item.description = ""

# Write out records
with open(logfile, "a") as fh2:
    print("# Writing data to:", OUTPUT, file=fh2)
#end with
fh2.close()

# Write
SeqIO.write(results, OUTPUT, "fasta")

"""
with open(OUTPUT, "a") as fh:
    #fh.write("")
    for item in results:
         fh.write(">" + str(item.description) + "\n" + str(item.seq) + "\n")
fh.close()
"""

# Report on no matches, this needs to be written to a log

# Send to log 
with open(logfile, "a") as fh2:
    print("--- The following had no matches", file=fh2)
    for item in no_match: 
        print(item, file=fh2)
#end with
fh2.close()


#print("--- The following had no matches")
#print("Total:", len(no_match))
#for item in no_match:
#    print(item)


# =============================================================================
# End of file    
# =============================================================================
