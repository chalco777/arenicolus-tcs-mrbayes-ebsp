#! /usr/bin/env python3
#
# Run ./collapse_cleanup.py file-in filetype(nexus or fasta) then optional sequence_cleaner flags
# Lauren Chan 06 July 2013
# Updated for Python 3: prints, input(), imports and safer file handling.

import sys

try:
    from Bio import SeqIO
except ImportError:
    sys.exit("Biopython is required (pip install biopython)")

def sequence_cleaner(filename, filetype, min_length=0, por_n=100):
    ft = filetype.lower()
    if "nex" in ft:
        print("Infile is in Nexus format")
        dat = SeqIO.parse(filename, "nexus")
    elif "fa" in ft:
        print("Infile is in fasta format")
        dat = SeqIO.parse(filename, "fasta")
    else:
        print("Unclear what filetype this is... quitting")
        sys.exit(0)
    
    ### POACHED FROM SEQUENCE_CLEANER.PY ###
    # create our hash table to add the sequences
    sequences = {}
 
    # Using the biopython parse we can read our input
    for seq_record in dat:
        # Take the current sequence
        sequence = str(seq_record.seq).upper()
        # Check if the current sequence is according to the user parameters
        if len(sequence) >= int(min_length) and (float(sequence.count("N")) / float(len(sequence))) * 100 <= float(por_n):
            # If the sequence passed the test and it isn't in the hash table,
            # the sequence and its id are going to be in the hash
            seq_id = getattr(seq_record, "id", None) or getattr(seq_record, "name", None) or ""
            if sequence not in sequences:
                sequences[sequence] = seq_id
            else:
                # If it is already in the hash table, concatenate the ID
                sequences[sequence] += "*" + seq_id
 
    # Write the clean sequences
    clear_filename = "clear_" + filename + ".fa"
    with open(clear_filename, "w", encoding="utf-8") as output_file:
        for sequence, ids in sequences.items():
            output_file.write(">{}\n{}\n".format(ids, sequence))
 
    print("CLEAN!!!\n{} holds your haplotypes.\n".format(clear_filename))
    ### -- END POACHING -- ###

    print("Now we'll clean up the taxon labels\n")
    outroot = input("Root for outfiles: ")
    make_nex = input("Make a Nexus file? y/n: ")
    reduce_hapnames(clear_filename, outroot, make_nex)
    

def reduce_hapnames(infile, outroot, make_nex):
    tosort = SeqIO.parse(infile, "fasta")  # Get data from fasta
    sorted_names = sorted(x.name for x in tosort)
    record_index = SeqIO.index(infile, "fasta")
    records = (record_index[id] for id in sorted_names) 
    
    # Create output files... need to decide how to name these.
    with open(outroot + ".haps.fasta", "w", encoding="utf-8") as fastaout, \
         open(outroot + ".hapkey.txt", "w", encoding="utf-8") as namesout:
    
        for index, dat in enumerate(records):   # For each instance:
            thenames = dat.name                 # Get the long name
            seqname = thenames.split('*')       # Split by *
            seqcount = len(seqname)             # Count the number of indiv with that haplotype
        
            fastaout.write('>{}_{}\n{}\n'.format(seqname[0], seqcount, str(dat.seq)))
            namesout.write('{} {}_{}: {}\n'.format(index, seqname[0], seqcount, seqname))
    
    # Dimensions of the datamatrix
    n = index + 1
    seqlen = len(str(dat.seq))   
    
    # Whether or not we make a nexus file
    if "y" in make_nex or "Y" in make_nex:
        write_nexus(outroot, n, seqlen, infile)
        print("Generating Nexus matrix with {} seqs and {} basebairs".format(n, seqlen))
    else:
        print("Not generating a Nexus file")


# Writing the nexus file
def write_nexus(outroot, ntax, nchar, infile):
    cleanedfasta = SeqIO.parse(outroot + ".haps.fasta", "fasta")
    with open(outroot + ".haps.nex", "w", encoding="utf-8") as nexusout:
        nexusout.write('#NEXUS\n')
        nexusout.write('[Converted from {} using reduce_convert.py, by LMC July 06, 2013]\n'.format(infile))
        nexusout.write('[{} holds the haplotype designations]\n\n'.format(outroot + '.hapkey.txt'))
        
        nexusout.write('Begin DATA;\n')
        nexusout.write('\tdimensions ntax={} nchar={};\n'.format(ntax, nchar))
        nexusout.write('\tformat datatype=DNA missing=? gap=-;\n\tmatrix\n')
        for i in cleanedfasta:
            nexusout.write('\t\t{}\t\t{}\n'.format(i.name, str(i.seq)))
        nexusout.write(';\nEND;')

# Running of script   
userParameters = sys.argv[1:]
 
try:
    if len(userParameters) == 1:
        sequence_cleaner(userParameters[0], "fasta")
    elif len(userParameters) == 2:
        sequence_cleaner(userParameters[0], userParameters[1])
    elif len(userParameters) == 3:
        sequence_cleaner(userParameters[0], userParameters[1], int(userParameters[2]))
    elif len(userParameters) == 4:
        sequence_cleaner(userParameters[0], userParameters[1], int(userParameters[2]), int(userParameters[3]))
    else:
        print("There is a problem!")
except Exception as e:
    print("There is a problem:", e)
    sys.exit(1)