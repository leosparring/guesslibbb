import subprocess
import csv
import argparse
import time
import numpy as np
import scipy.stats as stats

###--------- GLOBAL VARIABLES ---------###

# Specify from what line of the user's transcripts
# that sequences should be aligned from. Index from zero.
START_FROM_SEQUENCE_NR = 0

# Specify from how many paired alignments data
# should be collected from.
NR_OF_PAIRS = 100

# Specify the maximum nr of alignments to perform before giving up on determining
# the library type.
MAX_BLATS = 1000

# Specify the significance threshold for the p-values.
SIGNIFICANT_P = 0.05

# Dictionary of the different library types.
LIB_TYPE_DICT = {
	0 : "I",
	1 : "O"
}

###--------- TRAINED STATS ---------###

# This array contains typical data of an inwards library, the first element
# is the nr of reads that were inward oriented, the second is the total
# nr of collected pairs
INWARD_DATA = ([1000, 5])

# This array contains typical data of an outwards library, the first element
# is the nr of reads that were outward oriented, the second is the total
# nr of collected pairs
OUTWARD_DATA = ([900, 100])


###--------- FUNCTIONS ---------###

def guesslib_genomic_pair(ref, user_transcripts_f1, user_transcripts_f2):
	'''
	Finds NR_OF_PAIRS pairs and
	determines if the library is forward, reverse or unstranded
	'''

	# Initializing variables
	start_time = time.time()
	# Naming tmp fasta file according to input file name
	tmp_fasta_1 = "tmp/tmp_blat_inpt_1_" + user_transcripts_f1
	lib_type = "N/A"
	succesful_lib_determination = False
	read_start = START_FROM_SEQUENCE_NR * 4
	seqs_searched = 0
	collected_pairs = 0
	o_and_f = 0
	o_and_r = 0
	i_and_f = 0
	i_and_r = 0
	invalid_orientation = 0
	p_value = 0.99

	try:
		with open (user_transcripts_f1) as f_in:
			for i in range(read_start): # Skipping to START_FROM_SEQUENCE_NR
				next(f_in)
			while ((collected_pairs < NR_OF_PAIRS or p_value > SIGNIFICANT_P)
				and seqs_searched < MAX_BLATS):
				pair_type = "N/A" # Resetting variables for new search
				nr_of_results_R1 = 0
				nr_of_results_R2 = 0

				print(f"Collected pairs: {collected_pairs}. Collecting at least {NR_OF_PAIRS} pairs.")
				seq_id_R1 = write_tmp_fasta(tmp_fasta_1, f_in)
				(seq_start_R1, seq_end_R1, nr_of_results_R1,
					ref_transcript_id_R1) = run_blat(ref, tmp_fasta_1)
				seqs_searched += 1

				if (nr_of_results_R1 == 1):
					(seq_start_R2, seq_end_R2, nr_of_results_R2,
	            		ref_transcript_id_R2) = blat_corresponding_seq_in_f2(ref,
	                                                user_transcripts_f2,
	                                                seq_id_R1)

					if (nr_of_results_R2 == 1 and ref_transcript_id_R1 == ref_transcript_id_R2):
						pair_type = pair_analysis(seq_start_R1,
					                                seq_end_R1,
					                                seq_start_R2,
					                                seq_end_R2)
						collected_pairs += 1
					
						# Incrementing corresponding pair type
						if (pair_type == "OF"):
							o_and_f += 1
							print("Added to OF.\n")
						elif (pair_type == "OR" ):
							o_and_r += 1
							print("Added to OR.\n")
						elif (pair_type == "IF"):
							i_and_f += 1
							print("Added to IF.\n")
						elif (pair_type == "IR"):
							i_and_r += 1
							print("Added to IR.\n")
						else:
							invalid_orientation += 1
							print("Added to invalid orientation.\n")

						if (collected_pairs > NR_OF_PAIRS-1):
							(lib_type,
								p_value) = get_libtype_and_pvalue(o_and_f, o_and_r,
																	i_and_f, i_and_r,
																	collected_pairs)

	except StopIteration:
		pass

	subprocess.run(['rm', tmp_fasta_1]) # Removing the tmp FASTA file

	if (seqs_searched == MAX_BLATS or lib_type == "N/A"):
		print(f'Guesslib could not determine the library type.\n'
			f'Tried to align {MAX_BLATS} sequences. With this, '
			f'{collected_pairs} pairs were collected.\n'
			f'Of these, there were {o_and_f} OFs, {o_and_r} ORs, '
			f'{i_and_f} IFs, {i_and_r} IRs.\n'
			f'Possibly, the wrong reference sequence was used.')
	else:
		succesful_lib_determination = True

		print(f'\nThe library type determination was succesful\n'
			f'OF: {o_and_f}, OR: {o_and_r}, '
			f'IF: {i_and_f}, IR: {i_and_r}, '
			f'Neither: {invalid_orientation}.\n'
			f'The most likely lib type is {lib_type}.')

	end_time = time.time()
	analysis_time = round(((end_time - start_time)/60), 2)

	print(f'Guesslib took {int(analysis_time)} minutes and'
		f' {int((analysis_time-int(analysis_time))*60)} seconds to do the analysis.')

	return (lib_type, o_and_f, o_and_r, i_and_f, i_and_r, collected_pairs,
			succesful_lib_determination, analysis_time)

def get_libtype_and_pvalue(o_and_f, o_and_r, i_and_f, i_and_r, collected_pairs):
	'''
	This function makes a decision about the 
	library types based on the nr of different pair types.
	'''

	lib_type = "N/A"
	p_values = []
	p_value = 0.99

	# First we check inward, outward orientation
	inwards = i_and_f + i_and_r
	outwards = o_and_f + o_and_r

	p_values.append(stats.fisher_exact([INWARD_DATA,
		[inwards, collected_pairs-inwards]])[1])
	p_values.append(stats.fisher_exact([OUTWARD_DATA,
		[outwards, collected_pairs-outwards]])[1])

	# This is the lowest of the two p values.
	p_value = min(p_values)
	# Retrieving orientation corresponding to highest p_val, from dictionary
	lib_type = LIB_TYPE_DICT.get(p_values.index(max(p_values))) 
	print(f'The list of p-values ["I", "O"] = {p_values}\n'
		f'The library type appears to be: {lib_type}.\n'
		f'The second to highest p-value is: {p_value}')

	return lib_type, p_value

def pair_analysis(seq_start_R1, seq_end_R1, seq_start_R2, seq_end_R2):
	'''
	This function determines inward or outward and reverse
	or forward orientation based on the input of two reads.
	'''

	pair_type = "N/A"

	# This is the distance between the reads' midpoints, but it can be negative
	distance_r2_r1 = (((seq_start_R2 + seq_end_R2)/2) -
						((seq_start_R1 + seq_end_R1)/2))

	if (seq_start_R1 < seq_end_R1 and seq_start_R2 > seq_end_R2): # If R1 is F and R2 is R
		if (seq_start_R1 < seq_start_R2): # If R1 starts upstream of R2
			pair_type = "IF"
		elif (seq_start_R1 > seq_start_R2): # If R1 starts downstream of R2
			pair_type = "OF"
	elif (seq_start_R1 > seq_end_R1 and seq_start_R2 < seq_end_R2): # IF R1 is R and R2 is F
		if (seq_start_R1 > seq_start_R2): # If R1 starts downstream of R2
			pair_type = "IR"
		elif (seq_start_R1 < seq_start_R2): # If R1 starts upstream of R2
			pair_type = "OR"

	print(f'(Midpoint of R2) - (Midpoint of R1) is: {distance_r2_r1}')

	return pair_type

def blat_corresponding_seq_in_f2(ref, user_transcripts, seq_id_R1):
    '''
    This function finds the corresponding seq to the id of the input in the
    reference fastq file specified as ref. It then performs a blat and returns
    the start, end and e-value of this hit.
    '''

    seq_id_R2 = "N/A"
    iterations = 0
    ref_transcript_id = "N/A"
    tmp_in_2 = "tmp/tmp_blat_inpt_2_" + user_transcripts # Naming temp fasta file

    with open(user_transcripts) as f_in:
        first_line = '>' + f_in.readline()[1:] # First line of FASTA file
        seq_id_R2 = first_line.split(" ")[0]

        while (seq_id_R2 != seq_id_R1): # Looking for match to seq_id_R1
            next(f_in)
            next(f_in)
            next(f_in)
            first_line = '>' + f_in.readline()[1:] # First line of FASTA
            seq_id_R2 = first_line.split(" ")[0]
            iterations += 1

        # Writing out tmp FASTA
        with open(tmp_in_2, "w+") as f_out:
            f_out.write(first_line)
            f_out.write(f_in.readline())

    print("Now blatting the corresponding sequence in file 2.")

    (seq_start, seq_end, nr_of_results,
        ref_transcript_id) = run_blat(ref, tmp_in_2)

    subprocess.run(['rm', tmp_in_2]) # Removing the tmp file

    return (seq_start, seq_end,
            nr_of_results, ref_transcript_id)

def write_tmp_fasta(tmp_fasta, fastq):
	'''
	Writes a temporary fasta file, taking in a fastq-file.
	'''

	with open(tmp_fasta, "w+") as f_out: # Creating FASTA tmp
		first_line = '>' + fastq.readline()[1:]
		seq_id = first_line.split(" ")[0]
		print(seq_id)
		f_out.write(first_line)  # Writes the header into f_out.
		f_out.write(fastq.readline()) # Writes the sequence into f_out.
		next(fastq) # Skipping past the +, and phred lines.
		next(fastq)

	return seq_id

def run_blat(ref, seq):
    ''' This function calls blat with a single sequence seq.fa as
        input against the reference ref.fa and returns seq_identity,
        seq_start, seq_end'''

    nr_of_results = 0
    seq_identity = "0"
    seq_start = "0"
    seq_end = "0"
    second_hit_e_value = "0"
    e_value = "0"
    ref_transcript_id = "N/A"
    tmp_rslt = seq + "_rslt"
    ref = 'reference_sequences/' + ref

    # Blatting, and creating tmp blat result file
    subprocess.run(['./blat', ref, seq, '-out=blast8', tmp_rslt])
    for line in open(tmp_rslt):
        nr_of_results += 1 # Each row represents an alignment
    if (nr_of_results > 0):
        with open(tmp_rslt, 'r') as f_in:
            first_hit = next(csv.reader(f_in, delimiter='\t'))
            ref_transcript_id = first_hit[1]
            seq_identity = first_hit[2]
            seq_start = first_hit[8]
            seq_end = first_hit[9]
            e_value = first_hit[10]
            if (nr_of_results > 1):
                second_hit = next(csv.reader(f_in, delimiter='\t'))
                second_hit_e_value = second_hit[10]

    print(f'the identity is: {seq_identity}')
    print(f'the E value of the first hit is: {e_value}')
    print(f'the seq start is: {seq_start}')
    print(f'the seq end is: {seq_end}')
    print(f"the second hit's E value is: {second_hit_e_value}")
    print(f'the reference transcript id is: {ref_transcript_id}')
    print(f'the number of results for this blat is: {nr_of_results}\n')

    subprocess.run(['rm', tmp_rslt]) # Removing the tmp blat rslt file

    return (int(seq_start), int(seq_end),
            nr_of_results, ref_transcript_id)

###--------- MAIN ---------###

def main():
	
	parser = argparse.ArgumentParser()
	parser.add_argument("-r", "--reference", required=True, help="reference sequence")
	parser.add_argument("-f1", "--file_1", required=True, help="FASTQ with read 1 sequences")
	parser.add_argument("-f2", "--file_2", required=True, help="FASTQ with paired read 2 sequences")
	args = parser.parse_args()

	guesslib_genomic_pair(args.reference, args.file_1, args.file_2)
		
if __name__ == "__main__":
	main()
 
