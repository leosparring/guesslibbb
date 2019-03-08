import subprocess
import csv
import argparse
import time
import numpy as np
import scipy.stats as stats

###--------- GLOBAL VARIABLES ---------###

# Specify from what line of the user's transcripts
# that sequences should be aligned from.
START_FROM_SEQUENCE_NR = 0

# Specify from how many alignments that data
# should be collected from at the least.
NR_OF_READS = 10

# Specify the maximum nr of alignments to perform before giving up on determining
# the library type
MAX_BLATS = 1000

# Specify the significance threshold for the p-values.
SIGNIFICANT_P = 0.05

# Specify the significance threshold for p_value corresponding to
# the null hypotheses that we want to be true. collected data looks
# like corresponding trained data.
NULL_P = 0.01

# Dictionary of the different library types.
LIB_TYPE_DICT = {
	0 : "SF",
	1 : "SR",
	2 : "U"
}

###--------- TRAINING STATS ---------###

# This array contains typical data of an unstranded paired library
# First element represents forwards, the second one is 'other types'
UNSTRANDED_DATA = ([500,500])

# This array contains typical data of a stranded paired library
# First element represents forwards or reverses, the second one is 'other types'
STRANDED_DATA = ([995, 2])

###--------- FUNCTIONS ---------###

def guesslib_single(ref, user_transcripts):
	'''
	Finds NR_OF_READS pairs and
	determines if the library is forward, reverse or unstranded
	'''

	# Initializing variables
	start_time = time.time()
	# Naming tmp fasta file according to input file name
	tmp_fasta = "tmp/tmp_blat_inpt_" + user_transcripts
	lib_type = "N/A"
	succesful_lib_determination = False
	read_start = START_FROM_SEQUENCE_NR * 4
	seqs_searched = 0
	collected_reads = 0
	forward = 0
	reverse = 0
	p_value = 0.99

	try:
		with open (user_transcripts) as f_in:
			for i in range(read_start): # Skipping to START_FROM_SEQUENCE_NR
				next(f_in)
			while ((collected_reads < NR_OF_READS or p_value > SIGNIFICANT_P)
				and seqs_searched < MAX_BLATS):
				read_type = "N/A"
				nr_of_results = 0

				print(f"Collected reads: {collected_reads}. Collecting at least {NR_OF_READS} reads.")
				write_tmp_fasta(tmp_fasta, f_in)
				(seq_start, seq_end, nr_of_results,
					ref_transcript_id) = run_blat(ref, tmp_fasta)
				seqs_searched += 1
		
				if (nr_of_results == 1):
					read_type = orientation_analysis(seq_start, seq_end)
					if (read_type == "F"):
						forward += 1
					elif (read_type == "R"):
						reverse += 1
					collected_reads += 1

					if (collected_reads > (NR_OF_READS-1)):
						(lib_type,
							p_value) = get_libtype_and_pvalue(forward, reverse,
																collected_reads)
	except StopIteration:
		pass
	
	subprocess.run(['rm', tmp_fasta]) # Removing the tmp FASTA file

	if (seqs_searched == MAX_BLATS or lib_type == "N/A"):
		print(f'Guesslib could not determine the library type.\n'
			f'Tried to align {MAX_BLATS} sequences. With this, '
			f'{collected_reads} reads were collected.\n'
			f'Possibly, the wrong reference sequence was used.')
	else:
		succesful_lib_determination = True

		print(f'\nThe library type determination was succesful\n'
			f'Forwards: {forward}, reverses: {reverse}, '
			f'The most likely lib type is {lib_type}.')

	end_time = time.time()
	analysis_time = round(((end_time - start_time)/60), 2)
	print(f'Guesslib took {int(analysis_time)} minutes and'
		f' {int((analysis_time-int(analysis_time))*60)} seconds to do the analysis.')
	
	return (lib_type, forward, reverse, collected_reads,
			succesful_lib_determination, analysis_time)

def get_libtype_and_pvalue(forward, reverse, collected_reads):
	'''
	This function calculates exact p-values for the library types and
	returns the most likely library type and the p-value of data against the
	second to most likely library type.
	'''

	p_values = 0.99
	lib_type = "N/A"
	p_values = []

	p_values.append(stats.fisher_exact([STRANDED_DATA, [forward, collected_reads-forward]])[1])
	p_values.append(stats.fisher_exact([STRANDED_DATA, [reverse, collected_reads-reverse]])[1])
	p_values.append(stats.fisher_exact([UNSTRANDED_DATA, [forward, collected_reads-forward]])[1])

	if (max(p_values) > NULL_P):
		p_value = sorted(p_values, reverse=True)[1] # This is the second to highest p_value.

		lib_type = LIB_TYPE_DICT.get(p_values.index(max(p_values)))
	
	print(f'The list of p-values ["SF", "SR", "U"] = {p_values}\n'
		f'The library type appears to be: {lib_type}.\n'
		f'The second to highest p-value is: {p_value}\n') 

	return lib_type, p_value


def orientation_analysis(seq_start, seq_end):
	'''
	This function determines the orientation of one read.
	'''

	read_type = "N/A"

	print(f'The read start and end is ({seq_start}, {seq_end}).')
	if (seq_start < seq_end):
		read_type = "F"
		print("The read type is F.\n")
	elif (seq_start > seq_end):
		read_type = "R"
		print("The read type is R.\n")

	return read_type

def write_tmp_fasta(tmp_fasta, fastq):
	'''
	Writes a temporary fasta file, taking in a fastq-file
	'''

	with open(tmp_fasta, "w+") as f_out: # Creating FASTA tmp
		first_line = '>' + fastq.readline()[1:]
		seq_id = first_line.split(" ")[0]
		print(seq_id)
		f_out.write(first_line)  # Writes the header into f_out.
		f_out.write(fastq.readline()) # Writes the sequence into f_out.
		next(fastq) # Skipping past the +, and phred lines.
		next(fastq)

def run_blat(ref, seq):
    '''
    This function calls blat with a single sequence seq.fa as
    input against the reference ref.fa and returns seq_identity,
    seq_start, seq_end
    '''
        
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
	parser.add_argument("-f", "--user_transcripts", required=True, help="FASTQ library")
	args = parser.parse_args()

	guesslib_single(args.reference, args.user_transcripts)

if __name__ == "__main__":
	main()
