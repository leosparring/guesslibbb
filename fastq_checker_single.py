import argparse
import subprocess

###--------- GLOBAL VARIABLES ---------###

# Specify the least number of sequences in the out-file
NR_OF_SEQUENCES_THRESHOLD = 1000

# Specify the quality threshold you want these sequences to have
# If you don't get NR_OF_SEQUENCES seqs of this quality, the quality will
# lower in order to accomodate that
QUALITY_THRESHOLD_START = 40

# If this is set to True, the output file will replace the input file, the initial
# data will thus be lost. (be careful)
REPLACE_INPUT = True

###--------- FUNCTIONS ---------###

def proper_fastq_format(fastq_file):
	'''
	This function returns True if the fastq file is properly formatted.
	'''

	proper_fastq_format = False
	fastq_decision = 0
	counter = 0

	try:
		with open(fastq_file) as f_in:
			try:
				for line in f_in:
					if (counter % 4 == 0): # Checking first line
						if (line.startswith('@')):
							fastq_decision += 1
					if (counter % 4 == 1): # Checking second line
						second_line = line
						if all(c in 'ATGCNatgcn\n' for c in list(line)):
							fastq_decision += 1
					if (counter % 4 == 2): # Checking third line
						if (line.startswith('+')):
							fastq_decision += 1
					if (counter % 4 == 3): # Checking fourth line
						if (len(line) == len(second_line)):
							fastq_decision += 1
					counter += 1
				if (fastq_decision == counter):
					proper_fastq_format = True

			except UnicodeDecodeError:
				fastq_format = False
				print("You had a UnicodeDecodeError.")
	except IOError:
		fastq_format = False
		print("You had an IOError.")

	return proper_fastq_format


def calculate_phred_quality(phred_string):
	'''
	This function calculates the average Q-phred based on an ascii-string,
	using the phred ascii-33 system.
	'''

	phred_char_list = list(phred_string)
	phred_quality = 0

	for i in range(len(phred_char_list)):
		phred_quality = phred_quality + ord(phred_char_list[i])-33

	phred_quality = phred_quality/len(phred_char_list)
	
	return phred_quality


def check_format_and_remove_low_quality_reads_single(fastq_file):
	'''
	This function writes a file with reads over a certain
	average Q-threshold. It then checks if it has at least
	NR_OF_SEQUENCES_THRESHOLD sequences. If it doesnt, it lowers the
	threshold by 1, and loops through the file again to append the file
	with sequences of meeting this lower quality threshold.
	'''

	quality_threshold = QUALITY_THRESHOLD_START
	q_fastq_file = "singlechecked_" + fastq_file
	nr_of_sequences_out = 0
	fin_avg_quality = 0
	fout_avg_quality = 0
	nr_of_seqs_fin = 0
	seq_ids = [] # A list of sequences stored in the outfiles

	proper_format = proper_fastq_format(fastq_file)

	if proper_format:
		f = open(q_fastq_file, "w+") # Initializing the outfile.
		f.close()

		first_iteration = True
		while (nr_of_sequences_out < NR_OF_SEQUENCES_THRESHOLD):
			with open(fastq_file) as f_in:
				counter = 0
				for line in f_in:
					if (counter % 4 == 0): # Reading the first line of seq
						first_line = line
					if (counter % 4 == 1): # Reading the second line
						second_line = line
					if (counter % 4 == 2): # Reading the third line
						third_line = line
					if (counter % 4 == 3): # Reading th quality line
						fourth_line = line
						if (first_line not in seq_ids):
							phred_quality = calculate_phred_quality(fourth_line)
							if first_iteration:
								# Calculating the avg quality, the first loop
								fin_avg_quality += phred_quality
								nr_of_seqs_fin += 1
							if (phred_quality > quality_threshold):
								seq_ids.append(first_line) # adding to seq_ids
								fout_avg_quality += phred_quality
								nr_of_sequences_out += 1
								with open(q_fastq_file, "a") as f_out: 
									f_out.write(first_line)
									f_out.write(second_line)
									f_out.write(third_line)
									f_out.write(fourth_line)
						
					counter += 1 # Nr of lines read in the file
			quality_threshold = quality_threshold - 1 # Reducing quality threshold for 
			first_iteration = False 					# next iteration

		fin_avg_quality = fin_avg_quality / nr_of_seqs_fin
		fout_avg_quality = fout_avg_quality / nr_of_sequences_out

		if REPLACE_INPUT:
			subprocess.run(["mv", q_fastq_file, fastq_file])

		print(f'All of the sequences in the out file '
			f'have an average quality above {quality_threshold}.')

	return proper_format, nr_of_sequences_out, round(fout_avg_quality, 2), round(fin_avg_quality, 2)

###--------- MAIN ---------###

def main():

	parser = argparse.ArgumentParser()
	parser.add_argument("-f", "--file", help="FASTQ file")
	args = parser.parse_args()

	(proper_format,
		nr_of_sequences_fout,
		fout_avg_quality,
		fin_avg_quality) = check_format_and_remove_low_quality_reads_single(args.file)
	
	print(f'The input file was in proper fastq format: {proper_format}\n'
		f'It contains {nr_of_sequences_fout} nr of sequences '
		f'and average quality Q = {fout_avg_quality}.\n'
		f'The average quality of {args.file} Q = {fin_avg_quality}.')


if __name__ == "__main__":
	main()
