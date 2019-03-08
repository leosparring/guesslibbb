import subprocess
import argparse

###--------- GLOBAL VARIABLES ---------###

# Specify the least number of sequences in the out-files
NR_OF_SEQUENCES_THRESHOLD = 1000

# Specify the quality threshold you want these sequences to have
# If you don't get NR_OF_SEQUENCES seqs of this quality, the quality threshold will
# lower in order to accomodate that
QUALITY_THRESHOLD_START = 40

# If this is set to True, the output file will replace the input file, the initial
# data will thus be lost. (be careful)
REPLACE_INPUT = False

###--------- FUNCTIONS ---------###

def check_format_and_remove_low_quality_reads_pair(fastq_file_1, fastq_file_2):
	'''
	This function writes a file with reads over a certain
	average Q-threshold. It then checks if it has at least
	NR_OF_SEQUENCES_THRESHOLD sequences. If it doesnt, it lowers the
	threshold by 1, and loops through the file again to append the file
	with sequences of meeting this lower quality threshold. It also makes
	sure to only include a sequence in its out files if the corresponding
	sequence in the other file also meets the quality threshold.
	'''

	# Naming outfiles and initializing variables
	pairchecked_f1 = "pairchecked_" + fastq_file_1
	pairchecked_f2 = "pairchecked_" + fastq_file_2
	quality_threshold = QUALITY_THRESHOLD_START
	proper_formats = False
	fin1_avg_quality = 0
	fin2_avg_quality = 0
	fin_avg_quality = 0
	fout1_avg_quality = 0
	fout2_avg_quality = 0
	fout_avg_quality = 0
	nr_of_seqs_fout = 0
	nr_of_seqs_fin = 0
	seq_ids = [] # A list of sequences stored in the outfiles

	proper_format_f1 = proper_fastq_format(fastq_file_1) # Checking if f1 and f2 are
	proper_format_f2 = proper_fastq_format(fastq_file_2) # proper formatted fastq files

	if (proper_format_f1 and proper_format_f2):
		proper_formats = True

		fout_1 = open(pairchecked_f1, "w+") # Initializing outfile 1
		fout_1.close()
		fout_2 = open(pairchecked_f2, "w+") # Initializing outfile 2
		fout_2.close()

		first_iteration = True
		while (nr_of_seqs_fout < NR_OF_SEQUENCES_THRESHOLD):
			with open(fastq_file_1) as fin_1, open(fastq_file_2) as fin_2:
				counter = 0
				for line in fin_1:
					if (counter % 4 == 0): # Reading first line of seq in f1 and f2
						first_line_f1 = line
						first_line_f2 = fin_2.readline()
					if (counter % 4 == 1): # Reading second line in f1 and f2
						second_line_f1 = line
						second_line_f2 = fin_2.readline()
					if (counter % 4 == 2): # Reading third line in f1 and f2
						third_line_f1 = line
						third_line_f2 = fin_2.readline()
					if (counter % 4 == 3): # Reading fourth line in f1 and f2
						fourth_line_f1 = line
						fourth_line_f2 = fin_2.readline()
						if first_line_f1 not in seq_ids:  # Making sure we haven't added the seq already
						# Calculating the quality of the sequences
							phred_quality_f1 = calculate_phred_quality(fourth_line_f1)
							phred_quality_f2 = calculate_phred_quality(fourth_line_f2)
							if first_iteration:
								# Calculating the average quality of the in file seqs
								fin1_avg_quality += phred_quality_f1
								fin2_avg_quality += phred_quality_f2
								nr_of_seqs_fin += 1
							if (phred_quality_f1 > quality_threshold and phred_quality_f2 > quality_threshold):
								seq_ids.append(first_line_f1)
								# Calculating the average quality of the out file seqs
								fout1_avg_quality += phred_quality_f1
								fout2_avg_quality += phred_quality_f2
								nr_of_seqs_fout += 1
								with open(pairchecked_f1, "a") as fout_1: # Appending
									fout_1.write(first_line_f1)			# the sequence to
									fout_1.write(second_line_f1)		# outfile 1
									fout_1.write(third_line_f1)
									fout_1.write(fourth_line_f1)
								with open(pairchecked_f2, 'a') as fout_2: # Appending the
									fout_2.write(first_line_f2)				# sequence to
									fout_2.write(second_line_f2)			# outfile 2
									fout_2.write(third_line_f2)
									fout_2.write(fourth_line_f2)
					
					counter += 1 # Nr of lines read in file1 and file2
			quality_threshold = quality_threshold - 1 # Lowering quality
			first_iteration = False						# for next iteration

		if REPLACE_INPUT: # Replacing the input files with the output files
			subprocess.run(["mv", pairchecked_f1, fastq_file_1])
			subprocess.run(["mv", pairchecked_f2, fastq_file_2])

		# Calculating the average Q-values
		fin1_avg_quality = fin1_avg_quality / nr_of_seqs_fin
		fin2_avg_quality = fin2_avg_quality / nr_of_seqs_fin
		fin_avg_quality = (fin1_avg_quality + fin2_avg_quality)/2
		fout1_avg_quality = fout1_avg_quality / nr_of_seqs_fout
		fout2_avg_quality = fout2_avg_quality / nr_of_seqs_fout
		fout_avg_quality = (fout1_avg_quality + fout2_avg_quality)/2

		print(f'All of the sequences in the out files '
			f'have average qualities above {quality_threshold}.\n'
			f'Average quality of {fastq_file_1}: {fin1_avg_quality}\n'
			f'Average quality of corresponding outfile to f1: {fout1_avg_quality}\n'
			f'Average quality of {fastq_file_2}: {fin2_avg_quality}\n'
			f'Average quality of corresponding outfile to f2: {fout2_avg_quality}')

	return proper_formats, nr_of_seqs_fout, round(fout_avg_quality, 2), round(fin_avg_quality, 2)


def proper_fastq_format(fastq_file):
	'''
	This function returns True if the fastq file is properly formatted.
	'''

	fastq_format = False
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
					fastq_format = True

			except UnicodeDecodeError:
				fastq_format = False
				print("You had a UnicodeDecodeError.")
	except IOError:
		fastq_format = False
		print("You had an IOError.")

	return fastq_format


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


###--------- MAIN ---------###

def main():

	parser = argparse.ArgumentParser()
	parser.add_argument("-f1", "--file_1", required=True, help="FASTQ file one")
	parser.add_argument("-f2", "--file_2", required=True, help="corresponding FASTQ of paired sequences")
	args = parser.parse_args()

	(proper_formats, nr_of_seqs_fout, fout_avg_quality,
		fin_avg_quality) = check_format_and_remove_low_quality_reads_pair(args.file_1,
																			args.file_2)

	print(f'The two input files are properly formatted: {proper_formats}\n'
		f'The number of sequences in your subsetted files are: {nr_of_seqs_fout}')

if __name__ == "__main__":
	main()