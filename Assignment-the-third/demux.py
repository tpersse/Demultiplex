import argparse
import gzip

def getArgs():
	parser = argparse.ArgumentParser()
	# parser.add_argument(
	# 	"-i", help="input index filename", type=str, required=True
	# )
	parser.add_argument(
		"--index1", help="barcodes file 1", type=str, required=True
	)
	parser.add_argument(
		"--index2", help="barcodes file 2", type=str, required=True
	)
	parser.add_argument(
		"--read1", help="read file 1", type=str, required=True
	)
	parser.add_argument(
		"--read2", help="read file 2", type=str, required=True
	)

	return parser.parse_args()

args = getArgs()
# i = args.i
index1 = args.index1
index2 = args.index2
read1 = args.read1
read2 = args.read2

files = (index1, index2, read1, read2)

CompBaseDict = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}

########### DEFINING FUNCTIONS ###########

def reverse_comp(sequence):
	'''This function takes a DNA sequence as an argument, and returns its reverse complement'''
	revcomp = ''
	i = -1
	seqlist = sequence.split()
	for ch in seqlist[i]:
		revcomp += CompBaseDict.get(ch)
		i+=(-1)
	return (revcomp[::-1])

def convert_phred(letter: str, ch=33):
    """Converts a single character into a phred score"""
    return ord(letter) - ch

def mean_qual_score(phred_score):
    x = 0
    for ch in phred_score:
        x += (convert_phred(ch))
    return (x / len(phred_score))

########### ESTABLISHING GLOBAL VARIABLES, BUILDING UP THE BASE OF THE BODY TO COME ###########

index_list = ['GTAGCGTA', 'CGATCGAT', 'GATCAAGG', 'AACAGCGA', 'TAGCCATG', 'CGGTAATC', 'CTCTGGAT', 'TACCGGAT', 'CTAGCTCA', 'CACTTCAC', 'GCTACTCT', 'ACGATCAG', 'TATGGCAC', 'TGTTCCGT', 'GTCCTAAG', 'TCGACAAG', 'TCTTCGAC', 'ATCATGCG', 'ATCGTGGT', 'TCGAGAGT', 'TCGGATTC', 'GATCTTGC', 'AGAGTCCA', 'AGGATAGC']

# generating files
revcomp_barcodes = []
for barcode in index_list:
	revcomp_barcodes.append(reverse_comp(barcode))

r1_outfile_dict = {}
r2_outfile_dict = {}

for barcode in index_list:
	file = 'DemuxFilesOut/r1' + barcode + '.fastq'
	r1_outfile_dict(barcode) = open(file, 'w')
r1_outfile_dict[r1low] = open('DemuxFilesOut/r1lowQ.fastq')
r1_outfile_dict[r1hopped] = open('DemuxFilesOut/r1hopped.fastq')

for barcode in revcomp_barcodes:
	file = 'DemuxFilesOut/r2' + barcode + '.fastq'
	r2_outfile_dict(barcode) = open(file, 'w')
r2_outfile_dict[r2low] = open('DemuxFilesOut/r2lowQ.fastq')
r2_outfile_dict[r2hopped] = open('DemuxFilesOut/r2hopped.fastq')

########### SORTING AND BINNING FROM FILES ###########

first = True
i = 0
with gzip.open (index1, "rt") as i1, (index2, "rt") as i2, (read1, "rt") as r1, (read2, "rt") as r2:
	while True:
		i1line = i1.readline().strip()
		i2line = i1.readline().strip()
		r1line = i1.readline().strip()
		r2line = i1.readline().strip()
		i += 1
		if i1line == '': # end of file, means end the loop
			break
		if i%4 == 0:
			header_r1 = r1line
			header_r2 = r2line
		if i%4 == 1:
			bc1 = i1line
			read1 = r1line
			read2 = r2line
			bc2 = b2line
		if i%4 == 2:
			spacer1 = r1line
			spacer2 = r2line
		if i%4 == 3:
			qualscore_r1 = r1line
			qualscore_r2 = r2line
			qualscore_i1 = i1line
			qualscore_i2 = i2line
			if 'N' not in qualscore_i1 or qualscore_i2:

				# order here? two screens for quality, how to couple them or order them?

				if mean_qual_score(qualscore_i1) or mean_qual_score(qualscore_i2) =< 30:
					
					# checks to see if index 1 is equal to the reverse complement of index 2
					
					if bc1 == bc2(i2line):	
					r1_outfile_dict[bc1].write(header_r1, '\n', read1, '\n', spacer1, '\n', qualscore_r1, '\n')
					r2_outfile_dict[bc2].write(header_r2, '\n', read2, '\n', spacer2, '\n', qualscore_r2, '\n')
					
					# if index 1 is not equal to the reverse complement of index 2, adds each to respective index hopped files

					else:
					r1_outfile_dict[r1hopped].write(header_r1, '\n', read1, '\n', spacer1, '\n', qualscore_r1, '\n')
					r2_outfile_dict[r2hopped].write(header_r2, '\n', read2, '\n', spacer2, '\n', qualscore_r2, '\n')

			# if the reads do not adhere to either of the quality score checks (N present and greater than 30 on average), add reads to the low files. 

			else:
				r1_outfile_dict[r1low].write(header_r1, '\n', read1, '\n', spacer1, '\n', qualscore_r1, '\n')  # need to create these files
				r2_outfile_dict[r2low].write(header_r2, '\n', read2, '\n', spacer2, '\n', qualscore_r2, '\n')
	
	
	for file in r1_outfile_dict.values():
		file.close()
	for file in r2_outfile_dict.values():
		file.close()
	for file in [i1, r1, r2, i2]:
		file.close()
			

