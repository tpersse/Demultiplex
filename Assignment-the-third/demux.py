import argparse
import gzip

def getArgs():
	parser = argparse.ArgumentParser()
	parser.add_argument(
		"-Q", help="Average quaulity score cutoff for indices", type=float, required=True
	)
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
q = args.Q
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

# amount counters
total = 0
hopped = 0
proper = 0
lowQ = 0

index_list = ['GTAGCGTA', 'CGATCGAT', 'GATCAAGG', 'AACAGCGA', 'TAGCCATG', 'CGGTAATC', 'CTCTGGAT', 'TACCGGAT', 'CTAGCTCA', 'CACTTCAC', 'GCTACTCT', 'ACGATCAG', 'TATGGCAC', 'TGTTCCGT', 'GTCCTAAG', 'TCGACAAG', 'TCTTCGAC', 'ATCATGCG', 'ATCGTGGT', 'TCGAGAGT', 'TCGGATTC', 'GATCTTGC', 'AGAGTCCA', 'AGGATAGC']

barcodes_dictionary = {}
for barcode in index_list:
	barcodes_dictionary[barcode] = 0

# generating files
revcomp_barcodes = []
for barcode in index_list:
	revcomp_barcodes.append(reverse_comp(barcode))

r1_outfile_dict = {}
r2_outfile_dict = {}
index_counter = {}

for barcode in index_list:
	file = 'DemuxFilesOut/r1' + barcode + '.fastq'
	r1_outfile_dict[barcode] = (open(file, 'w'))
r1_outfile_dict['r1low'] = (open('DemuxFilesOut/r1lowQ.fastq', 'w'))
r1_outfile_dict['r1hopped'] = (open('DemuxFilesOut/r1hopped.fastq', 'w'))

for barcode in revcomp_barcodes:
	file = 'DemuxFilesOut/r2' + barcode + '.fastq'
	r2_outfile_dict[barcode] = (open(file, 'w'))
r2_outfile_dict['r2low'] = (open('DemuxFilesOut/r2lowQ.fastq', 'w'))
r2_outfile_dict['r2hopped'] = (open('DemuxFilesOut/r2hopped.fastq', 'w'))


########### SORTING AND BINNING FROM FILES ###########

i = 0
with gzip.open(index1, "rt") as i1, gzip.open(index2, "rt") as i2, gzip.open(read1, "rt") as r1, gzip.open(read2, "rt") as r2:
# with open(index1, "r") as i1, open(index2, "r") as i2, open(read1, "r") as r1, open(read2, "r") as r2:
	while True:
	
		#headers
		header_r1 = r1.readline().strip()
		header_r2 = r2.readline().strip()
		i1.readline()
		i2.readline()

		#seqs
		bc1 = i1.readline().strip()
		read1 = r1.readline().strip()
		read2 = r2.readline().strip()
		bc2 = i2.readline().strip()
		
		#plus
		i1.readline()
		i2.readline()
		spacer1 = r1.readline().strip()
		spacer2 = r2.readline().strip()

		#quality scores
		qualscore_r1 = r1.readline().strip()
		qualscore_r2 = r2.readline().strip()
		qualscore_i1 = i1.readline().strip()
		qualscore_i2 = i2.readline().strip()

		if header_r1 == '': # end of file, means end the loop
			break

		if bc1 in r1_outfile_dict and bc2 in r2_outfile_dict and (mean_qual_score(qualscore_i1) >= q and mean_qual_score(qualscore_i2) >= q):
				
			if bc1 == reverse_comp(bc2):
				r1_outfile_dict[bc1].write(header_r1 +'\n' + read1 +'\n' + spacer1 + '\n' + qualscore_r1 + '\n')
				r2_outfile_dict[bc2].write(header_r2 + '\n' + read2 +'\n' + spacer2 + '\n' + qualscore_r2 + '\n')
				barcodes_dictionary[bc1] += 1
				proper +=1

			else:
				r1_outfile_dict['r1hopped'].write(header_r1 +'\n' + read1 +'\n' + spacer1 + '\n' + qualscore_r1 + '\n')
				r2_outfile_dict['r2hopped'].write(header_r2 + '\n' + read2 + '\n' + spacer2 + '\n' + qualscore_r2 + '\n')
				hopped += 1

		# if the reads do not adhere to either of the quality score checks (N present and greater than 30 on average), add reads to the low files. 

		else:
			r1_outfile_dict['r1low'].write(header_r1 + '\n' + read1 + '\n' + spacer1 + '\n' + qualscore_r1 + '\n')  # need to create these files
			r2_outfile_dict['r2low'].write(header_r2 + '\n' + read2 +'\n' + spacer2 + '\n' + qualscore_r2 + '\n')
			lowQ +=1
		i+=1
	
	for file in r1_outfile_dict.values():
		file.close()
	for file in r2_outfile_dict.values():
		file.close()
	for file in [i1, r1, r2, i2]:
		file.close()

########### RUNNING STATISTICAL ANALYSES ON OUTPUT ###########

percent_hopped = (hopped / i) * 100

percent_toolow = (lowQ / i) * 100

percent_proper = (proper / i) * 100

with open('DemuxFilesOut/answers.md', 'w') as answers:
	answers.write('the total number of reads: ' + str(i) + '\n')
	answers.write('the percentage of proper reads: ' + str(percent_proper) + '%' + '\n')
	answers.write('the percentage of reads that hopped: ' + str(percent_hopped)+ '%' + '\n')
	answers.write('the percentage of reads that were too low quality: ' + str(percent_toolow)+ '%' + '\n')
	answers.write('Index Pair' + '\t' + 'Raw Count of Matched Reads' + '\t' + 'Percentage of the Whole' + '\n')
	for entry in barcodes_dictionary:
		answers.write(entry + reverse_comp(entry) + '\t' + str(barcodes_dictionary[entry]) + '\t' + str((barcodes_dictionary[entry]/i)*100) + '%' + '\n')
	answers.close()
