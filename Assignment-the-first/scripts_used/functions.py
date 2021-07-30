def reverse_comp(sequence:str)->str:
	'''This function takes a DNA sequence as an argument, and returns its reverse complement'''
	return

if __name__ == "__main__":
	assert reverse_comp('AAAAA') == 'TTTTT'
	assert reverse_comp('GAGA') == 'TCTC'


def single_N_fix(index:str)->str:
	'''This function takes the combined indices for a given read, and returns a corrected version wherein the N present is replaced by the correct nt'''
	return 
	
if __name__ == "__main__":
	assert single_N_fix('ANA+TGT') == 'ACA+TGT'
	assert single_N_fix('NGGG+CCCT') == 'AGGG+CCCT'

