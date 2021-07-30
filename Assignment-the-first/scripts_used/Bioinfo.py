def validate_base_seq(seq: str,RNAflag: bool=False) -> bool:
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    DNAbases = set('ATGCatcg')
    RNAbases = set('AUGCaucg')
    return set(seq)<=(RNAbases if RNAflag else DNAbases)

if __name__ == "__main__":
    assert validate_base_seq("AATAGAT") == True, "Validate base seq does not work on DNA"
    assert validate_base_seq("AAUAGAU", True) == True, "Validate base seq does not work on RNA"
    assert validate_base_seq("Hi there!") == False, "Validate base seq fails to recognize nonDNA"
    assert validate_base_seq("Hi there!", True) == False, "Validate base seq fails to recognize nonDNA"
    print("Passed DNA and RNA tests")

def qual_score(phred_score):
    x = 0
    for ch in phred_score:
        x += (convert_phred(ch))
    return (x / len(phred_score))

if __name__ == "main":
    assert qual_score('GGGGTTTTTGGGGGGG') == 42.0625
    assert qual_score('GOOGOOGAGA') == 40.0
    assert qual_score('!!!!!!!!!!') == 0.0

def convert_phred(letter: str, ch=33):
    """Converts a single character into a phred score"""
    return ord(letter) - ch

if __name__ =="main":
    assert convert_phred("I") == 40, "wrong phred score for 'I'"
    assert convert_phred("C") == 34, "wrong phred score for 'C'"
    assert convert_phred("2") == 17, "wrong phred score for '2'"
    assert convert_phred("@") == 31, "wrong phred score for '@'"
    assert convert_phred("$") == 3, "wrong phred score for '$'"


def gc_content(DNA):
    '''Returns GC content of a DNA sequence as a decimal between 0 and 1.'''
    DNA = DNA.upper()         #Make sure sequence is all uppercase
    Gs = DNA.count("G")       #count the number of Gs
    Cs = DNA.count("C")       #count the number of Cs
    return (Gs+Cs)/len(DNA)

if __name__ == "main":
    assert gc_content("GCGCGC") == 1
    assert gc_content("AATTATA") == 0
    assert gc_content("GCATGCAT") == 0.5
    print("correctly calculated GC content")

if __name__ == "main":
    from os import path
    assert oneline_fasta('onelinetest.txt', 'funtimeout.txt') == (4, 4)
    oneline_fasta('onelinetest.txt', 'funtimeout.txt')
    assert path.exists('funtimeout.txt')

RNA = set('AUGCaucg')

DNA = set('ATCGNatcgn')