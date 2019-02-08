#Bryan Fisk
#2/7/19

import os, sys, re
import argparse
from Bio.Seq import Seq

def get_options():
    #get options from command line input
    parser = argparse.ArgumentParser(description = 'Generate reverse complement for sequence input.')
    parser.add_argument('--output', '-o', help = 'Set output file.')
    parser.add_argument('--input', '-i', help = 'Set input file.')
    parser.add_argument('--force', '-f', help = 'Force input of non-AGCT characters.', action = 'store_true')
    parser.add_argument('--custom', '-c', help = 'Custom reverse complement algorithm. Default uses Biopython.')
    args = parser.parse_args()
    return args

def reverse(s):
    """Return the sequence string in reverse order."""
    '''# make a list of letters from string
    letterlist = list(s)
    
    # reverse the list
    reversedlist = [letterlist[k] for k in range(len(letterlist) - 1, -1, -1)]
    
    # join the letters of the list into string and return
    return ''.join(reversedlist)'''
    return s[::-1]
    
def complement(s):
    """Return the complementary sequence string."""
    # dictionary setup for complement
    complement_dict = {'A' : 'T', 'T' : 'A', 'G' : 'C', 'C' : 'G'}
    
    # make a list of letters from string
    letterlist = list(s)
    
    # for loop of the letters and call the base_complementary dictionary
    complementlist = [complement_dict[letter] for letter in letterlist]
    
    # join the letters of the list into string and return
    return ''.join(complementlist)
    
def acgt(c):
    """Return True if character isn't any of ACGT or False otherwise"""
    if c not in list("ACGT"):
        return False
    else:
        return True
    
def get_input():
    # check DNA letter (only ACGTacgt)
    pattern = re.compile('^[ACGT]*$')
    input_validation = False
    while input_validation == False:
        # get input sequence
        if args.input:
            dna_seq = args.input
        else:
            dna_seq = input('Type your DNA sequence : ').upper() # change it to upper case
        if pattern.search(dna_seq):
            input_validation = True
        elif args.force:
            input_validation = True
        else:
            #if any(x in dna_seq for x in ['-f', '--force']):
            if any(x in dna_seq for x in ['-F', '--FORCE']):
                dna_seq = dna_seq.replace('-F', '')
                dna_seq = dna_seq.replace('--FORCE', '')
                #[dna_seq.replace(dna_seq, x, '') for x in [' -F', ' --FORCE']]
                return dna_seq
            print("Invalid input. If you would like to pass non ACGT characters, use -f, --force after input.")
            input_validation = False
    return dna_seq

def reverse_complement(dna_seq):
    # call reverse function
    # call complement function
    # print output
    if args.custom:
        return reverse(complement(dna_seq))
    else:
        seq = Seq(dna_seq)
        return seq.reverse_complement()
    # exit the program
    #exit()

def output(seq):
    if args.output:
        file = open(args.output, 'w')
        file.write(seq)
        file.close()
    else:
        print(seq)

if __name__ == '__main__':
    #[print(k) for k in dir(locals)]
    args = get_options()
    dna = get_input()
    print(dna)
    rc = reverse_complement(dna)
    output(rc)