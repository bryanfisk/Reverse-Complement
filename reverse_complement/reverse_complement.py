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
    parser.add_argument('--force', '-f', help = 'Force input of non-AGCT characters. (Uses Biopython)', action = 'store_true')
    parser.add_argument('--custom', '-c', help = 'Custom reverse complement algorithm. Default uses Biopython.', action = 'store_true')
    args = parser.parse_args()
    return args

def reverse(s):
    """Return the sequence string in reverse order."""
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

def get_file_input():
    with open(args.input, 'r') as file:
        inputstring = file.read()
        file.close()

    pattern = re.compile('[^ACGT]')
    inputstring = inputstring.split('\n')
    inputstring = [inputstring[x].strip() for x in range(len(inputstring))]
    inputstring = list(filter(None, inputstring))
    headers = []
    strings = []
    for x in inputstring:
        if x[0] == '>':
            headers.append(x)
            strings.append('')
        else:  
            strings[len(headers) - 1] += x
    if any([pattern.search(x) for x in strings]) and args.custom:
        f = input("Invalid input. Would you like to use non-ACGT characters? (Biopython)")
        if f in ['y', 'Y']:
            args.custom = False
        else:
            exit()
    return headers, strings

def get_input():
    # check DNA letter (only ACGTacgt)
    pattern = re.compile('^[ACGT]+$')
    input_validation = False
    while input_validation == False:
        # get input sequence
        dna_seq = input('Type your DNA sequence : ').upper() # change it to upper case
        if pattern.search(dna_seq) or args.force:
            input_validation = True
        elif dna_seq == None:
            input_validation = False
        else:
            f = input("Invalid input. Would you like to use non-ACGT characters? (Biopython)")
            if f in ['y', 'Y']:
                args.custom = False
                return dna_seq
            input_validation = False
    return dna_seq

def reverse_complement(dna_seq):
    # call reverse function
    # call complement function
    # print output
    if args.custom:
        comp = [''.join(list(map(complement, dna_seq)))]
        rev = ''.join(list(map(reverse, comp)))
        return rev
    else:
        seq = Seq(dna_seq)
        return seq.reverse_complement()
    # exit the program
    #exit()

def output(strings):
    if args.output:
        file = open(args.output, 'w')
        file.write(strings)
        file.close()
    else:
        print(strings)

def compile_output(strings, headers):
    output = ''
    for x in range(len(headers)):
        output += headers[x]
        output += '\n'
        for y in range(len(strings[x]) // 70 + 1):
            output += strings[x][0 + 70 * y:70 + 70 * y]
            output += '\n'
        output += '\n' 
    return output

if __name__ == '__main__':
    args = get_options()
    if args.input:
        headers, strings = get_file_input()
        strings = list(map(str, list(map(reverse_complement, strings))))
        rc = compile_output(strings, headers)
    else:
        dna = get_input()
        rc = str(reverse_complement(dna))
    output(rc)