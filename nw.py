#importing libraries
import math, sys, re

#chek if fasta file is provided in command line
if len(sys.argv) == 2:
    filename = sys.argv[1]
else:
    #if not, ask for user's input
    filename = input("Enter fasta file's name: ") #fasta_for_alignment.txt

#reading fasta file with 2 sequences    
def read_fasta(filename):
    seq1 = ""
    seq2 = ""
    count_seq = 0
    try:
        infile = open(filename, "r")
        for line in infile:
            if line.startswith(">"):
                count_seq += 1
            elif count_seq == 2:
                if re.search(r'^\D+$', line[:-1]):
                    seq2 += line[:-1]
            elif re.search(r'^\D+$', line[:-1]):
                seq1 += line[:-1]
        infile.close()
        
        if len(seq2) == 0 or len(seq1) == 0:
            raise ValueError("File is missing a sequence.")
            
        return seq1, seq2
    
    except IOError as io:
        print("Can't open file, reason:", str(io))
        sys.exit(1)
        
seq1, seq2 = read_fasta(filename)

#Implementing NW algorithm
def needleman_wunsch(seq1, seq2, match_score = 1, mismatch_score = -1, gap_penaly = -2):
    m = len(seq1)
    n = len(seq2)
    
    #zeros list of lists (matrix)
    matrix = [[0] * (m + 1) for _ in range(n+1)]
    
    #matrix 1st row and 1st column filled with gap penalties
    for i in range(len(matrix[0])):
        matrix[0][i] = i * gap_penaly
    for i in range(len(matrix)):
        matrix[i][0] = i * gap_penaly
    
    #craete matrix with scores
    for i in range(1, len(matrix)):
    #iterate through all the rows of the matrix starting from row 1, 0 row is skipped as it is filled with gap penalties.
        for j in range(1, len(matrix[0])):
            #check for number from the top
            top_number = matrix[i-1][j] + gap_penaly

            #check for number from the left
            left_number = matrix[i][j-1] + gap_penaly

            #check for number from the diagonal
            #row index = seq2 index, col index = seq1 index
            #check if letters match
            if seq2[i-1] == seq1[j-1]:
                diagonal_number = matrix[i-1][j-1] + match_score

            # letters don’t match
            else:
                diagonal_number = matrix[i-1][j-1] + mismatch_score

            max_number = max(top_number, left_number, diagonal_number)

            matrix[i][j] = max_number
    
    score = matrix[i][j]
    
    #tracing back alignment
    i, j = len(seq2), len(seq1)
    alignment_seq1 = ""
    alignment_seq2 = ""

    while i > 0 or j > 0:
        if matrix[i][j] == matrix[i-1][j-1] + (match_score if seq1[j-1] == seq2[i-1] else mismatch_score):
            alignment_seq1 = seq1[j-1] + alignment_seq1
            alignment_seq2 = seq2[i-1] + alignment_seq2
            i -= 1
            j -= 1
        elif matrix[i][j] == matrix[i][j-1] + gap_penaly:
            alignment_seq1 = seq1[j-1] + alignment_seq1
            alignment_seq2 = "_" + alignment_seq2
            j -= 1
        else:
            alignment_seq1 = "_" + alignment_seq1
            alignment_seq2 = seq2[i-1] + alignment_seq2
            i -= 1
            
    return alignment_seq1, alignment_seq2, score

alignment_seq1, alignment_seq2, score = needleman_wunsch(seq1, seq2, match_score = 1, mismatch_score = -1, gap_penaly = -2)

print("Alignment for 1st sequence:", alignment_seq1)
print("Alignment for 2nd sequence:", alignment_seq2)
print("Score:", score)
