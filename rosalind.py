# rosalind.info

# Counting DNA Nucleotides

# Problem when reading: because the text document in only a string, one apparently needs to double the backlashes in the file path to "escape the backlashes in the string"

# file_path = "C:\\Users\\saxel\\Downloads\\rosalind_dna.txt"
# with open(file_path, 'r') as file:
#     s = file.read()

# nt = [0, 0, 0, 0]
# for i in range(len(s)):
#     if s[i] == "A":
#         nt[0] += 1
#     elif s[i] == "C":
#         nt[1] += 1
#     elif s[i] == "G":
#         nt[2] += 1
#     elif s[i] == "T":
#         nt[3] += 1
# print(nt[0], " ", nt[1], " ", nt[2], " ", nt[3])


# Transcribing DNA into RNA

# file_path = "C:\\Users\\saxel\\Downloads\\rosalind_rna.txt"
# with open(file_path, 'r') as file:
#     s = file.read()

def dna_to_rna(s: 'str'):
    t = ""
    for i in range(len(s)):
        if s[i] != "T":
            t += s[i]
        else:
            t += "U"
    return t

# print(dna_to_rna(s))


# Complementing a Strand of DNA

# file_path = "C:\\Users\\saxel\\Downloads\\rosalind_revc (2).txt"
# with open(file_path, 'r') as file:
#     s = file.read()

def dna_revcomp(s: 'str'): # constructs the reverse complement of a gen
    sc = ""
    for c in reversed(s):
        if c == "A":
            sc += "T"
        elif c == "T":
            sc += "A"
        elif c == "C":
            sc += "G"
        elif c == "G":
            sc += "C"
    return sc

# print(dna_revcomp(s))


# Rabbits and Recurrence Relations

# n = 36
# k = 5
# i = 2
# f = [0] * n
# f[0] = 1
# f[1] = 1
# while i < n:
#     f[i] = f[i - 1] + k * f[i - 2]
#     i += 1
# print(f)


# Computing GC Content

# keys = []
# gena = []
# curr = ""
# started = False
# file_path = "C:\\Users\\saxel\\Downloads\\rosalind_gc (1).txt"
# with open(file_path, 'r') as file:
#     for line in file:
#         if not started:
#             keys.append(line[1:-1])
#             started = True
#         else:
#             if line[0] == ">":
#                 gena.append(curr)
#                 keys.append(line[1:-1])
#                 curr = ""
#             else:
#                 curr += line[:-1]
# gena.append(curr)

# cg = []
# max_cg = 0
# max_ind = 0
# for i in range(len(gena)):
#     counter = 0
#     for c in gena[i]:
#         if c in {"C", "G"}:
#             counter += 1
#     curr_cg = counter / len(gena[i])
#     cg.append(curr_cg)
#     if curr_cg > max_cg:
#         max_cg = curr_cg
#         max_ind = i
# print(keys, gena)
# print(keys[max_ind], "\n", 100 * cg[max_ind])



# Counting Point Mutations

# gena = []
# file_path = "C:\\Users\\saxel\\Downloads\\rosalind_hamm.txt"
# with open(file_path, 'r') as file:
#     for line in file:
#         gena.append(line)

# s = gena[0][:-1]
# t = gena[1][:-1]

# counter = 0
# for i in range(len(s)):
#     if s[i] != t[i]:
#         counter += 1
# print(counter)


# Mendel's First Law

# file_path = "C:\\Users\\saxel\\Downloads\\rosalind_iprb.txt"
# with open(file_path, 'r') as file:
#     s = file.read()

# s = s.split()
# for i in range(len(s)):
#     temp = 0
#     for c in s[i]:
#         temp = temp * 10
#         temp += ord(c) - ord("0")
#     s[i] = temp
# k = s[0]
# m = s[1]
# n = s[2]

def over(a: 'int', b: 'int') -> 'int': # computes the binomial number "a over b"
    if a < b:
        raise ValueError("Error: the first number has to be greater or equal than the second")
    curr = 1
    c = a - b
    while a > c:
        curr = curr * a
        a -= 1
    while b > 0:
        curr = curr // b
        b -= 1
    return curr

# p1 = over(k, 2) + over(k, 1) * over(m, 1) + over(k, 1) * over(n, 1)
# p2 = over(m, 2) * (3 / 4)
# p3 = over(m, 1) * over(n, 1) * (1 / 2)
# print((p1 + p2 + p3) / over(k + m + n, 2))



# Translating RNA into Protein

# file_path = "C:\\Users\\saxel\\Downloads\\rosalind_prot (2).txt"
# with open(file_path, 'r') as file:
#     s = file.read()

def rna_to_prot(s: "str"):
    n = len(s)
    i = 0
    prot = ""
    while i < n:
        codon = s[i: i + 3]
        if codon in {"UUU", "UUC"}:
            prot += "F"
        elif codon in {"UUA", "UUG", "CUU", "CUC", "CUA", "CUG"}:
            prot += "L"
        elif codon in {"UCU", "UCC", "UCA", "UCG", "AGU", "AGC"}:
            prot += "S"
        elif codon in {"UAU", "UAC"}:
            prot += "Y"
        elif codon in {"UGU", "UGC"}:
            prot += "C"
        elif codon in {"UGG"}:
            prot += "W"
        elif codon in {"CCU", "CCC", "CCA", "CCG"}:
            prot += "P"
        elif codon in {"CAU", "CAC"}:
            prot += "H"
        elif codon in {"CAA", "CAG"}:
            prot += "Q"
        elif codon in {"CGU", "CGC", "CGA", "CGG", "AGA", "AGG"}:
            prot += "R"
        elif codon in {"AUU", "AUC", "AUA"}:
            prot += "I"
        elif codon in {"AUG"}:
            prot += "M"
        elif codon in {"ACU", "ACC", "ACA", "ACG"}:
            prot += "T"
        elif codon in {"AAU", "AAC"}:
            prot += "N"
        elif codon in {"AAA", "AAG"}:
            prot += "K"
        elif codon in {"GUU", "GUC", "GUA", "GUG"}:
            prot += "V"
        elif codon in {"GCU", "GCC", "GCA", "GCG"}:
            prot += "A"
        elif codon in {"GAU", "GAC"}:
            prot += "D"
        elif codon in {"GAA", "GAG"}:
            prot += "E"
        elif codon in {"GGU", "GGC", "GGA", "GGG"}:
            prot += "G"
        elif codon in {"UAA", "UAG", "UGA"}:
            return prot
        i += 3
    return prot

# print(rna_to_prot(s))


# Finding a Motif in DNA

# file_path = "C:\\Users\\saxel\\Downloads\\rosalind_subs.txt"
# with open(file_path, 'r') as file:
#     r = file.read()

# r = r.split()
# s = r[0]
# t = r[1]

# n = len(s)
# m = len(t)
# i = 0
# repeats = []
# while i < n - m + 1:
#     if s[i:i + m] == t:
#         repeats.append(i + 1)
#     i += 1
# print(repeats)



# Consensus and Profile

def fasta_read(file_path: 'str'): # Reads gens with names in fasta format, with the caveat that the last string should end in a space or line jump, otherwise it does not detect the last base
    keys = []
    gena = []
    curr = ""
    started = False
    with open(file_path, 'r') as file:
        for line in file:
            if not started:
                keys.append(line[1:-1])
                started = True
            else:
                if line[0] == ">":
                    gena.append(curr)
                    keys.append(line[1:-1])
                    curr = ""
                else:
                    curr += line[:-1]
    gena.append(curr)
    return keys, gena

# file_path = "C:\\Users\\saxel\\Downloads\\rosalind_cons.txt"
# keys, gena = fasta_read(file_path)

# n = len(gena[0])
# profile_matrix = [[0] * n, [0] * n, [0] * n, [0] * n]
# consensus_str = ""
# for gen in gena:
#     for i in range(n):
#         if gen[i] == "A":
#             profile_matrix[0][i] += 1
#         elif gen[i] == "C":
#             profile_matrix[1][i] += 1
#         elif gen[i] == "G":
#             profile_matrix[2][i] += 1
#         elif gen[i] == "T":
#             profile_matrix[3][i] += 1

# for i in range(n):
#     col = [profile_matrix[j][i] for j in range(4)]
#     curr_max = 0
#     curr_ind = 0
#     k = 0
#     while k < 4:
#         if col[k] > curr_max:
#             curr_ind = k
#             curr_max = col[k]
#         k += 1
#     if curr_ind == 0:
#         consensus_str += "A"
#     if curr_ind == 1:
#         consensus_str += "C"
#     if curr_ind == 2:
#         consensus_str += "G"
#     if curr_ind == 3:
#         consensus_str += "T"

# print(consensus_str)
# print("A:", end="")
# for i in range(n):
#     print(" ", profile_matrix[0][i], end="")
# print("\n", end="")
# print("C:", end="")
# for i in range(n):
#     print(" ", profile_matrix[1][i], end="")
# print("\n", end="")
# print("G:", end="")
# for i in range(n):
#     print(" ", profile_matrix[2][i], end="")
# print("\n", end="")
# print("T:", end="")
# for i in range(n):
#     print(" ", profile_matrix[3][i], end="")



# Mortal Fibonacci Rabbits

def num_list_read(file_path: 'str'):
    with open(file_path, 'r') as file:
        s = file.read()
    s = s.split()
    for i in range(len(s)):
        temp = 0
        for c in s[i]:
            temp = temp * 10
            temp += ord(c) - ord("0")
        s[i] = temp
    return s


# file_path = "C:\\Users\\saxel\\Downloads\\rosalind_fibd.txt"
# s = num_list_read(file_path)

# n = s[0]
# m = s[1]

# adults = [0] * n
# children = [0] * n
# adults[0] = 0
# adults[1] = 1
# children[0] = 1
# children[1] = 0
# i = 2
# while i < n:
#     adults[i] += adults[i - 1] + children[i - 1]
#     if i >= m:
#         adults[i] -= children[i - m]
#     children[i] += adults[i - 1]
#     i += 1
# print(adults[n - 1] + children[n - 1])



# Overlap Graphs

# file_path = "C:\\Users\\saxel\\Downloads\\rosalind_grph.txt"
# keys, gena = fasta_read(file_path)

# for i in range(len(gena)):
#     suff = gena[i][-3:]
#     for j in range(len(gena)):
#         if i != j:
#             pref = gena[j][:3]
#             if suff == pref:
#                 print(keys[i], keys[j])



# Calculating Expected Offspring

# file_path = "C:\\Users\\saxel\\Downloads\\rosalind_iev.txt"
# s = num_list_read(file_path)

# p = s[0] + s[1] + s[2] + s[3] * (3 / 4) + s[4] * (1 / 2)
# print(2 * p)



# Finding a Shared Motif

# file_path = "C:\\Users\\saxel\\Downloads\\rosalind_lcsm (1).txt"
# keys, gena = fasta_read(file_path)

# # first find the matching substrings of the first two strings
# s1 = gena[0]
# s2 = gena[1]
# n = len(s1)
# m = len(s2)
# k = 1
# match = set()
# while k < n:
#     print(k)
#     found_k = False
#     for i in range(n - k + 1):
#         curr = s1[i:i + k]
#         j = 0
#         found = False
#         while found == False and j + k <= m:
#             if s2[j: j + k] == curr:
#                 found = True
#             j += 1
#         if found == True:
#             match.add(curr)
#             found_k = True
#     if found_k == False:
#         k = n
#     k += 1

# # then check if these matches appear also in all the rest
# l = 2
# g = len(gena)
# while l < g:
#     print(l)
#     match_list = list(match)
#     s3 = gena[l]
#     m = len(s3)
#     for i in range(len(match_list)):
#         curr = match_list[i]
#         k = len(curr)
#         j = 0
#         found = False
#         while found == False and j + k <= m:
#             if s3[j:j + k] == curr:
#                 found = True
#             j += 1
#         if found == False:
#             match.remove(curr)
#     l += 1

# longest = ""
# for ele in list(match):
#     if len(ele) > len(longest):
#         longest = ele
# print(longest)



# Independent Alleles

# file_path = "C:\\Users\\saxel\\Downloads\\rosalind_lia.txt"
# s = num_list_read(file_path)

# k = s[0]
# N = s[1]

# independently of your genes, if you mate with a Aa, you have 1/2 probability of getting Aa offspring
# this means, you have 1/4 probability of getting AaBb offspring by Mendel's second law
# and now the probability of having j occurrences of AaBb is the probability of a binomial with p = 1/4 and n = 2 ** k, the members of the generation:
# we go through each member and have 1/4 probability of that member being AaBb, and the probabilities of different members are independent
# moreover, we need to take into account all permutations of these members, thus the binomial coefficient counts these different orders

# p = 0
# n = 2 ** k
# i = N
# while i <= n:
#     p += over(n, i) * ((1 / 4) ** i) * ((3 / 4) ** (n - i))
#     i += 1
# print(p)



# Finding a Protein Motif

from urllib.request import urlopen

def fasta_uniprot_read(file_path: 'str'):
    with open(file_path, 'r') as file:
        s = file.read()
    s = s.split()

    prots = [""] * len(s)

    for i in range(len(s)):
        temp = ""
        j = 0
        while j < len(s[i]) and s[i][j] != "_":
            temp += s[i][j]
            j += 1
        prots[i] = temp
    seqs = []
    curr = ""
    for i in range(len(prots)):
        r = urlopen("http://www.uniprot.org/uniprot/" + prots[i] + ".fasta").read()
        r = r.decode("utf-8")
        i = 0
        while r[i:i + 3] != "SV=":
            i += 1
        i += 4
        curr = r[i:]
        curr = "".join(curr.split())
        seqs.append(curr)
    return s, prots, seqs

# file_path = "C:\\Users\\saxel\\Downloads\\rosalind_mprt (2).txt"
# s, prots, seqs = fasta_uniprot_read(file_path)

# Nglyc = [] * len(s)
# for i in range(len(s)):
#     j = 0
#     positions = []
#     while j < len(seqs[i]) - 3:
#         if seqs[i][j] == "N":
#             if seqs[i][j + 1] != "P" and seqs[i][j + 2] in {"S", "T"} and seqs[i][j + 3] != "P":
#                 positions.append(j + 1) # remember they count positions with 1-indexing
#         j += 1 # NNTST would give two N-glycosylation motifs: NNTS and NTST !
#     Nglyc.append(positions)
# for i in range(len(s)):
#     if Nglyc[i]:
#         print(s[i])
#         for num in Nglyc[i]:
#             print(num, end=" ")
#         print("")



# Inferring mRNA from Protein

# file_path = "C:\\Users\\saxel\\Downloads\\rosalind_mrna (1).txt"
# with open(file_path, 'r') as file:
#     s = file.read()
# combos = 3
# for c in s:
#     if c in {"F", "Y", "C", "H", "Q", "K", "N", "D", "E"}:
#         combos = combos * 2
#     elif c == "I":
#         combos = combos * 3
#     elif c in {"P", "T", "V", "A", "G"}:
#         combos = combos * 4
#     elif c in {"L", "S", "R"}:
#         combos = combos * 6
#     while combos > 1000000:
#         combos -= 1000000
# print(combos)



# Open Reading Frames

# file_path = "C:\\Users\\saxel\\Downloads\\rosalind_orf (2).txt"
# keys, gena = fasta_read(file_path)
# gen = gena[0]

def dna_to_prot(gen: 'str'):
    n = len(gen)
    i = 0
    prot = ""
    while i < n:
        codon = gen[i: i + 3]
        if codon in {"TTT", "TTC"}:
            prot += "F"
        elif codon in {"TTA", "TTG", "CTT", "CTC", "CTA", "CTG"}:
            prot += "L"
        elif codon in {"TCT", "TCC", "TCA", "TCG", "AGT", "AGC"}:
            prot += "S"
        elif codon in {"TAT", "TAC"}:
            prot += "Y"
        elif codon in {"TGT", "TGC"}:
            prot += "C"
        elif codon in {"TGG"}:
            prot += "W"
        elif codon in {"CCT", "CCC", "CCA", "CCG"}:
            prot += "P"
        elif codon in {"CAT", "CAC"}:
            prot += "H"
        elif codon in {"CAA", "CAG"}:
            prot += "Q"
        elif codon in {"CGT", "CGC", "CGA", "CGG", "AGA", "AGG"}:
            prot += "R"
        elif codon in {"ATT", "ATC", "ATA"}:
            prot += "I"
        elif codon in {"ATG"}:
            prot += "M"
        elif codon in {"ACT", "ACC", "ACA", "ACG"}:
            prot += "T"
        elif codon in {"AAT", "AAC"}:
            prot += "N"
        elif codon in {"AAA", "AAG"}:
            prot += "K"
        elif codon in {"GTT", "GTC", "GTA", "GTG"}:
            prot += "V"
        elif codon in {"GCT", "GCC", "GCA", "GCG"}:
            prot += "A"
        elif codon in {"GAT", "GAC"}:
            prot += "D"
        elif codon in {"GAA", "GAG"}:
            prot += "E"
        elif codon in {"GGT", "GGC", "GGA", "GGG"}:
            prot += "G"
        elif codon in {"TAA", "TAG", "TGA"}:
            return prot
        i += 3
    if i == n:
        return "not a protein"              

# i = 0
# g = len(gen)
# prots = set()
# while i < g - 2:
#     if gen[i: i + 3] == "ATG":
#         prot = dna_to_prot(gen[i:])
#         if prot != "not a protein":
#             prots.add(prot)
#         else:
#             i = g - 3
#     i += 1

# revcom_gen = ""
# for i in range(g):
#     if gen[i] == "C":
#         revcom_gen += "G"
#     elif gen[i] == "G":
#         revcom_gen += "C"
#     elif gen[i] == "A":
#         revcom_gen += "T"
#     elif gen[i] == "T":
#         revcom_gen += "A"
# revcom_gen = revcom_gen[::-1]

# i = 0    
# while i < g - 2:
#     if revcom_gen[i: i + 3] == "ATG":
#         prot = dna_to_prot(revcom_gen[i:])
#         if prot != "not a protein":
#             prots.add(prot)
#         else:
#             i = g - 3
#     i += 1

# for prot in list(prots):
#     if prot: # avoid getting Nones
#         print(prot)



# Enumerating Gene Orders

import itertools

# file_path = "C:\\Users\\saxel\\Downloads\\rosalind_perm (1).txt"
# s = num_list_read(file_path)

# with open("output.txt", "a") as f:
#     n = s[0]

#     temp = 1
#     i = 2
#     while i <= n:
#         temp = temp * i
#         i += 1
#     print(temp, file=f)

#     numbers = set()
#     for i in range(n):
#         numbers.add(i + 1)

#     for subset in list(itertools.permutations(numbers)):
#         subset = list(subset)
#         for i in range(len(subset)):
#             if i != len(subset) - 1:
#                 print(subset[i], end=" ", file=f)
#             else:
#                 print(subset[i], file=f)



# Calculating Protein Mass

# file_path = "C:\\Users\\saxel\\Downloads\\rosalind_prtm.txt"
# with open(file_path, 'r') as file:
#     s = file.read()

# weight = 0
# for c in s:
#     if c == "A":
#         weight += 71.03711
#     elif c == "C":
#         weight += 103.00919
#     elif c == "D":
#         weight += 115.02694
#     elif c == "E":
#         weight += 129.04259
#     elif c == "F":
#         weight += 147.06841
#     elif c == "G":
#         weight += 57.02146
#     elif c == "H":
#         weight += 137.05891
#     elif c == "I":
#         weight += 113.08406
#     elif c == "K":
#         weight += 128.09496
#     elif c == "L":
#         weight += 113.08406
#     elif c == "M":
#         weight += 131.04049   
#     elif c == "N":
#         weight += 114.04293
#     elif c == "P":
#         weight += 97.05276
#     elif c == "Q":
#         weight += 128.05858
#     elif c == "R":
#         weight += 156.10111
#     elif c == "S":
#         weight += 87.03203
#     elif c == "T":
#         weight += 101.04768
#     elif c == "V":
#         weight += 99.06841
#     elif c == "W":
#         weight += 186.07931
#     elif c == "Y":
#         weight += 163.06333 
# print(weight)



# Locating Restriction Sites

# file_path = "C:\\Users\\saxel\\Downloads\\rosalind_revp (1).txt"
# keys, gena = fasta_read(file_path)

# gen = gena[0]
# n = len(gen)

# with open("output.txt", "a") as f:
#     conj = {"AT", "TA", "CG", "GC"}
#     i = 0
#     while i < n - 1:
#         if gen[i:i + 2] in conj:
#             j = i - 1
#             k = i + 2
#             l = 2
#             while j >= 0 and k < n and gen[j] + gen[k] in conj:
#                 j -= 1
#                 k += 1
#                 l += 2
#             while l >= 4:
#                 # print(gen[j + 1: j + 1 + l], file=f)
#                 print(j + 2, l, file=f)
#                 l -= 2
#                 j += 1
#         i += 1



# RNA Splicing

# file_path = "C:\\Users\\saxel\\Downloads\\rosalind_splc.txt"
# keys, gena = fasta_read(file_path)

# gen = gena[0]
# n = len(gen)
# introns = []
# introns_len = []
# for i in range(len(gena) - 1):
#     introns.append(gena[i + 1])
#     introns_len.append(len(gena[i + 1]))

# curated= ""
# i = 0
# while i < n:
#     cut = False
#     for j in range(len(introns)):
#         if i + introns_len[j] <= n and gen[i:i + introns_len[j]] == introns[j]:
#             i += introns_len[j]
#             cut = True
#     if cut == False:
#         curated += gen[i]
#         i += 1
#     cut = False

# curated = dna_to_rna(curated)
# curated = rna_to_prot(curated)
# print(curated)


# Enumerating k-mers Lexicographically

# file_path = "C:\\Users\\saxel\\Downloads\\rosalind_lexf.txt"
# with open(file_path, 'r') as file:
#     s = file.read()

# s = s.split()
# n = len(s)
# num = 0
# temp = s[n - 1]
# for c in temp:
#     num = num * 10
#     num += ord(c) - ord("0")
# letters = []
# for i in range(n - 1):
#     letters.append(s[i])
# # letters.sort() # this option is only if the order is not given but derived from the usual alphabet

# for option in itertools.product(letters, repeat = num):
#     for let in option:
#         print(let, end="")
#     print("")   



# Longest Increasing Subsequence

# file_path = "C:\\Users\\saxel\\Downloads\\rosalind_lgis (1).txt"
# s = num_list_read(file_path)
# n = s[0]

# nums = [s[i + 1] for i in range(n)]

# # For ascending

# # We construct a new array of the same length
# # and store in each position the longest increasing
# # subsequence until that point: for that we look
# # at previous numbers that are smaller than the
# # current one and among them select the one with
# # the longest increasing subsequence ending in it

# best = 0
# k = 0

# aux = [1 for _ in range(n)]
# for i in range(n):
#     curr = nums[i]
#     for j in range(i):
#         if nums[i] > nums[j] and aux[i] <= aux[j]:
#             aux[i] = aux[j] + 1
#         if aux[i] > best:
#             k = i
#             best = aux[i]

# # The longest increasing subsequence will be the one
# # which corresponding greatest number in aux;
# # from there on, we traverse the list bacwards by
# # looking at the numbers whose aux is always 1 lower.
# # This will give us the ascending sequence in reverse
# # order, so we also need to print it from the end back

# i = k
# ascen = [nums[i]]
# curr = best
# while curr > 1:
#     if aux[i] == curr - 1:
#         ascen.append(nums[i])
#         curr = curr - 1
#     i -= 1

# m = len(ascen)
# for i in range(m - 1):
#     print(ascen[m - i - 1], end=" ")
# print(ascen[0])

# # For descending

# best = 0
# k = 0

# aux = [1 for _ in range(n)]
# for i in range(n):
#     curr = nums[i]
#     for j in range(i):
#         if nums[i] < nums[j] and aux[i] <= aux[j]: # only change wrt ascending is the sign in the first condition here
#             aux[i] = aux[j] + 1
#         if aux[i] > best:
#             k = i
#             best = aux[i]


# i = k
# descen = [nums[i]]
# curr = best
# while curr > 1:
#     if aux[i] == curr - 1:
#         descen.append(nums[i])
#         curr = curr - 1
#     i -= 1

# m = len(descen)
# for i in range(m - 1):
#     print(descen[m - i - 1], end=" ")
# print(descen[0])



# Genome Assembly as Shortest Superstring

# file_path = "C:\\Users\\saxel\\Downloads\\rosalind_long.txt"
# keys, gena = fasta_read(file_path)

def matching_dna(gen1, gen2):
    # tells you if the end of gen1 and the start of gen2 overlap, where the overlaps are more than half the length of the genes
    # if they overlap, it tells you from which position on in the first they overlap
    # otherwise, it returns -1
    n = len(gen1)
    m = len(gen2)
    k = max(m, n)
    k += k % 2
    k = k // 2
    i = n - k - 1 # -1 because we say "more" than half of their length overlap
    found = False
    while i >= 0 and found == False:
        j1 = i
        j2 = 0
        while j1 < n and j2 < m and gen1[j1] == gen2[j2]:
            j1 += 1
            j2 += 1
        if j1 == n:
            found = True
        else:
            i -= 1
    if found == True:
        return i
    else:
        return -1 

# with open("output.txt", "a") as f:
#     # for each string, look at the one it should be pasted to (j) and at which position (a)
#     sorting = [[-1, -1]] * len(gena)
#     used = set()
#     for i in range(len(gena)):
#         gen1 = gena[i]
#         j = 0
#         while j < len(gena):
#             if i != j and j not in used:
#                 a = matching_dna(gen1, gena[j])
#                 if a != -1:
#                     sorting[i] = [j, a]
#                     used.add(j)
#                     j = len(gena)
#             j += 1

#     # let us look for the start of the dna sequence
#     i = 0
#     while i in used:
#         i += 1

#     # now patch them together
#     output = gena[i][:sorting[i][1]]
#     i = sorting[i][0]
#     counter = 1
#     while counter < len(gena) - 1:
#         output += gena[i][:sorting[i][1]]
#         i = sorting[i][0]
#         counter += 1
#     output += gena[i]
#     print(output, file=f)




# Perfect Matchings and RNA Secondary Structures

# file_path = "C:\\Users\\saxel\\Downloads\\rosalind_pmch.txt"
# keys, gena = fasta_read(file_path)

# # Since the problem allows for adjacent edges to be matched, as seen in Figure 4,
# # we just need to count the appearences of A and U and of C and G, and the number
# # of perfect matchings will be (#A)! * (#C)!, assuming A and U appear the same number
# # of times, and C and G too

# gen = gena[0]
# AU = 0
# CG = 0
# for let in gen:
#     if let == "A":
#         AU += 1
#     elif let == "C":
#         CG += 1
# matchings = 1
# while AU > 0:
#     matchings = matchings * AU
#     AU -= 1
# while CG > 0:
#     matchings = matchings * CG
#     CG -= 1
# print(matchings)



# Partial Permutations

# # The total number of permutations of k-subsets inside n-sets can be counted by
# # first choosing the k elements in the subset from the set and then computing their 
# # permutations. Thus, the formula should be P(n,k) = (n over k) * k!

def factorial(n: 'int') -> 'int': # Computes the factorial of a number
    if n <= 0:
        raise ValueError("Error: the number has to be strictly greater than 0")
    temp = 1
    while n > 0:
        temp = temp * n
        n -= 1
    return temp


# file_path = "C:\\Users\\saxel\\Downloads\\rosalind_pper.txt"
# s = num_list_read(file_path)

# n = s[0]
# k = s[1]

# print(over(n, k) * factorial(k) % (10 ** 6))   



# Introduction to Random Strings

# file_path = "C:\\Users\\saxel\\Downloads\\rosalind_prob.txt"
# with open(file_path, 'r') as file:
#     s = file.read()
# s = s.split()
# for i in range(len(s) - 1):
#     temp = 0
#     j = 2
#     while j < len(s[i + 1]):
#         temp = temp * 10
#         temp += ord(s[i + 1][j]) - ord("0")
#         j += 1
#     s[i + 1] = temp / (10 ** (len(s[i + 1]) - 2))
# gen = s.pop(0)

# # To compute the probability of the random string with GC-content x being exactly
# # a particular one, we first need to count the AT and CG appearances. The probability
# # of obtaining a AT in a position is then (1 - x) / 2, and that of a CG is x / 2.
# # Thus, the probability of matching the given string is (((1 - x) / 2) ** AT) * ((x / 2) ** CG)
# # Finally, we take the common logarithm of that probability for the exercise: to compute
# # the common logarithm only using the natural logarithm we notice that log_10 x = y <==>
# # x = 10 ** y <==> ln x = y * ln 10 <==> y = (ln x) / (ln 10)

# from math import log

# AT = 0
# CG = 0
# for let in gen:
#     if let in {"A", "T"}:
#         AT += 1
#     elif let in {"C", "G"}:
#         CG += 1

def round3dec(x): # Rounds a decimal number to 3 digits
    x_str = str(x)
    if "." not in x_str: # If it has no decimals digits, it is already rounded
        return x
    temp = 0
    i = 0
    negative = False
    if x_str[i] == "-":
        negative = True
        i += 1
    while x_str[i] != ".":
        temp = temp * 10
        temp += ord(x_str[i]) - ord("0")
        i += 1
    if len(x_str[i:]) <= 4: # if it has exactly 3 decimal digits, there is no rounding necessary
        return x
    i += 1
    for _ in range(3):
        temp = temp * 10
        temp += ord(x_str[i]) - ord("0")
        i += 1
    check = ord(x_str[i]) - ord("0")
    if check >= 5:
        temp += 1
    if negative == True:
        temp = - temp
    return temp / 1000   

# prob = [log((((1 - x) / 2) ** AT) * ((x / 2) ** CG))/(log(10))  for x in s]
# for num in prob:
#     print(round3dec(num), end=" ")    



# Enumerating Oriented Gene Orderings

# file_path = "C:\\Users\\saxel\\Downloads\\rosalind_sign.txt"
# s = num_list_read(file_path)

# n = s[0]

# # To count all signed permutations, we first need to choose which of the integers are negative
# # This corresponds to choosing a subset of the set {1, ..., n}, inclusive the empty set, which
# # corresponds to none being negative. There are 2 ** n subsets in an n-set.
# # For each of these, we just count the possible permutations of n elements, which are n!


# # To generate the possible combinations, we use the numbers from 0 to 2 ** (n - 1) in binary:
# # they can be written as a list of n 0's and 1's, and the 0's will mean we take the number as
# # positive and 1 as negative. Then, we construct their permutations.

# nums = [i + 1 for i in range(n)]

def dec2bin(x, n): # turns an integer number in decimal representation into a binary list of length n corresponding to each binary representation in inverse order (works only for integers 0 <= x < 2 ** n)
    xbin = [0 for _ in range(n)]
    i = 0
    while x > 0:
        xbin[i] = x % 2
        x = (x - x % 2) // 2
        i += 1
    return xbin

# with open("output.txt", "a") as f:
#     print((2 ** n) * factorial(n), file=f)
#     for j in range(2 ** n):
#         jbin = dec2bin(j, n)
#         for i in range(n):
#             if jbin[i] == 1:
#                 nums[i] = -abs(nums[i])
#             elif jbin[i] == 0:
#                 nums[i] = abs(nums[i])
#         for subset in itertools.permutations(nums):
#             subset = list(subset)
#             for i in range(len(subset)):
#                 if i != len(subset) - 1:
#                     print(subset[i], end=" ", file=f)
#                 else:
#                     print(subset[i], file=f)



# Finding a Spliced Motif

# file_path = "C:\\Users\\saxel\\Downloads\\rosalind_sseq.txt"
# keys, gena = fasta_read(file_path)

# s = gena[0]
# t = gena[1]
# k = len(t)

# i = 0
# j = 0
# pos = [0 for _ in range(k)]
# while j < len(t):
#     if s[i] == t[j]:
#         pos[j] = i + 1
#         j += 1
#     i += 1
# for i in range(k):
#     if i != k - 1:
#         print(pos[i], end=" ")
#     else:
#         print(pos[i])




# Transitions and Transversions

# file_path = "C:\\Users\\saxel\\Downloads\\rosalind_tran.txt"
# keys, gena = fasta_read(file_path)

# s1 = gena[0]
# s2 = gena[1]
# k = len(s1)

# trs = 0 # Transitions
# trv = 0 # Transversions

# for i in range(k):
#     if (s1[i] == "A" and s2[i] == "G") or (s1[i] == "G" and s2[i] == "A") or (s1[i] == "C" and s2[i] == "T") or (s1[i] == "T" and s2[i] == "C"):
#         trs += 1
#     elif (s1[i] in {"A", "G"} and s2[i] in {"C", "T"}) or (s2[i] in {"A", "G"} and s1[i] in {"C", "T"}):
#         trv += 1
# print(trs/trv)




# Completing a Tree

# # A graph on n nodes is a tree iff it has n - 1 edges and no cycles

# file_path = "C:\\Users\\saxel\\Downloads\\rosalind_tree.txt"
# s = num_list_read(file_path)

# n = s[0]
# k = (len(s) - 1) // 2
# print(n - 1 - k)



# Catalan Numbers and RNA Secondary Structures

# file_path = "C:\\Users\\saxel\\Downloads\\rosalind_cat (3).txt"
# keys, gena = fasta_read(file_path)

# gen = gena[0]

def catalan(n): # Constructs the n-th Catalan number
    # The following will be the dynamic programming solution, but we can derive a closed formula from the recursion and this will be easier to compute
    # if n == 0:
    #     return 1
    # temp = 0
    # i = 1
    # while i <= n:
    #     temp += (catalan(i - 1) * catalan(n - i))
    #     i += 1
    # return temp

    # Using the closed formula:
    if n == 0:
        return 1
    return factorial(2 * n) // (factorial(n) * factorial(n + 1))

# # # Conjecture: if several A's or U's appear in a row, we can change them to just one
# # # As a result, we obtain an alternating sequence of A's and U's, so bases in odd positions
# # # can be paired with bases in even positions, i.e. we can apply Catalan numbers to that
# # # This means we just need to count the changes from A to U in the string (analogous for C and G)

# # AU_swaps = 0
# # CG_swaps = 0
# # AU_last = ""
# # CG_last = ""
# # for c in gen:
# #     if c == "A":
# #         if AU_last == "":
# #             AU_last = "A"
# #         elif AU_last == "U":
# #             AU_swaps += 1
# #             AU_last = "A"
# #     elif c == "U":
# #         if AU_last == "":
# #             AU_last = "U"
# #         elif AU_last == "A":
# #             AU_swaps += 1
# #             AU_last = "U"
# #     elif c == "C":
# #         if CG_last == "":
# #             CG_last = "C"
# #         elif CG_last == "G":
# #             CG_swaps += 1
# #             CG_last = "C"
# #     elif c == "G":
# #         if CG_last == "":
# #             CG_last = "G"
# #         elif CG_last == "C":
# #             CG_swaps += 1
# #             CG_last = "G"

# # # If the number of swaps is even, this means we started and ended with the same type of base:
# # # if we were to reorder the string to avoid that, we would have gotten one swap less, so we get rid of this extra one
# # # However, we want to count one extra node, the one we start with, so in practice we just add 1 to the argument of
# # # catalan numbers, i.e. to alternating knots, only if the number of swaps is odd
# # AU_cat = catalan((AU_swaps + AU_swaps % 2) // 2) % (10 ** 6)
# # CG_cat = catalan((CG_swaps + CG_swaps % 2) // 2) % (10 ** 6)
# # print((AU_cat * CG_cat) % (10 ** 6))
# # # This solution allows for crossings between AU and CG pairings, which is also not wanted!  

# # The solution without AU-CG crossings is the following:
    




# Error Correction in Reads

# file_path = "C:\\Users\\saxel\\Downloads\\rosalind_corr.txt"
# keys, gena = fasta_read(file_path)

def Hamming_dist(s: str, t: str): # computes the Hamming distance between gen-strings of equal length
    if len(s) != len(t):
        raise ValueError("Error: the gens need to have the same length as strings")
    dist = 0
    for i in range(len(s)):
        if s[i] != t[i]:
            dist += 1
    return dist

# # We first sort the genes between those repeated and those not repeated: the not repeated ones are the defect reads
# correct_gena = []
# incorrect_gena = []
# while len(gena) >= 2:
#     curr_gen = gena.pop(0)
#     curr_gen_revcomp = dna_revcomp(curr_gen) # we consider equivalence modulo reverse complement
#     matching = 0 # keeps count of whether we find the same gene read again
#     i = 0
#     while i < len(gena):
#         if gena[i] == curr_gen or gena[i] == curr_gen_revcomp:
#             gena.pop(i)
#             matching += 1
#             i -= 1
#         i += 1
#     if matching > 0:
#         correct_gena.append(curr_gen)
#     else:
#         incorrect_gena.append(curr_gen)
# if len(gena) == 1: # if only one is left, there is non to compare it to
#     incorrect_gena.append(gena.pop(0))

# # We compare the defect reads with the correct ones to find which among them has Hamming distance one to the defect
# with open("output.txt", "a") as f:
#     for gen in incorrect_gena:
#         j = 0
#         while j < len(correct_gena):
#             if Hamming_dist(gen, correct_gena[j]) == 1:
#                 print(gen + "->" + correct_gena[j], file=f)
#                 j = len(correct_gena)
#             elif Hamming_dist(gen, dna_revcomp(correct_gena[j])) == 1:
#                 print(gen + "->" + dna_revcomp(correct_gena[j]), file=f)
#                 j = len(correct_gena)
#             j += 1
    




# Counting Phylogenetic Ancestors

# # This is more of combinatorics-on-trees problem rather than a coding exercise
# # We use that the sum of degrees of nodes on the tree has to be equal to two times the number of edges,
# # since each edge is adjacent to exactly two nodes. Moreover, in trees, the number of edges is equal to
# # the number of nodes minus 1. Moreover, we can distinguish between internal nodes, which have degree 3,
# # and external nodes (leaves), which have degree 1 (there are no other options in an unrooted binary tree,
# # as defined). This gives the formula: internal nodes = leaves - 2

# file_path = "C:\\Users\\saxel\\Downloads\\rosalind_inod.txt"
# s = num_list_read(file_path)

# n = s[0]

# print(n - 2)





# k-Mer Composition

# file_path = "C:\\Users\\saxel\\Downloads\\rosalind_kmer.txt"
# keys, gena = fasta_read(file_path)

# gen = gena[0]

# # We first construct all 4-mers in lexicographical order by using the solution to the question "Enumerating k-mers Lexicographically"

# letters = ["A", "C", "G", "T"]
# num = 4

# comp = [0 for _ in range(len(letters) ** num)] # this range describes the possible num-mers with the alphabet in letters
# i = 0
# for option in itertools.product(letters, repeat = num):
#     curr = ""
#     for let in option:
#         curr += let # curr represents each of the 4-mers, and they are generated in lexicographical order
#     # Then, we check for each of the 4-mers when it appears in the gene, and store the number of occurrences in comp
#     j = 0
#     temp = 0
#     while j < len(gen) - num + 1:
#         if gen[j:j + num] == curr:
#             temp += 1
#         j += 1
#     comp[i] = temp
#     print(temp,end=" ")




# Speeding Up Motif Finding

# # The Knuth-Morris-Pratt algorithm is mentioned here, would be nice to implement it at some point!

# file_path = "C:\\Users\\saxel\\Downloads\\rosalind_kmp.txt"
# keys, gena = fasta_read(file_path)

# gen = gena[0]
# n = len(gen)

# failure = [0 for _ in range(n)]
# left = 1
# print(0, end=" ") # this represents failure[0]
# while left < n:
#     span = 1
#     temp = 0
#     while gen[left:left + span] == gen[:span]: # look for the longest substring starting at left which equals a prefix of gen
#         span += 1
#     for i in range(span - 1): # this longest substring indicates that at left + i there is a string of length i + 1 corresponding to a prefix, given i < span - 1; this might or might not be the optimal one for failure[left + i]
#         failure[left + i] = max(failure[left + i], i + 1)
#     left += 1
# with open("output.txt", "a") as f:
#     for fail in failure:
#         print(fail, end=" ", file=f)



# Finding a Shared Spliced Motif

# # This is the usual problem of finding the longest common subsequence of two strings

# We readjust the recursion limit so that the computer does not complain about making too many recursion steps
import sys
sys.setrecursionlimit(10 ** 4)

def LCS(s: str, t: str, i: int, j: int, memo: list[list[int]]): # finds the longest common subsequence between s[:i] and t[:j] by recursion aka DP
    if i > len(s) or j > len(t):
        raise Exception("One of the indices is larger than the size of the corresponding string")
    if i == 0 or j == 0:
        return ""
    if memo[i][j] != "":
        return memo[i][j]
    if s[i - 1] == t[j - 1]:
        memo[i][j] = LCS(s, t, i - 1, j - 1, memo) + s[i - 1]
        return memo[i][j]
    else:
        a = LCS(s, t, i - 1, j, memo)
        b = LCS(s, t, i, j - 1, memo)
        if len(a) >= len(b):
            memo[i][j] = a
            return memo[i][j]
        else:
            memo[i][j] = b
            return memo[i][j]
            

# file_path = "C:\\Users\\saxel\\Downloads\\rosalind_lcsq (1).txt"
# keys, gena = fasta_read(file_path)

# s = gena[0]
# t = gena[1]

# m = len(s)
# n = len(t)

# memo = [["" for _ in range(n + 1)] for _ in range(m + 1)]

# with open("output.txt", "a") as f:
#     print(LCS(s, t, m, n, memo), file=f)




# Ordering Strings of Varying Length Lexicographically

# file_path = "C:\\Users\\saxel\\Downloads\\rosalind_lexv.txt"
# with open(file_path, 'r') as file:
#     s = file.read()

# s = s.split()
# n = len(s)
# num = 0
# temp = s[n - 1]
# for c in temp:
#     num = num * 10
#     num += ord(c) - ord("0")
# letters = []
# for i in range(n - 1):
#     letters.append(s[i])

# done = set() # keep track of words already printed
# with open("output.txt", "a") as f:
#     for option in itertools.product(letters, repeat = num):
#         curr = ""
#         for let in option:
#             curr += let
#         # we just need to print the prefixes of each word we would have done in the case of strings of fixed length equals the maximum length allowed here, without repetition
#         for i in range(len(curr)):
#             temp = curr[:i + 1]
#             if temp not in done:
#                 done.add(temp)
#                 print(temp, file=f)




# Maximum Matchings and RNA Secondary Structures

# file_path = "C:\\Users\\saxel\\Downloads\\rosalind_mmch.txt"
# keys, gena = fasta_read(file_path)

# s = gena[0]

# A = 0
# C = 0
# G = 0
# U = 0

# for c in s:
#     if c == "A":
#         A += 1
#     elif c == "C":
#         C += 1
#     elif c == "G":
#         G += 1
#     elif c == "U":
#         U += 1

# # For the matchings between A and U, we can count the number of each of the bases and then if
# # wlog #A <= #U, we just need to choose how to much each of the A's with U's. For the first A
# # we have #U possibilities; for the next, #U - 1, and so on. Same with C and G

# if (A == 0 or U == 0) and (C == 0 or G == 0): # there might be no matching possible
#     print(0)
# else:
#     temp = 1
#     n = min(A, U)
#     m = max(A, U)
#     k = min(C, G)
#     l = max(C, G)
#     while n > 0:
#         temp = temp * m
#         m -= 1
#         n -= 1
#     while k > 0:
#         temp = temp * l
#         l -= 1
#         k -= 1
#     print(temp)




# Creating a Distance Matrix

# # We use our Hamming distance function and just compute things brute force

# file_path = "C:\\Users\\saxel\\Downloads\\rosalind_pdst (1).txt"
# keys, gena = fasta_read(file_path)

# n = len(gena)
# k = len(gena[0]) # we assume all genes have the same length

# dist_mat = [[0.0 for _ in range(n)] for _ in range(n)]

# with open("output.txt", "a") as f:
#     for i in range(n):
#         for j in range(n):
#             if i == j:
#                 print("{:.3f}".format(0.0), end=" ", file=f) # we use format to present exactly 3 digits after the coma
#             else:
#                 dist_mat[i][j] = Hamming_dist(gena[i], gena[j]) / k
#                 print("{:.3f}".format(dist_mat[i][j]), end=" ", file=f)
#         print("", file=f)





# Reversal Distance

# # We first define some relevant functions and extract the data from the file

def perm_inv(p: list[int]): # computes the inverse of a permutation
    n = len(p)
    q = [0 for _ in range(len(p))]
    for i in range(len(p)):
        q[p[i] - 1] = i + 1
    return q

def perm_comp(p: list[int], q: list[int]): # computes the composition of two permutations, where p is applied first and q second
    n = len(p)
    if n != len(q):
        raise ValueError("The permutations need to have the same length")
    return [q[p[i] - 1] for i in range(n)]

# file_path = "C:\\Users\\saxel\\Downloads\\rosalind_rear.txt"
# s = num_list_read(file_path)

# n = len(s)
# k = n // 10

# ps = [] # collection of the first permutations in each pair
# qs = [] # collection of the second permutations in each pair

# for i in range(k // 2):
#     ps.append(s[20 * i: 20 * i + 10])
#     qs.append(s[20 * i + 10: 20 * i + 20])

def perm_rev(p: list[int], i: int, j: int): # performs a reversal of the permutation p between the positions i and j (0-indexed and inclusive)
    q = [p[i] for i in range(len(p))]
    k = 0
    while i + k <= j:
        q[i + k] = p[j - k]
        k += 1
    return q

# # I could follow the literature on this problem, which finds algorithms in polynomial time.
# # However, I would first like to come up with my own solution: it is going to be brute force for now
# # I start by noticing that the reversal distance between p and q equals the reversal distance
# # between (p ** -1) * q and the identity, so we just need to compute the reversal distance from a permutation to the id

# # I decided to work with strings because lists are not hashable and start by 0 instead of 1 to only have 1-digit numbers

# p0 = [i + 1 for i in range(10)] # the identity permutation

def perm2str(p: list[int]): # takes a permutation as a list of integers between 1 and 10 and outputs the string (where each number is lowered by 1 so as to have only 1-digit numbers in the string)
    s = ""
    for num in p:
        s += str(num - 1)
    return s

def str2perm(s: str): # takes a string with a certain ordering of the digits 0 to 9 and constructs a list out of it, where each element in the list correspond to a digit + 1 in the order they appear in the string
    p = [ord(c) + 1 - ord("0") for c in s]
    return p

# # This is going to be very ugly code, I should improve it when I find time and ways to do it

# # all_perm = set()

# # dist_0 = set()
# # dist_0.add(perm2str(p0))
# # all_perm.add(perm2str(p0))

# # with open("perm_dist_10.txt", "a") as f:
# #     dist_1 = set()
# #     i = 0
# #     while i < 9:
# #         j = i + 1
# #         while j < 10:
# #             p = perm2str(perm_rev(p0, i, j))
# #             dist_1.add(p)
# #             all_perm.add(p)
# #             j += 1
# #         i += 1

# #     print("dist_1",file=f)
# #     print(dist_1,file=f)
# #     print("done")

# #     dist_2 = set()
# #     for s in list(dist_1):
# #         q = str2perm(s)
# #         i = 0
# #         while i < 9:
# #             j = i + 1
# #             while j < 10:
# #                 p = perm2str(perm_rev(q, i, j))
# #                 if p not in all_perm:
# #                     dist_2.add(p)
# #                     all_perm.add(p)
# #                 j += 1
# #             i += 1

# #     print("dist_2",file=f)
# #     print(dist_2,file=f)
# #     print("done")

# #     dist_3 = set()
# #     for s in list(dist_2):
# #         q = str2perm(s)
# #         i = 0
# #         while i < 9:
# #             j = i + 1
# #             while j < 10:
# #                 p = perm2str(perm_rev(q, i, j))
# #                 if p not in all_perm:
# #                     dist_3.add(p)
# #                     all_perm.add(p)
# #                 j += 1
# #             i += 1

# #     print("dist_3",file=f)
# #     print(dist_3,file=f)
# #     print("done")

# #     dist_4 = set()
# #     for s in list(dist_3):
# #         q = str2perm(s)
# #         i = 0
# #         while i < 9:
# #             j = i + 1
# #             while j < 10:
# #                 p = perm2str(perm_rev(q, i, j))
# #                 if p not in all_perm:
# #                     dist_4.add(p)
# #                     all_perm.add(p)
# #                 j += 1
# #             i += 1

# #     print("dist_4",file=f)
# #     print(dist_4,file=f)
# #     print("done")

# #     dist_5 = set()
# #     for s in list(dist_4):
# #         q = str2perm(s)
# #         i = 0
# #         while i < 9:
# #             j = i + 1
# #             while j < 10:
# #                 p = perm2str(perm_rev(q, i, j))
# #                 if p not in all_perm:
# #                     dist_5.add(p)
# #                     all_perm.add(p)
# #                 j += 1
# #             i += 1

# #     print("dist_5",file=f)
# #     print(dist_5,file=f)
# #     print("done")

# #     dist_6 = set()
# #     for s in list(dist_5):
# #         q = str2perm(s)
# #         i = 0
# #         while i < 9:
# #             j = i + 1
# #             while j < 10:
# #                 p = perm2str(perm_rev(q, i, j))
# #                 if p not in all_perm:
# #                     dist_6.add(p)
# #                     all_perm.add(p)
# #                 j += 1
# #             i += 1

# #     print("dist_6",file=f)
# #     print(dist_6,file=f)
# #     print("done")

# #     dist_7 = set()
# #     for s in list(dist_6):
# #         q = str2perm(s)
# #         i = 0
# #         while i < 9:
# #             j = i + 1
# #             while j < 10:
# #                 p = perm2str(perm_rev(q, i, j))
# #                 if p not in all_perm:
# #                     dist_7.add(p)
# #                     all_perm.add(p)
# #                 j += 1
# #             i += 1

# #     print("dist_7",file=f)
# #     print(dist_7,file=f)
# #     print("done")

# #     dist_8 = set()
# #     for s in list(dist_7):
# #         q = str2perm(s)
# #         i = 0
# #         while i < 9:
# #             j = i + 1
# #             while j < 10:
# #                 p = perm2str(perm_rev(q, i, j))
# #                 if p not in all_perm:
# #                     dist_8.add(p)
# #                     all_perm.add(p)
# #                 j += 1
# #             i += 1

# #     print("dist_8",file=f)
# #     print(dist_8,file=f)
# #     print("done")

# #     dist_9 = set()
# #     for s in list(dist_8):
# #         q = str2perm(s)
# #         i = 0
# #         while i < 9:
# #             j = i + 1
# #             while j < 10:
# #                 p = perm2str(perm_rev(q, i, j))
# #                 if p not in all_perm:
# #                     dist_9.add(p)
# #                     all_perm.add(p)
# #                 j += 1
# #             i += 1

# #     print("dist_9",file=f)
# #     print(dist_9,file=f)
# #     print("done")

# #     dist_10 = set()
# #     for s in list(dist_9):
# #         q = str2perm(s)
# #         i = 0
# #         while i < 9:
# #             j = i + 1
# #             while j < 10:
# #                 p = perm2str(perm_rev(q, i, j))
# #                 if p not in all_perm:
# #                     dist_10.add(p)
# #                     all_perm.add(p)
# #                 j += 1
# #             i += 1

# #     print("dist_10",file=f)
# #     print(dist_10,file=f)
# #     print("done")


# # Now we just literally compare it manually

# for i in range(len(ps)):
#     print(perm2str(perm_comp(perm_inv(qs[i]), ps[i])))

# # (does not work properly)






# Matching Random Motifs

# 
# The probability that at least one of the strings equals s is the complement of none equaling s. Since they are random and independent,
# # the probability that none is equal to s is the product of the probabilities that each is not equal to s.
# # The probability that a string is not equal to s is the probability that one of their characters is different, so 1 - probability that none of
# # the characters are different
# # The probability that no characters are different is the probability of constructing precisely that string, which given the GC content is straightforward


# file_path = "C:\\Users\\saxel\\Downloads\\rosalind_rstr.txt"
# with open(file_path, 'r') as file:
#     s = file.read()

# s = s.split()

def str2int(s: str): # converts a string into an integer
    temp = 0
    for c in s:
        temp = temp * 10
        temp += ord(c) - ord("0")
    return temp

def str2dec(s: str): # converts a string into a decimal number
    after_comma = False
    decimals = 0
    temp = 0
    for c in s:
        if after_comma == True:
            decimals += 1
        if c == ".":
            after_comma = True
        else:
            temp = temp * 10
            temp += ord(c) - ord("0")
    return temp / (10 ** decimals)


# N = str2int(s[0])
# x = str2dec(s[1])
# gen = s[2]

# # Compute the probability of a random string being constructed with that GC content:

# gc = x / 2 # probability of a nucleotide being G or C
# at = (1 - x) / 2 # probability of a nucleotide being A or T

# p = 1
# for c in gen:
#     if c in {"A", "T"}:
#         p = p * at
#     elif c in {"C", "G"}:
#         p = p * gc

# # Probability that a string is not equal to s
# p = 1 - p

# # Probability that all strings are not equal to s
# p = p ** N

# # Probability that at least one string is equal to s
# p = 1 - p
# print(p)



# Counting Subsets

# file_path = "C:\\Users\\saxel\\Downloads\\rosalind_sset.txt"
# s = num_list_read(file_path)

# n = s[0]

# print((2 ** n) % (10 ** 6))





# Introduction to Alternative Splicing

# file_path = "C:\\Users\\saxel\\Downloads\\rosalind_aspc.txt"
# s = num_list_read(file_path)

# n = s[0]
# m = s[1]

# temp = 0
# while m <= n:
#     temp += over(n, m) % (10 ** 6)
#     m += 1
# print(temp % (10 ** 6))




# Edit Distance

# file_path = "C:\\Users\\saxel\\Downloads\\rosalind_edit (2).txt"
# keys, gena = fasta_read(file_path)

# s = gena[0]
# t = gena[1]

# # This problem is equivalent to computing the length of the longest common substring. 
# # That substring stays fixed and we edit operations around it to get from one string to the other

# m = len(s)
# n = len(t)

# # Wlog we take s to be the longest string

# if n > m:
#     r = t
#     t = s
#     s = r
#     m = len(s)
#     n = len(t)

# memo = [["" for _ in range(n + 1)] for _ in range(m + 1)]

# u = LCS(s, t, m, n, memo)

# diff = 0

# # Add what it takes to get from u to the longest
# diff += max(m, n) - len(u)
# print(diff)

# # We also count the spaces between appearances of characters of u in both s and t, because if the spaces are larger in t, we have to eliminate the extra ones
# # before swapping the rest to those of s

# gaps_s = []
# gaps_t = []

# i = 0
# j = 0
# for c in u:
#     temp_s = 0
#     temp_t = 0
#     while s[i] != c:
#         temp_s += 1
#         i += 1
#     while t[j] != c:
#         temp_t += 1
#         j += 1
#     gaps_s.append(temp_s)
#     gaps_t.append(temp_t)
#     i += 1
#     j += 1
# gaps_s.append(m - i)
# gaps_t.append(n - j)

# for i in range(len(gaps_s)):
#     if gaps_s[i] < gaps_t[i]:
#         diff += gaps_t[i] - gaps_s[i]
#         print(diff)

# print(gaps_s)
# print(gaps_t)

# # probably not the correct idea





# Expected Number of Restriction Sites

# # Suppose m = len(t). Then the string t has n - m + 1 possibilities of appearing as a substring of a string of length n
# # For each of these possibilities, the probability of it being the string t is given by the GC content as in previous exercises

# file_path = "C:\\Users\\saxel\\Downloads\\rosalind_eval.txt"
# with open(file_path, 'r') as file:
#     s = file.read()

# s = s.split()
# n = str2int(s.pop(0))
# t = s.pop(0)

# B = [0.0 for _ in range(len(s))]
# for i in range(len(s)):
#     s[i] = str2dec(s[i])
#     gc = s[i] / 2
#     at = (1 - s[i]) / 2
#     # Compute the probability that a string of length len(t) is exactly t with the specific GC content given
#     p = 1
#     for c in t:
#         if c in {"A", "T"}:
#             p = p * at
#         elif c in {"C", "G"}:
#             p = p * gc
#     # The expected number of times t is a substring will be the sum of the probabilities of the string t appearing in position 1, 2, 3, 4, ..., n - len(t) + 1 (positions are 1-indexed)
#     # Each appearance is independent of each other (that is why we add them) and they all have the same probability just computed
#     B[i] = (n - len(t) + 1) * p
#     print(B[i], end=" ")




# Distances in Trees

# file_path = "C:\\Users\\saxel\\Downloads\\rosalind_nwck (2).txt"
# with open(file_path, 'r') as file:
#     s = file.read()

# s = s.split()
# n = len(s) // 3
# trees = [s[3 * i] for i in range(n)]
# starts = [s[3 * i + 1] for i in range(n)]
# ends = [s[3 * i + 2] for i in range(n)]

# # The different levels of the tree are separated by "(" and ")": the deepest level is the list separated by commas within the most nested expression between parentheses
# # Elements in the same level share a "parent", so they are at distance 2; they are at distance one with their children or parents, e.g. with nodes within a level for which they stay at the right of their closing ")"
# # or if they are in such a parenthesis, then to the element at the right of their closing parenthesis; so we have to count these parenthesis

# # Another way, maybe useful at some point, would be to just have a function translate a tree in Fenwick notation to the usual binary tree data structure

# for i in range(n):
#     tree = trees[i]
#     start = starts[i]
#     end = ends[i]
#     # We first cut the tree to only contain the part of the Newick string between start and end: the rest is irrelevant
#     i = 0
#     pruned_tree = ""
#     interesting = False
#     while tree[i] != ";":
#         if tree[i] in {"(", ")", ","}:
#             if interesting == True:
#                 pruned_tree += tree[i]
#             i += 1
#         else:
#             curr = ""
#             while i < len(tree) and tree[i] not in {"(", ")", ",", ";"}:
#                 curr += tree[i]
#                 i += 1
#             if curr in {start, end}:
#                 interesting = not interesting
#     dist = 0
#     # Now we interpret this info:
#     # If we find paired parentheses (...), these are irrelevant, since they refer to other branchings
#     # For the unpaired parentheses, ")" means what is at the left is child of what its at the right, so each of them increases the distance by 1
#     # Otherwise, "(" will always come with commas before, e.g. "A,(D" , and means A and the parent of D have the same parent, so their distance is 3
#     # However, if more than one of these "(" come in a row, e.g. "A,((D", then it means A and the parent of the parent of D have the same parent, so at that point the distance is 4:
#     # the general rule is then that the first opening parenthesis adds distance 3 and the subsequent ones add each distance 1
#     # If at the end there is a comma, this means that what is before the comma and the node after the comma have the same parent, meaning they are at distance 2, so we add that. Other commas are irrelevant since they affect sibling-relationships between nodes where the start or end node might be contained, but not the relationship between the start and end node
#     parentheses = []
#     for c in pruned_tree:
#         if c == "(":
#             parentheses.append("(")
#         elif c == ")":
#             if parentheses and parentheses[-1] == "(":
#                 parentheses.pop(-1)
#             else:
#                 parentheses.append(")")
#     for i in range(len(parentheses)):
#         if parentheses[i] == ")":
#             dist += 1
#         elif parentheses[i] == "(":
#             if i > 0 and parentheses[i - 1] == "(":
#                 dist += 1
#             else:
#                 dist += 3
#     if pruned_tree[-1] == "," and (not parentheses or parentheses[-1] == ")"):
#         dist += 2
#     print(dist, end=" ")





# Interleaving Two Motifs

# file_path = "C:\\Users\\saxel\\Downloads\\rosalind_scsp.txt"
# with open(file_path, 'r') as file:
#     r = file.read()

# r = r.split()
# s = r[0]
# t = r[1]

# def SCS(s: str, t: str): # Constructs the shortest common supersequence of s and t
#     m = len(s)
#     n = len(t)
#     # Look for the longest common subsequence
#     memo = [["" for _ in range(n + 1)] for _ in range(m + 1)]
#     u = LCS(s, t, m, n, memo)
#     # Suppose for a certain k that s[i] == t[j] == u[k] and s[i'] == t[j'] == u[k + 1]. Then between u[k] and u[k + 1] add first all the s[i < l < i'] and then all the t[j < l' < j']
#     # This produces a shortest common subsequence, since if there was a coincidence between an s[l] and a t[l'] in those ranges, then it is common to s and t and thus should be in the longest common subsequence
#     i = 0
#     j = 0
#     k = 0
#     scs = ""
#     while k < len(u):
#         while s[i] != u[k]:
#             scs += s[i]
#             i += 1
#         while t[j] != u[k]:
#             scs += t[j]
#             j += 1
#         scs += u[k]
#         k += 1
#         i += 1
#         j += 1
#     # Add also the ends of the strings if they are left to add
#     while i < m:
#         scs += s[i]
#         i += 1
#     while j < n:
#         scs += t[j]
#         j += 1
#     return scs

# with open("output.txt", "a") as f:
#     print(SCS(s, t), file=f)





# Introduction to Set Operations

# file_path = "C:\\Users\\saxel\\Downloads\\rosalind_seto.txt"
# with open(file_path, 'r') as file:
#     s = file.read()

# s = s.split()
# n = str2int(s[0])
# A = set()
# B = set()
# A.add(str2int(s[1][1:-1])) # the reason for using this form is because we might have numbers with more than one digit and we just want to get rid of "{", "," or "}"
# i = 2
# while s[i][-1] != "}":
#     A.add(str2int(s[i][:-1]))
#     i += 1
# A.add(str2int(s[i][:-1]))
# B.add(str2int(s[i + 1][1:-1]))
# i += 2
# while s[i][-1] != "}":
#     B.add(str2int(s[i][:-1]))
#     i += 1
# B.add(str2int(s[i][:-1]))

# X = set() # set of n elements
# for i in range(n):
#     X.add(i + 1)

# def set_union(A: set, B: set): # constructs the union of the set A and B
#     C = set()
#     for a in list(A):
#         C.add(a)
#     for b in list(B):
#         C.add(b)
#     return C

# def set_intersection(A: set, B: set): # constructs the intersection of the sets A and B
#     C = set()
#     for a in list(A):
#         if a in B:
#             C.add(a)
#     return C

# def set_diff(A: set, B: set): # constructs the difference of sets A - B
#     C = set()
#     for a in list(A):
#         if a not in B:
#             C.add(a)
#     return C

# def set_comp(A: set, X: set): # constructs the complement of the set A as a subset of X
#     return set_diff(X, A)

# with open("output.txt", "a") as f:
#     print(set_union(A, B), file=f)
#     print(set_intersection(A, B), file=f)
#     print(set_diff(A, B), file=f)
#     print(set_diff(B, A), file=f)
#     print(set_comp(A, X), file=f)
#     print(set_comp(B, X), file=f)




# Inferring Protein from Spectrum

# # We just need to compute the difference between two consecutive numbers, figure which isotope corresponds to that mass and add it to the resulting string 

# file_path = "C:\\Users\\saxel\\Downloads\\rosalind_spec.txt"
# with open(file_path, 'r') as file:
#     s = file.read()
# file_path = "C:\\Users\\saxel\\Downloads\\rosalind_monoisot.txt"
# with open(file_path, "r") as file:
#     r = file.read()

# r = r.split()
# isotopes = []
# masses = []
# for i in range(len(r) // 2):
#     isotopes.append(r[2 * i])
#     masses.append(str2dec(r[2 * i + 1]))

# s = s.split()
# s[0] = str2dec(s[0])
# protein = ""

# for i in range(len(s) - 1):
#     s[i + 1] = str2dec(s[i + 1])
#     diff = s[i + 1] - s[i]
#     j = 0
#     while abs(diff - masses[j]) > 0.0001:
#         j += 1
#     protein += isotopes[j]
# print(protein)





# Introduction to Pattern Matching

# file_path = "C:\\Users\\saxel\\Downloads\\rosalind_trie.txt"
# with open(file_path, 'r') as file:
#     s = file.read()

# s = s.split()
# trie = [[i + 1, i + 2, s[0][i]] for i in range(len(s[0]))]
# for j in range(len(s) - 1):
#     k = 0 # position in the trie
#     i = 0 # character in current gene
#     curr = 1 # parent node
#     # Check which path in the tree to follow for the new gene; we look at the parents and the letters at the edges in each step, and the next parent is the child of the previous one
#     while k < len(trie):
#         if [trie[k][0], trie[k][2]] == [curr, s[j + 1][i]]: # s[j + 1] because we already dealt with s[0]
#             curr = trie[k][1]
#             i += 1
#         else:
#             k += 1
#     # When we already found all coincidences, we just need to add the new nodes, taking into account in which one we landed last and adding the rest as new with indexes related to the lenght of the trie
#     while i < len(s[j + 1]):
#         trie.append([curr, len(trie) + 2, s[j + 1][i]])
#         curr = len(trie) + 1
#         i += 1
# with open("output.txt", "a") as f:
#     for t in trie:
#         print(t[0], end=" ", file=f)
#         print(t[1], end=" ", file=f)
#         print(t[2], file=f)





# Comparing Spectra with the Spectral Convolution

# file_path = "C:\\Users\\saxel\\Downloads\\rosalind_conv.txt"
# with open(file_path, 'r') as file:
#     s = file.read()

# s = s.splitlines()
# s1 = s[0].split()
# s2 = s[1].split()
# for i in range(len(s1)):
#     s1[i] = str2dec(s1[i])
# for i in range(len(s2)):
#     s2[i] = str2dec(s2[i])

# # We just need to construct the spectral convolution and see which element repeats the most, up to 4 significant figures    
# spec_conv = []
# for n1 in s1:
#     for n2 in s2:
#         spec_conv.append(n1 - n2)
# spec_conv.sort()
# i = 1
# curr = spec_conv[1]
# mult = 1
# x = 0
# max_mult = 0
# while i < len(spec_conv):
#     if abs(spec_conv[i] - curr) < 0.0001 :
#         mult += 1
#     else:
#         if mult > max_mult:
#             max_mult = mult
#             x = curr
#         curr = spec_conv[i]
#         mult = 1
#     i += 1
# print(max_mult)
# print(x)




# Creating a Character Table

file_path = "C:\\Users\\saxel\\Downloads\\rosalind_check.txt"
with open(file_path, 'r') as file:
    s = file.read()

# In the Newick format, if we picture the rooted tree as having nodes separated in horizontal lines, where the root is at height 0, and the children of nodes of height i are at height i + 1,
# then these heights correspond to the number of opening parentheses that are unclosed

# For example, in the sample dataset, dog and cat appear after 1 unclosed opening parenthesis, and so are at height 1; robot is at height 2 and elephant and mouse at height 3
# We can use this height structure to keep track of what gets separated when cutting a branch of the tree
# We save the positions after which the parentheses appear and how many of them are unclosed

i = 0
openings = [] # these just tells you the height of each taxon
siblings = [] # these merge together the taxa that are siblings (since these need to have the same digit)
curr = 0
while i < len(s):
    if s[i] == "(":
        curr += 1
    elif s[i] == ")":
        curr -= 1
    elif i > 0 and s[i] not in {"(", ")", ",", ";"} and s[i - 1] in {"(", ")", ",", ";"}:
        openings.append([i, curr])
        if i == 1 or s[i - 2] in {"(", ")", ",", ";"}:
            siblings.append([i, curr])
    i += 1

n = len(openings) # the number of taxa
print(openings)
print(siblings)
for option in itertools.product({0, 1}, repeat = n):
    k = sum(option)
    # if k > 1 and k < n - 1:
        # still to be done
        # do not forget to order the taxa lexicographically when printing!
