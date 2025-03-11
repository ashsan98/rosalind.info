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
# letters.sort()

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

# def catalan(n): # Constructs the n-th Catalan number
#     # The following will be the dynamic programming solution, but we can derive a closed formula from the recursion and this will be easier to compute
#     # if n == 0:
#     #     return 1
#     # temp = 0
#     # i = 1
#     # while i <= n:
#     #     temp += (catalan(i - 1) * catalan(n - i))
#     #     i += 1
#     # return temp

#     # Using the closed formula:
#     if n == 0:
#         return 1
#     return factorial(2 * n) // (factorial(n) * factorial(n + 1))

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

file_path = "C:\\Users\\saxel\\Downloads\\rosalind_corr.txt"
keys, gena = fasta_read(file_path)

def Hamming_dist(s, t): # computes the Hamming distance between gen-strings of equal length
    if len(s) != len(t):
        raise ValueError("Error: the gens need to have the same length as strings")
    dist = 0
    for i in range(len(s)):
        if s[i] != t[i]:
            dist += 1
    return dist

# We first sort the genes between those repeated and those not repeated: the not repeated ones are the defect reads
correct_gena = []
incorrect_gena = []
while len(gena) >= 2:
    curr_gen = gena.pop(0)
    curr_gen_revcomp = dna_revcomp(curr_gen) # we consider equivalence modulo reverse complement
    matching = 0 # keeps count of whether we find the same gene read again
    i = 0
    while i < len(gena):
        if gena[i] == curr_gen or gena[i] == curr_gen_revcomp:
            gena.pop(i)
            matching += 1
            i -= 1
        i += 1
    if matching > 0:
        correct_gena.append(curr_gen)
    else:
        incorrect_gena.append(curr_gen)
if len(gena) == 1: # if only one is left, there is non to compare it to
    incorrect_gena.append(gena.pop(0))

# We compare the defect reads with the correct ones to find which among them has Hamming distance one to the defect
with open("output.txt", "a") as f:
    for gen in incorrect_gena:
        j = 0
        while j < len(correct_gena):
            if Hamming_dist(gen, correct_gena[j]) == 1:
                print(gen + "->" + correct_gena[j], file=f)
                j = len(correct_gena)
            elif Hamming_dist(gen, dna_revcomp(correct_gena[j])) == 1:
                print(gen + "->" + dna_revcomp(correct_gena[j]), file=f)
                j = len(correct_gena)
            j += 1
    