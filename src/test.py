#!/usr/bin/env python3
import sys
import Bio

"""
This module is written for the final exam of python for genomic data science on coursera.
"""

# process the fasta file, return a dictionary 
def fasta_to_dic(filename):
	''' Process the fasta file, extract id and sequence, return as a dictionary'''
	try:
		f=open(filename)
	except IOError:
		print('File %s not found' % filename)

	seqs = {}
	for line in f:
		line=line.rstrip()
		if line.startswith('>'):
			name=line.split()[0][1:]
			#print(name)
			seqs[name] = ''
		else:
			seqs[name] = seqs[name]+line
	return seqs
	f.close()


# make a length dictionary 
def count_seq_lenth(seqs):
	length_dict = {}
	for key,value in seqs.items():
		length_dict[key] = len(value)
	return length_dict



def find_ORF(dna):
	"""find all possible ORF given a dna sequence. return a list of 
	[frame, start position, stop position, ORF length, ORF seq, strand];
	dna : a DNA sequence """
	start_coden = ['ATG','atg']
	stop_coden = ['TAG', 'TAA', 'TGA','tag','taa','tga']
	frames = [1,2,3]
	orf_list = []

	# forwar strand
	for frame in frames:
		#print('Frame: ', frame)
		for i in range(frame-1, len(dna), 3):
			if dna[i:i+3] in start_coden:
				start = i
				#print("start: ", start)
				for j in range(i+3, len(dna),3):
					if dna[j:j+3] in stop_coden:
						stop = j+3
						#print("stop: ", stop)
						ORF=dna[start:stop]
						length = len(ORF)
						strand = 'forward'
						new_orf = [frame, start, stop, length, ORF, strand]
						#print('new_orf: ', new_orf)
						orf_list.append(new_orf)
						break # when find the first stop coden, then break the loop. 
	
	# reverse strand
	from Bio.Seq import Seq
	dna=str(Seq.reverse_complement(Seq(dna)))
	for frame in frames:
		#print('Frame: ', frame)
		for i in range(frame-1, len(dna), 3):
			if dna[i:i+3] in start_coden:
				start = i
				#print("start: ", start)
				for j in range(i+3, len(dna),3):
					if dna[j:j+3] in stop_coden:
						stop = j+3
						#print("stop: ", stop)
						ORF=dna[start:stop]
						length = len(ORF)
						strand = 'reverse'
						new_orf = [frame, start, stop, length, ORF, strand]
						#print('new_orf: ', new_orf)
						orf_list.append(new_orf)
						break # when find the first stop coden, then break the loop. 

	return orf_list



def find_repeats(dna, n):
	""" find all repeats of length n in the dna sequence, return a dictionary of repeats and counts of the repeat"""
	if n >= float(len(dna))/2 :
		print("Repeat length is more than half of the dna length. Not able to find repeats")
	else:
		# 1. get all possiber n-mer
		nmer = set()
		for i in range(len(dna)-n+1):
			nmer.add(dna[i:i+n])
		
		# 2. for each unique n-mer, search one by one through the sequence
		repeats={}
		for item in nmer:
			count=0
			for i in range(len(dna)-n +1):
				if dna[i:i+n] == item:
					count += 1
				repeats[item] = count 
		return(repeats)



























