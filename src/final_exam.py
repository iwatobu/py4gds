#!/usr/bin/env python3

import test
import pandas as pd

testfile='../data/dna2.fasta'
seqs = test.fasta_to_dic(testfile)

'''
How many records are in the multi-FASTA file?
What is the length of the longest sequence in the file?
What is the length of the shortest sequence in the file?
'''

print('records_number in the file : %d' % len(seqs))
# bash grep -e '^>' dna.example.fasta| wc -l

lengths = test.count_seq_lenth(seqs)
print(lengths.values())
l=max(lengths.values())
s=min(lengths.values())
print("longest length is %d" % l)
print("shortest length is %d" % s)
for key,value in lengths.items():
	if value == l:
		print("longest seq id is %s" % key)
	elif value == s:
		print("shortest seq id is %s" % key)


'''
What is the length of the longest ORF appearing in reading frame 2 of any of the sequences?
What is the starting position of the longest ORF in reading frame 3 in any of the sequences? 
The position should indicate the character number where the ORF begins. 
For instance, the following ORF:
> sequence1
ATGCCCTAG
starts at position 1.

What is the length of the longest ORF appearing in any sequence and in any forward reading frame?
What is the length of the longest forward ORF that appears in the sequence with the identifier  gi|142022655|gb|EQ086233.1|16?
'''


orf_df = pd.DataFrame(columns = ['frame', 'start', 'stop', 'length', 'ORF','strand','identifier'])
for key,value in seqs.items():
	temp = pd.DataFrame(test.find_ORF(value), columns=['frame', 'start', 'stop', 'length', 'ORF','strand'])
	temp['identifier'] = key
	orf_df = pd.concat([orf_df, temp])
#print(orf_df.head())
#orf_df can be used to identify the longest ORF for each identifer 
forward = orf_df[orf_df.strand=='forward']

print("\nThe length of the longest ORF in frame 2 of any seq is: \n")
print(max(orf_df[orf_df.frame==2].length)) 
print(max(forward[forward.frame==2].length))

print("\nThe starting position of the longest ORF in reading frame 3 in any of the sequences: \n")
t1=orf_df[orf_df.frame==3].sort_values(['identifier','length']).groupby('identifier').last().reset_index().sort_values('length')
print(t1[['frame', 'start', 'stop', 'length','strand','identifier']].tail()) # get the start position + 1


print("\nthe length of the longest ORF appearing in any sequence and in any forward reading frame is : \n")
print(max(forward.length))

print("\n What is the length of the longest forward ORF that appears in the sequence with the identifier  gi|142022655|gb|EQ086233.1|16?? \n")
print(max(forward[forward.identifier == 'gi|142022655|gb|EQ086233.1|16'].length) )


#Find the most frequently occurring repeat of length 6 in all sequences. How many times does it occur in all?
repeat_df = pd.DataFrame(columns = ['repeat_seq', 'repeat_count','identifier'])
N=6
for key,value in seqs.items():
	dic = test.find_repeats(dna=value, n=N)
	temp = pd.DataFrame({'repeat_seq' : list(dic.keys()), 'repeat_count' : list(dic.values())})
	temp['identifier'] = key
	repeat_df = pd.concat([repeat_df, temp])

#print(repeat_df.dtypes)
repeat_df.repeat_count = pd.to_numeric(repeat_df.repeat_count)
#print(repeat_df.dtypes)
#print(repeat_df.head())
# filter repeat count >1
repeat_df = repeat_df[repeat_df.repeat_count > 1]
print(repeat_df.shape)

most_freqent_repeat = pd.DataFrame({'Occurence_time': repeat_df.groupby('repeat_seq').repeat_count.sum()}).sort_values('Occurence_time', ascending=False)
print(most_freqent_repeat.head())
most_freqent_repeat_seq = list(most_freqent_repeat.index)[0]
print(most_freqent_repeat_seq)

print("Find the most frequently occurring repeat of length 6 in all sequences. How many times does it occur in all? \n")
print(sum(repeat_df[repeat_df.repeat_seq == most_freqent_repeat_seq].repeat_count))


#Find all repeats of length 12 in the input file. Let's use Max to specify the number of copies of the most frequent repeat of length 12.  
#How many different 12-base sequences occur Max times?
repeat_df = pd.DataFrame(columns = ['repeat_seq', 'repeat_count', 'identifier'])
N=12
for key,value in seqs.items():
	dic = test.find_repeats(dna=value, n=N)
	temp = pd.DataFrame({'repeat_seq' : list(dic.keys()), 'repeat_count' : list(dic.values())})
	temp['identifier'] = key
	repeat_df = pd.concat([repeat_df, temp])

#print(repeat_df.dtypes)
repeat_df.repeat_count = pd.to_numeric(repeat_df.repeat_count)
#print(repeat_df.dtypes)
#
# filter repeat count >1
repeat_df = repeat_df[repeat_df.repeat_count > 1]
print(repeat_df.shape)
print(repeat_df.head())

most_freqent_repeat = pd.DataFrame({'Occurence_time': repeat_df.groupby('repeat_seq').repeat_count.sum()}).sort_values('Occurence_time', ascending=False)
print(most_freqent_repeat.head())

max_num = most_freqent_repeat.Occurence_time[0]
print(most_freqent_repeat[most_freqent_repeat.Occurence_time == max_num].shape)


#Which one of the following repeats of length 7 has a maximum number of occurrences?
repeat_df = pd.DataFrame(columns = ['repeat_seq', 'repeat_count', 'identifier'])
N=7
for key,value in seqs.items():
	dic = test.find_repeats(dna=value, n=N)
	temp = pd.DataFrame({'repeat_seq' : list(dic.keys()), 'repeat_count' : list(dic.values())})
	temp['identifier'] = key
	repeat_df = pd.concat([repeat_df, temp])

#print(repeat_df.dtypes)
repeat_df.repeat_count = pd.to_numeric(repeat_df.repeat_count)

# filter repeat count >1
repeat_df = repeat_df[repeat_df.repeat_count > 1]
print(repeat_df.shape)
print(repeat_df.head())


most_freqent_repeat = pd.DataFrame({'Occurence_time': repeat_df.groupby('repeat_seq').repeat_count.sum()}).sort_values('Occurence_time', ascending=False)
print(most_freqent_repeat.head())




