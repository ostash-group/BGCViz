from Bio import SeqIO
import sys
import os
import re
import pandas as pd
import argparse



def convert_gbff(seq_file):
	if seq_file.split('.')[-1] == 'gbff':
		if os.path.exists(os.path.splitext(os.path.basename(seq_file))[0]+'.gbk'):
			converted = os.path.splitext(os.path.basename(seq_file))[0]+'.gbk'
		else:
			file_name = os.path.basename(seq_file)
			with open(seq_file, "r") as f:
				for index, record in enumerate(SeqIO.parse(f, "genbank")):
					if index > 0:
						break
					else:
						SeqIO.write(record, os.path.splitext(os.path.basename(seq_file))[0]+'.gbk', "genbank")
			converted = os.path.splitext(os.path.basename(seq_file))[0]+'.gbk'
	else:
		converted = seq_file
	return converted

def solve_incomplete_CDS(start_new, end_new, loci):
	for feat in loci:
		if feat.location.start.position <= start_new and feat.location.end.position > start_new:
			start_new = feat.location._start.position
		if feat.location.start.position <= end_new and feat.location.end.position > end_new:
			end_new = feat.location.end.position
	return start_new, end_new

def write_gbs(data,seq_file, fl_name):
	counter = 0
	fl = convert_gbff(seq_file)
	with open(fl, "r") as handle:
		record =  SeqIO.read(handle, "genbank")
		loci = [feat for feat in record.features if feat.type == "CDS"]
		locus_to_separate=data['separate_before'].split(',')
		loci_len = len(locus_to_separate)
		for to_split in range(loci_len):
			for loc in range(len(loci)):
				if loci[loc].qualifiers['locus_tag'][0] == locus_to_separate[to_split]:
					if to_split == 0:
						start = 0
						stop = loci[loc-1].location._end.position
						subrecord = record[start:stop]
						SeqIO.write(subrecord, "dissected/"+os.path.splitext(fl_name)[0]+"_subcluster_"+str(counter+1)+".gb", "genbank")
						counter +=1
					if to_split == (loci_len-1):
						start = loci[loc].location._start.position
						stop = loci[len(loci)-1].location._end.position
					else:
						start = loci[loc].location._start.position
						stop = [loci[loc-1].location._end.position for loc in range(len(loci)) if loci[loc].qualifiers['locus_tag'][0] == locus_to_separate[to_split+1]]
					subrecord = record[start:stop]
					SeqIO.write(subrecord, "dissected/"+os.path.splitext(fl_name)[0]+"_subcluster_"+str(counter+1)+".gb", "genbank")
					counter +=1


					
def split_gb_files(csv_file):
	data_to_search = {
		"antismash" :  "Antismash",
		"deepbgc" : "DeepBGC" ,
		"prism" : "PRISM",
		"sempi" :  "SEMPI",
		"gecco" :  "GECCO"
	}
	data = pd.read_csv(csv_file)
	for index, row in data.iterrows():
		if pd.isna(row['group']) : continue
		group_name = str('group_' + str(int(row['group'])))
		if os.path.isdir(group_name):
			pass
		else:
			print("Didn't found group: "+ group_name)
		try:
			file_to_open = data_to_search[row['by_software']]
		except KeyError():
			print("Could not find software for "+str(group_name)+" : "+str(row['by_software']))
		r = re.compile(file_to_open+".*")
		hit = list(filter(r.match, os.listdir(group_name)))[0]
		file_to_open=group_name+"/"+hit
		write_gbs(row, file_to_open, hit)


def main():
	if not os.path.isdir('dissected'):
		os.mkdir('dissected')
	# Parsing arguments
	parser = argparse.ArgumentParser(description='Small helper script for BGCViz')
	required = parser.add_argument_group('Required arguments')
	required.add_argument("-i", "--input", help=".csv file with clusters to separate", required = True)
	args = parser.parse_args()

	# Run grouping for gb files
	split_gb_files(args.input)
	# Bye message
	print("Separation of clusters finished successfuly!")
	print("Please find genbank files in 'dissected' folder")


if __name__ == "__main__":
	main()


