from Bio import SeqIO
import sys
import os
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

def write_gbs(group_by, data, label, seq_file, args):
	counter = 0 
	fl = convert_gbff(seq_file)
	with open(fl, "r") as handle:
		record =  SeqIO.read(handle, "genbank")
		loci = [feat for feat in record.features if feat.type == "CDS"]
		for  index, row in pd.DataFrame(group_by[label].dropna()).iterrows():
			counter += 1
			if (counter >= len(pd.DataFrame(group_by[label].dropna()).index)):
					break
			else:
				start = []
				stop = []
				list_l = row[label].split(",")
				for i in range(len(list_l)):
					list_l[i] = int(list_l[i])
					start.append(data[list_l[i] == data['Cluster']].Start.item())
					stop.append(data[list_l[i] == data['Cluster']].Stop.item())
					group = group_by.Group[index]
					if os.path.isdir(group):
						pass
					else:
						os.mkdir(group)
					if not args.quiet:
						print("Working on: "+ label+"_"+"cluster_"+str(list_l[i])+"_"+str(group))
					if os.path.exists(group+"/"+label+"_"+"cluster_"+str(list_l[i])+"_"+str(group)+".gb") and not args.force:
						if counter <=1: 
							print("Files exist! Please use --force option to override them")
						continue
					start_new, end_new = solve_incomplete_CDS(int(start[i]), int(stop[i]), loci)
					subrecord = record[start_new:end_new]
					annotation={"molecule_type":"DNA"}
					subrecord.annotations = annotation
					SeqIO.write(subrecord, group+"/"+label+"_"+"cluster_"+str(list_l[i])+"_"+str(group)+".gb", "genbank")

					
def group_gb_files(group_by, seq_file, args):
	data_to_search = {
		"antiSMASH" : ("anti_biocircos.csv", "Antismash"),
		"DeepBGC" : ("deep_biocircos.csv","DeepBGC" ),
		"PRISM" : ("prism_biocircos.csv", "PRISM"),
		"RREfinder" : ("rre_biocircos.csv", "RRE-Finder"),
		"SEMPI" : ("sempi_biocircos.csv", "SEMPI"),
		"ARTS" : ("arts_biocircos.csv", "ARTS"),
		"PRISM supplement" : ("prism_supp_biocircos.csv", "PRISM-Supp"),
		"GECCO" : ("gecco_biocircos.csv", "GECCO")
	}
	for k,v in data_to_search.items():
		if not args.quiet:
			print("Searching for "+str(k)+" files...")
		if os.path.exists(v[0]):
			if not args.quiet:
				print("Found!")
			data = pd.read_csv(v[0])
			label = v[1]
			write_gbs(group_by, data, label, seq_file, args)


def run_clinker(group_by, args):
	for  index, row in pd.DataFrame(group_by["Group"]).iterrows():
		group = group_by.Group[index]
		if os.path.isdir("clinker_plots"):
			pass
		else:
			os.mkdir("clinker_plots")
		os.system(f"clinker {group} --plot clinker_plots/{group}.html  -i 0.9 -j {args.jobs}")

def main():
	# Reading data
	group_by = pd.read_csv("group_by.csv", dtype = str)

	# Parsing arguments
	parser = argparse.ArgumentParser(description='Small helper script for BGCViz')
	parser.add_argument("-i", "--input", help="Input .gb/.gbk/.gbff file. One record per file will be used (as one genome). Required")
	parser.add_argument("--force", help="Force overwrite calculated results. [default = False]",action=argparse.BooleanOptionalAction)
	parser.add_argument("--quiet", help="Run silently. Clinker will run as usual. [default = False]",action=argparse.BooleanOptionalAction)
	parser.add_argument("-cl", "--run_clinker", help="Automatically runs clinker on groups. Results are stored in 'clinker_plots' folder. [default = False]", 
	action=argparse.BooleanOptionalAction)
	parser.add_argument("-j", "--jobs", help="Number of threads for clinker analysis (0=all). [default = 0] ", default=0)
	args = parser.parse_args()

	# Run grouping for gb files
	group_gb_files(group_by, args.input, args)

	#Run clinker
	if args.run_clinker:
		run_clinker(group_by, args)

	# Bye message
	print("Analysis finished successfuly!")


if __name__ == "__main__":
	main()


