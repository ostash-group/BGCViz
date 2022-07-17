from Bio import SeqIO
import sys
import os
import pandas as pd


def write_gbs(group_by, data, label, seq_file):
	counter = 0
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
				file =  SeqIO.parse(open(seq_file), "genbank")
				print("Working on: "+ label+"_"+"cluster_"+str(list_l[i])+"_"+str(group))
				for record in file:
					loci = [feat for feat in record.features if feat.type == "CDS"]
					start_new = int(start[i])
					end_new = int(stop[i])
					subrecord = record[start_new:end_new]
					annotation={"molecule_type":"DNA"}
					subrecord.annotations = annotation
					SeqIO.write(subrecord, group+"/"+label+"_"+"cluster_"+str(list_l[i])+"_"+str(group)+".gb", "genbank")

def main():
	group_by = pd.read_csv("group_by.csv", dtype = str)
	seq_file = sys.argv[1]

	print("Searching for antismash files...")
	if os.path.exists("anti_biocircos.csv"):
		print("Found!")
		data = pd.read_csv("anti_biocircos.csv")
		label = "Antismash"
		write_gbs(group_by, data, label, seq_file)

	print("Searching for deepbcg files...")
	if os.path.exists("deep_biocircos.csv"):
		print("Found!")
		data = pd.read_csv("deep_biocircos.csv")
		label = "DeepBGC"
		write_gbs(group_by, data, label, seq_file)

	print("Searching for prism files...")
	if os.path.exists("prism_biocircos.csv"):
		print("Found!")
		data = pd.read_csv("prism_biocircos.csv")
		label = "PRISM"
		write_gbs(group_by, data, label, seq_file)

	print("Searching for rre-finder files...")
	if os.path.exists("rre_biocircos.csv"):
		print("Found!")
		data = pd.read_csv("rre_biocircos.csv")
		label = "RRE-Finder"
		write_gbs(group_by, data, label, seq_file)

	print("Searching for sempi files...")
	if os.path.exists("sempi_biocircos.csv"):
		print("Found!")
		data = pd.read_csv("sempi_biocircos.csv")
		label = "SEMPI"
		write_gbs(group_by, data, label, seq_file)

	print("Searching for ARTS files...")
	if os.path.exists("arts_biocircos.csv"):
		print("Found!")
		data = pd.read_csv("arts_biocircos.csv")
		label = "ARTS"
		write_gbs(group_by, data, label, seq_file)

	print("Searching for PRISM supplement files...")
	if os.path.exists("prism_supp_biocircos.csv"):
		print("Found!")
		data = pd.read_csv("prism_supp_biocircos.csv")
		label = "PRISM-Supp"
		write_gbs(group_by, data, label, seq_file)

	print("Searching for GECCO files...")
	if os.path.exists("gecco_biocircos.csv"):
		print("Found!")
		data = pd.read_csv("gecco_biocircos.csv")
		label = "GECCO"
		write_gbs(group_by, data, label, seq_file)



if __name__ == "__main__":
	main()

