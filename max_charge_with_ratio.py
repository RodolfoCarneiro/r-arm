import sys

# Iterface updating the percentage of the ORFs read
def percent(index, protein_length, window, lenlista): 
	proteins = 100.0*(float(index))/float(protein_length)
	aminoacids = 100.0*float(window)/float(lenlista)

	proteins2 = int(proteins)
	aminoacids2 = int(aminoacids)
	print str(proteins2), '% of proteins complete.   ', str(aminoacids2), '% of current gene complete.    \r',

def evaluate(lista, hcharge, window, i, lista_length, min_ratio, fsequence):

	max_charge = -99999 # The highest charge found in the protein
	max_ratio = -99999 # The highest ratio found in the protein
	best_sequence = '' # sequenceuence of the highest charge or ratio

	while window <= len(lista): # Try every possible window size

		# Browse every window in the protein
		for k in range(len(lista)-window+1):

			percent(i, lista_length, window, len(lista))
			sequence = '' # Current sequence
			line_length = 0 # Line length (fasta format)
			charge = 0 # Current charge

			# Browse every aminoacid in protein
			for ii in range(window):

				# Count charges
				if lista[k+ii] == 'K' or lista[k+ii] == 'R':
					charge += 1
				if lista[k+ii] == 'E' or lista[k+ii] == 'D':
					charge -= 1
				if lista[k+ii] == 'H':
					charge += hcharge

				# If current line has more than 60 characters, go to next line
				if line_length >= 60:
					sequence += '\n'
					line_length = 0

				# Write current sequence
				sequence += lista[k+ii]
				# Update line length
				line_length += 1

			# Calculate ratio
			ratio = float(charge)/float(window)

			# If current ratio is higher than min ratio, or if current ratio is the highest in the protein 
			if ratio >= min_ratio or ratio >= max_ratio:
				# If charge is highest in the protein:
				if charge > max_charge:
					# Define current parameters as the highest ones
					max_ratio = ratio
					max_charge = charge
					best_sequence = sequence
					# Define position of the sequence in the protein
					startpos = k+1
					endpos = k+window
					
		
		window += 1

	# Write seuence in sequence archive
	fsequence.write(best_sequence)

	return str(max_charge) + ';' + str(startpos) + ';' + str(endpos) + ';' + str(max_ratio)

def main():
	if len(sys.argv) < 2:
		print ("Use: <program file name> <input file name>")
		sys.exit(0)

	# Open input file
	fin = open(sys.argv[1], 'r') 
	contents = fin.read()

	# Open output file
	fout = open(str.split(sys.argv[1], '.')[0]+'_charge.csv', 'w')
	fsequence = open(str.split(sys.argv[1], '.')[0]+'_sequence.fasta', 'w')

	# Define hcharge
	hcharge = -2
	while hcharge < -1:
		try:
			hcharge = float(input("Histidine charge: "))
		except:
			hcharge = -2
		if hcharge > 1:
			hcharge = -2

	# Define Field Size
	window = 0
	while window <= 0:
		try:
			window = int(input("Field Initial Size: "))
		except:
			window = 0

	# Define Min Ratio
	min_ratio = 0
	while min_ratio <= 0 or min_ratio > 1:
		try:
			min_ratio = float(input("Ratio of positive charges: "))
		except:
			min_ratio = 0


	# Split the proteins
	lista1 = str.split(contents, '>')

	# Split each line
	for i in range(1, len(lista1)):

		lista2 = lista1[i].replace('\r','')
		lista2 = str.split(lista2,'\n')

		# Concatenate only the aminoacid sequenceuence in 'lista3'
		lista3 = []
		for j in range(1, len(lista2)):

			# Split each character (aminoacid) in 'lista4'
			lista4 = list(lista2[j])

			lista3 += lista4

		# Write protein name and charge (from evaluate)
		fsequence.write('>'+lista2[0]+'\n')
		fout.write('>'+lista2[0]+";"+evaluate(lista3, hcharge, window, i, len(lista1), min_ratio, fsequence))
		if i < len(lista1):
			fout.write("\n")
			fsequence.write('\n')

	fin.close()
	fout.close()
	print '100% of proteins complete.                                    '
main()