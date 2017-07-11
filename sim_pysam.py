import pysam

#from the bam file convert regions "MT", 3242, 3243 into a sam file with header "wh"
#this includes 8213 reads in adition to the header

def bam_to_MTsam(bam_file, outfile):
    samfile = pysam.AlignmentFile(bam_file, 'rb')
    ouf = pysam.AlignmentFile(outfile, "w", template=samfile)
    for read in samfile.fetch("MT", 3242, 3243):
        ouf.write(read)
    return outfile

def bam_to_sam(bam_file, outfile):
    samfile = pysam.AlignmentFile(bam_file, 'rb')
    ouf = pysam.AlignmentFile(outfile, "wbh", template=samfile)
    for read in samfile.fetch():
	ouf.write(read)
    return outfile

def sam_to_bam(sam_file, outfile):
    samfile = pysam.AlignmentFile(sam_file, 'r')
    ouf = pysam.AlignmentFile(outfile, "wb", template=samfile)
    for read in samfile.fetch():
        ouf.write(read)
    return outfile

#use the MT bam file file to identify the alleles at the desired 3243 coordinate
#create a dictionary in the output positions dictionary that contains 'total'
#return positions dictionary
def pysam_allele_identifier(sam_file_path):
    samfile = pysam.AlignmentFile(sam_file_path, 'rb')
    positions = {'A' : [], 'T' : [], 'C' : [], 'G' : [], 'N' : []}
    total = 0
    for pcolumn in samfile.pileup("MT", 3242, 3243):
        for pread in pcolumn.pileups:
            if pcolumn.pos == 3242:
		#append read_id to base
		positions[pread.alignment.query_sequence[int(pread.query_position)]].append(pread.alignment.query_name)
		total += 1
    positions['total'] = total
    samfile.close()
    return positions

#function to handle logic of more or less mutants required and pass appropriate arguments to the 
#mutator function
def pysam_pileup(bam_file_path, mutant_name, mosaic, positions, ref=None, alt=None):
    total = positions['total']
    print 'total = ' + str(total)
    print 'mosacisism level required = ' + str(mosaic) + '%'
    number_alts_needed = (total/100) * mosaic
    print 'number ALT allele needed = ' + str(number_alts_needed)
    if alt:
        print 'number ALT present = ' + str(len(positions[alt]))
        mosaic_dict = { "number_alts_needed" : number_alts_needed, "number_alts_have" : len(positions[alt]), "mosic_level" : mosaic}
        print mosaic_dict
    if number_alts_needed < len(positions[alt]):
        print 'you need a lower level of moscaicism than is in this sample!'
        if not ref and not alt:
            print 'you must include reference and alternative alleles to continue'
        elif ref and alt:
            mutants = mutator(bam_file_path, mutant_name, mosaic_dict, positions, ref, alt, higher=True)
	    print 'mutator_activated'
    elif number_alts_needed > len(positions[alt]):
        print 'you need a higher level of moscaicism!'
    return mutants

#the mutator funtion :
# return a hash with an incrementing integer_unmutated_sequence_sequence_id as a key and the mutated sequence as the value
#{integer_AATTCCGG_read_id :  mutated_sequence]}
def mutator(bam_file_path, mutant_name, mosaic_dict, positions, ref=None, alt=None, higher=False):
    count = 0
    alt_to_ref_count = 0
    samfile = pysam.AlignmentFile(bam_file_path, 'rb')
    #{integer_AATTCCGG_read_id :  mutated_sequence]}
    sequences_to_substitute = {}
    duplicate_keys = {}
    if higher:
        decrease_by = mosaic_dict["number_alts_have"] - mosaic_dict["number_alts_needed"]
	print 'you need to decrease by ' + str(decrease_by)
    for pcolumn in samfile.pileup("MT", 3242, 3243):
        for pread in pcolumn.pileups:
            if not pread.is_del and not pread.is_refskip and pcolumn.pos == 3242:
	        count += 1
                if higher:
		    base = pread.alignment.query_sequence[int(pread.query_position)]
		    if base == alt and alt_to_ref_count < decrease_by:
			altered_list = list(pread.alignment.query_sequence)
			altered_list[int(pread.query_position)] = 'A'
			altered_sequence = ''.join(altered_list)
			#{integer_AATTCCGG_read_id :  mutated_sequence]}
			key = pread.alignment.query_sequence +'_'+ pread.alignment.query_name
			if key in sequences_to_substitute.keys():
			    duplicate_keys[key] = altered_sequence
			    alt_to_ref_count += 1
			elif key not in sequences_to_substitute.keys():
			    sequences_to_substitute[key] = altered_sequence
			    alt_to_ref_count += 1
    samfile.close()
    return [sequences_to_substitute, mutant_name, duplicate_keys]

#add mutant sequences to sam file adding the integer_unmutated_sequence_sequence_id to an array to ensure all are accounted for
def add_mutants_to_sam(sequences_to_substitute, mutant_name, MT_sam, duplicate_sequences_to_subs):
    match_array = []
    non_match_array = []
    other_match_array = []
    print len(duplicate_sequences_to_subs)
    with open(MT_sam, 'r') as inf, open(mutant_name, 'w') as ouf:
        for line in inf:
            spl = line.split('\t')
	    inf_key = spl[9] + '_' + spl[0]
            for k,v in sequences_to_substitute.iteritems():
		if inf_key == k and k not in match_array:
		    spl[9] = v
		    ouf.write('\t'.join(map(str,spl)))
		    match_array.append(k)
	    for k,v in duplicate_sequences_to_subs.iteritems():
		if inf_key == k and k not in other_match_array:
		    spl[9] = v
		    ouf.write('\t'.join(map(str,spl)))
		    other_match_array.append(k)
	    if inf_key not in sequences_to_substitute.keys():
		ouf.write(line)
		non_match_array.append(inf_key)
    print len(match_array)
    print len(other_match_array)
    print len(non_match_array)

#append the large headerless sam file to the mutated reads

#sort the sam file and convert to bam file

if __name__ == '__main__':
    #isolate m3243 reads
    #bam_to_MTsam = bam_to_MTsam("v501_1169_FEMALE.realigned.bam", "v501_1169_FEMALE.realigned.MT3243.sam")
    #generate sam with all reads 
    #bam_to_sam = bam_to_sam("v501_1169_FEMALE.realigned.bam", "v501_1169_FEMALE.realigned.sam")
    #externally remove target reads using bam_to_MTsam output and picard. Use FilterSamReads from picard using a headerless desired MT regions fil

    #target_reads = pysam_allele_identifier("v501_1169_FEMALE.realigned.bam")
    #mutate_reads = pysam_pileup("v501_1169_FEMALE.realigned.bam", "10_v501_1169_FEMALE.realigned.MT.sam", 10, target_reads, ref='A', alt='G')
    #add_mutants_to_sam(mutate_reads[0], mutate_reads[1], "v501_1169_FEMALE.realigned.MT3243.sam")

    #sam_to_bam("NO_MT3243_v501_1169_FEMALE.sam", "10_MT3243_v501_1169_FEMALE.bam")
    #mutated_bam = pysam_allele_identifier("10_MT3243_v501_1169_FEMALE.sorted.bam")
    #for k,v in mutated_bam.iteritems():
#	try:
#	    print k, len(v)
#	except:
#	    print k, v

    #this tests that the target reads have been removed 
    #sam_to_bam("NO_MT3243_v501_1169_FEMALE.sam", "NO_MT3243_v501_1169_FEMALE.bam")
    
    positions = pysam_allele_identifier("/mnt/Data4/working_directory/stuart/python-2-7-10/scripts/mitochondial_sensitivity/mutated_files/10_MT3243_v501_1169_FEMALE.sorted.bam")
    for k,v in positions.iteritems():
	try:
    	    print k, len(v)
	except:
	    print k, v

   #create samfile with only reads that map to m3243 (A)
   #create sam file with all reads and header
   #create sam file subrtacting the reads that map to m3243 using picard FilterSamReads and headerless (A) as readfile (B)
   #pull out target reads using a pileup
   #mutate the reads depending on required level of mosaicism - first n are mutated
   #append mutated read to sam file (B)
   #convert sam_file(B) back to bam format
   #sort mutated bam by coordinate


   #----------create 5% 
   #target_reads = pysam_allele_identifier("v501_1169_FEMALE.realigned.bam")
   #mutate_reads = pysam_pileup("v501_1169_FEMALE.realigned.bam", "05_v501_1169_FEMALE.realigned.MT.sam", 5, target_reads, ref='A', alt='G') 
   #add_mutants_to_sam(mutate_reads[0], mutate_reads[1], "v501_1169_FEMALE.realigned.MT3243.sam", mutate_reads[2])
   #sam_to_bam("05_v501_1169_FEMALE.realigned.sam", "05_v501_1169_FEMALE.realigned.bam")

   #check for duplicate
   #bam_to_MTsam("v501_1169_FEMALE.realigned.bam", "outfile")
