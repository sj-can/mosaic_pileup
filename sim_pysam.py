import pysam, random

    #with open("target_reads", 'w') as ouf:
        #for item in read_names:
	#    ouf.write(item + '\n')


def pysam_read_identifier(bam_file_path):
        samfile = pysam.AlignmentFile(bam_file_path, 'rb')
	positions = {'A' : [], 'T' : [], 'C' : [], 'G' : [], 'N' : [] 'all_reads' : []}
        total = 0
	all_reads = []
        for pcolumn in samfile.pileup("MT", 3242, 3243):
	    for pread in pcolumn.pileups:
	        if pcolumn.pos == 3242:
		    read = pread.alignment.query_name
		    base = pread.alignment.query_sequence[int(pread.query_position)]
		    positions[base].append(read)
		    positions['all_reads'].append(pread)
		    total += 1
	positions['total'] = total
        samfile.close()
	print 'read_count ' + str(positions['total'])
        return positions

#takes target bam, pileup a subset and outputs a bam
def pysam_pileup(bam_file_path, mutant_name, mosaic, positions, ref=None, alt=None):
    total = positions['total']
    print 'total = ' + str(total)
    number_alts_needed = total/mosaic
    print 'mosacisism level required = ' + str(mosaic) + '%'
    print 'number ALT allele needed = ' + str(number_alts_needed)
    if alt:
        print 'number ALT present = ' + str(len(positions[alt]))
        mosaic_dict = { "number_alts_needed" : total/mosaic, "number_alts_have" : len(positions[alt]), "mosic_level" : mosaic}
        print mosaic_dict
    if number_alts_needed < len(positions[alt]):
	print 'you need a lower level of moscaicism!'
	if not ref and not alt:
	    print 'you must include reference and alternative alleles to continue'
	elif ref and alt:
	    mutator(bam_file_path, mutant_name, mosaic_dict, positions, ref, alt, higher=True)
    elif number_alts_needed > len(positions[alt]):
	print 'you need a higher lvel of moscaicism!'


def mutator(bam_file_path, mutant_name, mosaic_dict, positions, ref=None, alt=None, higher=False):
    print '------------------------'
    print 'MUTATOR ACTIVATED!'    
    print '------------------------'
    count = 0
	for pread in pcolumn.pileups:
	    if not pread.is_del and not pread.is_refskip and pcolumn.pos == 3242:
	    	count += 1
		if higher:
		    alt_to_ref_count = 0
		    #randomly iterate through every read, if its a 'G'/alt alter it to reference until alt_to_ref_count == 
		    
	    if pcolumn.pos != 3242:
		pass
    print count
    samfile.close()
    mutant_samfile.close()
        

#def pysam_pileup(bam_file_path):
#    samfile = pysam.AlignmentFile(bam_file_path, 'rb')
#    for pcolumn in samfile.pileup("MT", 3242, 3200):
#        for pread in pcolumn.pileups:
#            if not pread.is_del and not pread.is_refskip and pcolumn.pos == 3242:
#	        #print len(pread.alignment.get_reference_positions())
#		print pread.alignment.reference_start
#		if pread.alignment.reference_start == 3242:
#		    print pread.alignment.query_alignment_sequence
#		#print len(pread.alignment.query_alignment_sequence)
#
 #               #print pcolumn.n
  #              #print pread.query_position
   #             #print pcolumn.pos
    #            read = pread.alignment.query_name
     #           base = pread.alignment.query_sequence[int(pread.query_position)]
      #          #print type(base)
       #         #print read, base
		#print pread
#                #positions[base].append(read)
 #   #positions['total'] = total
 #   samfile.close()

#def pysam_iterator(input_bam, percentage_heteroplasmy):
#    mutant_name = str(percentage_heteroplasmy) + '_' + input_bam
#    samfile = pysam.AlignmentFile(bam_file_path, 'rb', template=samfile)
#    mutant_samfile = pysam.AlignmentFile(mutant_name, 'wb')
#    #iterate through all mitochondiral reads 
#    for read in samfile.fetch("MT"):
#	if samfile.pileup("MT", 3242, 3243):


if __name__ == '__main__':
    target_reads = pysam_read_identifier("v501_1169_FEMALE.realigned.bam")
    pysam_pileup("v501_1169_FEMALE.realigned.bam", "mutant_v501_1169_FEMALE.realigned.bam", 10, target_reads, ref='A', alt='G')
