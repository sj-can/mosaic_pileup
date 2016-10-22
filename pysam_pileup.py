import pysam, subprocess

def get_tngs_ids(input_file):
        tngs_ids = []
        with open(input_file, 'r') as inf:
            for line in inf:
                altered_line = line.strip(' \n\r')
                tngs_ids.append(altered_line)
        return tngs_ids

def search_bam_list_file(tngs_ids, bam_find_shell):
    matching_ids = {}
    non_matching_ids = []
    bam_file_paths = []
    command = ["bash", bam_find_shell]
    process = subprocess.Popen(command, stdout=subprocess.PIPE)
    bam_file_paths = process.communicate()[0]
    all_bams = bam_file_paths.split('\n')
    for item in tngs_ids:
	match_count = 0
        for bam_file in all_bams:
            if item in bam_file and 'realigned' in bam_file:
                match_count += 1
                if match_count == 1:
                     #create_new_dictionary entry with id as key
                     matching_ids[item] = [bam_file]
                if match_count > 1:
                     #add to existing key
                     matching_ids[item].append(bam_file)
        if match_count == 0:
            non_matching_ids.append(item)
#    with open("summary", 'a') as s:
#	s.write('tngs_id_count = ' + str(len(tngs_ids)) + '\n')
#	s.write('matching_count = ' + str(len(matching_ids)) + '\n')
#	s.write('non_mathcing_count = ' + str(len(non_matching_ids)) + '\n')
#	s.write(str(non_matching_ids) + '\n')
    return matching_ids

def only_one_bam(matching_ids):
    for k,v in matching_ids.iteritems():
	if len(v) > 1:
	    for bam_file in v:
		if 'test' in bam_file:
		    v.remove(bam_file)
		    matching_ids[k] = v
		elif 'rerun' in bam_file:
		    v = [bam_file]
		else:
		    matching_ids[k] = v[0]
    return matching_ids

#takes key, value as {tngs_id : ['one_bam_file_path']
def pysam_pileup(tngs_id, bam_file_path):
    positions = {'A' : [], 'T' : [], 'C' : [], 'G' : []}
    samfile = pysam.AlignmentFile(bam_file_path, 'rb')
    total = 0
    for pcolumn in samfile.pileup("MT", 3242, 3243):
        for pread in pcolumn.pileups:
            #if not pread.is_del and not pread.is_refskip and pcolumn.pos == 3242:
	    if pcolumn.pos == 3242:
		total += 1
        	#print pcolumn.n
        	#print pread.query_position
    	        #print pcolumn.pos
	    	read = pread.alignment.query_name
	    	base = pread.alignment.query_sequence[pread.query_position]
		#print read, base
		positions[base].append(read)
	        #print base + '\t' + read
    positions['total'] = total
    samfile.close()
    for k, v in positions.iteritems():
	print k
	try:
	    print len(v)
	except:
	    print v
    print

def starting_file(filename):
    with open(filename, 'w') as f:
        f.write('tngs_id' + '\t' + 'A' + '\t' + 'T' + '\t' + 'C' + '\t' 'G' + '\t' + 'Total' + '\t' + 'bam_file_path')
    return filename

#def append_iteration_to_file(file_name, position_dict):


if __name__ == "__main__":
    b = 'sample_lists/'
    #tngs_id_list = [b +"tngs_id_positive_controls.txt", b + "undiagnosed_controls.txt", b + "assumed_negative_control.txt"]
    test_list = [tngs_id_positive_controls.txt]
    for item in test_list:
	tngs_ids = get_tngs_ids(item)
        bam_dict = search_bam_list_file(tngs_ids, "tngs_id_bam_paths")
	one_bam = only_one_bam(bam_dict)
	starting_file(item)
#	for k,v in one_bam.iteritems():
#	    if k == 'v501_1169':
#		print v[0]
#		pysam_pileup(str(k), str(v[0]))	    
