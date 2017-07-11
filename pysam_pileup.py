import pysam, subprocess, os

def get_tngs_ids(input_file):
        tngs_ids = []
        with open(input_file, 'r') as inf:
            for line in inf:
                altered_line = line.strip(' \n\r')
                tngs_ids.append(altered_line)
        return tngs_ids

def tngs_id_with_bam_path(input_file):
    bam_dict = {}
    with open(input_file, 'r') as inf:
        for line in inf:
	    altered_line = line.strip(' \n\r')
	    stripped_line = altered_line.split('\t')
	    bam_dict[stripped_line[0]] = [stripped_line[1]]
    return bam_dict

def search_data1(tngs_ids, bam_find_shell):
    matching_ids = {}
    non_matching_ids = []
    bam_file_paths = []
    command = ["bash", bam_find_shell]
    process = subprocess.Popen(command, stdout=subprocess.PIPE)
    bam_file_paths = process.communicate()[0]
    all_bams = bam_file_paths.split('\n')
    for item in tngs_ids:
        match_count = 0
	sample_code = item[-3:]
	for bam_file in all_bams:
	    if sample_code in bam_file and 'realigned' in bam_file and 'P5' in bam_file:
	        match_count += 1
		if match_count == 1:
		    matching_ids[item] = [bam_file]
		if match_count > 1:
                     #add to existing key
                     matching_ids[item].append(bam_file)
	if match_count == 0:
	    non_matching_ids.append(item)
    return [matching_ids, non_matching_ids]

def search_bam_list_file(tngs_ids, bam_find_shell):
    print len(tngs_ids)
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
		    v = bam_file
		elif 'run2' in bam_file:
		    matching_ids[k] = bam_file
		else:
		    matching_ids[k] = v[0]
    for key, value in matching_ids.iteritems():
        if type(value) == str:
	    matching_ids[key] = [value]
    return matching_ids

#takes key, value as {tngs_id : ['one_bam_file_path']}
def pysam_pileup(tngs_id, bam_file_path):
    positions = {'A' : [], 'T' : [], 'C' : [], 'G' : [], 'N' : []}
    samfile = pysam.AlignmentFile(bam_file_path, 'rb')
    total = 0
    for pcolumn in samfile.pileup("MT", 3140, 3345):
        for pread in pcolumn.pileups:
            if not pread.is_del and not pread.is_refskip and pcolumn.pos == 3242:
		total += 1
        	#print pcolumn.n
        	#print pread.query_position
    	        #print pcolumn.pos
	    	read = pread.alignment.query_name
	    	base = pread.alignment.query_sequence[int(pread.query_position)]
		#print type(base)
		#print read, base
		positions[base].append(read)
    positions['total'] = total
    samfile.close()
    return positions

def starting_file(filename):
    new_name = filename.rsplit(".", 1)[0] + "r02_undiagnosed_easy_mathces_ac"
    if not os.path.isfile(new_name):
        with open(new_name, 'w') as f:
            f.write('tngs_id' + '\t' + 'A' + '\t' + 'T' + '\t' + 'C' + '\t' 'G' + '\t' + 'Total' + '\t' + 'bam_file_path' + '\n')
    return new_name

def append_iteration_to_file(file_name, tngs_id, bam_path, d):
    with open(file_name, 'a') as ouf:
	ouf.write(tngs_id +'\t'+ str(len(d['A'])) +'\t'+ str(len(d['T'])) +'\t'+ str(len(d['C'])) +'\t'+ str(len(d['G'])) +'\t'+ str(d['total']) +'\t' + bam_path + '\n')

def output_list(list, file_name):
    with open(file_name, 'w') as f:
        for item in list:
	    f.write(item + '\n')

def check_for_dir_name(non_matching_ids):
    dir_name = ["001", "048", "049", "096", "097", "132", "144", "145", "216", "217", "240", "241", "336", "337", "384", "385", "432", "433", "480", "481", "504"]
    print len(non_matching_ids)
    with open("r02_in_batch_dir_name", 'w') as f:
        for id in non_matching_ids:
            for item in dir_name:
	        if item in id:
		    non_matching_ids.remove(id)
		    f.write(id + '\n')
    return non_matching_ids

if __name__ == "__main__":
    b = 'sample_lists/'
    test_list = ["undiagnosed_sc.txt"]
    for item in test_list:
	tngs_ids = get_tngs_ids(item)
	not_matching = search_bam_list_file(tngs_ids, "tngs_id_bam_paths")
        print len(not_matching)

   #use this to search data 1
	#rm_dir_names = check_for_dir_name(not_matching)
	#match_non_match_list = search_data1(rm_dir_names, "data1_bam_paths")
	#print len(match_non_match_list[0])
	#output_list(match_non_match_list[1], "P5_positives_latest_matches")
        #one_bam_p5 = only_one_bam(match_non_match_list[0])
        #out_file = starting_file(item)
	#for k,v in one_bam_p5.iteritems():
	#    print k,v 
	#    this_bam_pile = pysam_pileup(k, v[0])
	#    append_iteration_to_file(out_file, k, str(v[0]), this_bam_pile)

   #use this if there are no exceptions
        bam_dict = search_bam_list_file(tngs_ids, "tngs_id_bam_paths")
	#print bam_dict
	print len(bam_dict)
	one_bam = only_one_bam(bam_dict)
	out_file = starting_file(item)
	for k,v in one_bam.iteritems():
	    print k, v
	    this_bam_pile = pysam_pileup(k, v[0])
	    append_iteration_to_file(out_file, k, str(v[0]), this_bam_pile)
	   
    #use this if bam file already found
	#tngs_and_bam_dict = tngs_id_with_bam_path(item)
	#out_file = starting_file(item)
	#for k,v in tngs_and_bam_dict.iteritems():
	    #print k,v 
	    #this_bam_pile = pysam_pileup(k, v[0])
	    #append_iteration_to_file(out_file, k, str(v[0]), this_bam_pile)
	
		
