
import pysam, subprocess, os

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
		    v = bam_file
		else:
		    matching_ids[k] = v[0]
	if type(v) == str:
	    matching_ids[k] = [v]
    return matching_ids

#takes key, value as {tngs_id : ['one_bam_file_path']
def pysam_pileup(tngs_id, bam_file_path):
    positions = {'A' : [], 'T' : [], 'C' : [], 'G' : [], 'N' : []}
    samfile = pysam.AlignmentFile(bam_file_path, 'rb')
    total = 0
    for pcolumn in samfile.pileup("MT", 3242, 3243):
        for pread in pcolumn.pileups:
            if not pread.is_del and not pread.is_refskip and pcolumn.pos == 3242:
		total += 1
        	#print pcolumn.n
        	print pread.query_position
    	        #print pcolumn.pos
	    	read = pread.alignment.query_name
	    	base = pread.alignment.query_sequence[int(pread.query_position)]
		print type(base)
		print read, base
		positions[base].append(read)
	        #print base + '\t' + read
    positions['total'] = total
    samfile.close()
#    for k, v in positions.iteritems():
#	print k
#	try:
#	    print len(v)
#	except:
#	    print v
    return positions

def starting_file(filename):
    new_name = filename.rsplit(".", 1)[0] + "_ac"
    if not os.path.isfile(new_name):
        with open(new_name, 'w') as f:
            f.write('tngs_id' + '\t' + 'A' + '\t' + 'T' + '\t' + 'C' + '\t' 'G' + '\t' + 'Total' + '\t' + 'bam_file_path' + '\n')
    return new_name

def append_iteration_to_file(file_name, tngs_id, bam_path, d):
    with open(file_name, 'a') as ouf:
	ouf.write(tngs_id +'\t'+ str(len(d['A'])) +'\t'+ str(len(d['T'])) +'\t'+ str(len(d['C'])) +'\t'+ str(len(d['G'])) +'\t'+ str(d['total']) +'\t' + bam_path + '\n')

if __name__ == "__main__":
    b = 'sample_lists/'
    #tngs_id_list = [b +"tngs_id_positive_controls.txt", b + "undiagnosed_controls.txt", b + "assumed_negative_control.txt"]
    test_list = ["tngs_id_positive_controls.txt", "assumed_negative_control.txt"]
    for item in test_list:
	tngs_ids = get_tngs_ids(item)
        bam_dict = search_bam_list_file(tngs_ids, "tngs_id_bam_paths")
	one_bam = only_one_bam(bam_dict)
	out_file = starting_file(item)
	for k,v in one_bam.iteritems():
	    this_bam_pile = pysam_pileup(k, v[0])
	    append_iteration_to_file(out_file, k, str(v[0]), this_bam_pile)
	    
	    
		
