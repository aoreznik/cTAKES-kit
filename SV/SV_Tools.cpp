
using namespace BamTools;

const unsigned int SORT_DEFAULT_MAX_BUFFER_COUNT  = 500000; 
const unsigned int SORT_DEFAULT_MAX_BUFFER_MEMORY = 1024; 




using namespace std;

{
			long ai, bi;
			ai = strtol(pa, &pa, 10);
			bi = strtol(pb, &pb, 10);
			if (ai != bi) return ai<bi? -1 : ai>bi? 1 : 0;
		} else {
			++pa; ++pb;
		}
	}
		return (pa-a) < (pb-b)? -1 : (pa-a) > (pb-b)? 1 : 0;

}

int read_histo_file(string file_name,
{
	int file_len = 0;
	char line[LINE_MAX];

	if (file == NULL) {
		fprintf(stderr, "Could not open %s\n", file_name.c_str());
	}

	while ( fgets(line, LINE_MAX, file) )
		++file_len;

	fclose(file);

	file = fopen(file_name.c_str(), "r");
	int i = 0;
	bool first = true;
	while ( fgets(line, LINE_MAX, file) ) {
		unsigned int tmp = atoi(strtok(line, "\t"));
		if (first == true) {
			first = false;
		}
		++i;
	}


	fclose(file);

	return file_len;
}

int read_distro_file(string file_name,
{
	int file_len = 0;
	char line[LINE_MAX];

	if (file == NULL) {
		fprintf(stderr, "Could not open %s\n", file_name.c_str());
	}

	while ( fgets(line, LINE_MAX, file) )
		++file_len;

	fclose(file);


	file = fopen(file_name.c_str(), "r");
	int i = 0;
	bool first = true;
	int last;
	while ( fgets(line, LINE_MAX, file) ) {
		int v = atof(strtok(line, "\t"));
		if (first == true) {
			first = false;
		}
		last = v;
                double d = atof(strtok(NULL, "\t"));
		++i;
	}


	fclose(file);

	return file_len;
}

struct inter_chrom_sort
{
    bool operator()( const BamAlignment& l, const BamAlignment& r ) const {

	
	int32_t l_primary_ref = l.RefID,
			l_primary_pos = l.Position,
			l_secondary_pos = l.MatePosition,
			l_secondary_ref = l.MateRefID;
	
	int32_t r_primary_ref = r.RefID,
			r_primary_pos = r.Position,
			r_secondary_ref = r.MateRefID,
			r_secondary_pos = r.MatePosition;

	if (l.RefID > l.MateRefID) {
		l_primary_ref = l.MateRefID;
		l_primary_pos = l.MatePosition;
		l_secondary_ref = l.RefID;
		l_secondary_pos = l.Position;
	}

	if (r.RefID > r.MateRefID) {
		r_primary_ref = r.MateRefID;
		r_primary_pos = r.MatePosition;
		r_secondary_ref = r.RefID;
		r_secondary_pos = r.Position;
	}


        if (l_primary_ref < r_primary_ref) {
            return true;
        } else if (l_primary_ref > r_primary_ref) {
            return false;
        } else if (l_primary_ref == r_primary_ref) {
            if (l_secondary_ref < r_secondary_ref) {
                return true;
            } else if (l_secondary_ref > r_secondary_ref) {
                return false;
            } else if (l_secondary_ref == r_secondary_ref) {
                if (l_primary_pos < r_primary_pos) {
                    return true;
                } else if (l_primary_pos > r_primary_pos) {
                    return false;
                } else if (l_primary_pos == r_primary_pos) {
                    if (l_secondary_pos < r_secondary_pos) {
                        return true;
                    } else if (l_secondary_pos > r_secondary_pos) {
                        return false;
                    } else if (l_secondary_pos == r_secondary_pos) {
                        return true;
                    }
                }
            }
        }

        exit(1); 

	return  ( (l_primary_ref < r_primary_ref) ||

			  ((l_primary_ref == r_primary_ref) && 
					(l_secondary_ref < r_secondary_ref)) ||

			  ((l_primary_ref == r_primary_ref) && 
					(l_secondary_ref == r_secondary_ref) &&
					(l_primary_pos < r_primary_pos)) ||

			  ((l_primary_ref == r_primary_ref) && 
					(l_secondary_ref == r_secondary_ref) &&
					(l_primary_pos == r_primary_pos) && 
					(l_secondary_pos == r_secondary_pos)) );
    }
};

struct inter_chrom_rev_sort
{
    bool operator()( const BamAlignment& l, const BamAlignment& r ) const {

	
	int32_t l_primary_ref = l.RefID,
			l_primary_pos = l.Position,
			l_secondary_pos = l.MatePosition,
			l_secondary_ref = l.MateRefID;
	
	int32_t r_primary_ref = r.RefID,
			r_primary_pos = r.Position,
			r_secondary_ref = r.MateRefID,
			r_secondary_pos = r.MatePosition;

	if (l.RefID > l.MateRefID) {
		l_primary_ref = l.MateRefID;
		l_primary_pos = l.MatePosition;
		l_secondary_ref = l.RefID;
		l_secondary_pos = l.Position;
	}

	if (r.RefID > r.MateRefID) {
		r_primary_ref = r.MateRefID;
		r_primary_pos = r.MatePosition;
		r_secondary_ref = r.RefID;
		r_secondary_pos = r.Position;
	}


        if (l_primary_ref < r_primary_ref) {
            return false;
        } else if (l_primary_ref > r_primary_ref) {
            return true;
        } else if (l_primary_ref == r_primary_ref) {
            if (l_secondary_ref < r_secondary_ref) {
                return false;
            } else if (l_secondary_ref > r_secondary_ref) {
                return true;
            } else if (l_secondary_ref == r_secondary_ref) {
                if (l_primary_pos < r_primary_pos) {
                    return false;
                } else if (l_primary_pos > r_primary_pos) {
                    return true;
                } else if (l_primary_pos == r_primary_pos) {
                    if (l_secondary_pos < r_secondary_pos) {
                        return false;
                    } else if (l_secondary_pos > r_secondary_pos) {
                        return true;
                    } else if (l_secondary_pos == r_secondary_pos) {
                        return false;
                    }
                }
            }
        }

        exit(1); 

	return  ( (l_primary_ref < r_primary_ref) ||

			  ((l_primary_ref == r_primary_ref) && 
					(l_secondary_ref < r_secondary_ref)) ||

			  ((l_primary_ref == r_primary_ref) && 
					(l_secondary_ref == r_secondary_ref) &&
					(l_primary_pos < r_primary_pos)) ||

			  ((l_primary_ref == r_primary_ref) && 
					(l_secondary_ref == r_secondary_ref) &&
					(l_primary_pos == r_primary_pos) && 
					(l_secondary_pos == r_secondary_pos)) );
    }
};

bool sort_inter_chrom_bam(string in_file_name,
						  string out_file_name)
{
    BamReader reader;
    if ( !reader.Open(in_file_name) ) {
        cerr << "sort ERROR: could not open " << 
			in_file_name << " for reading... Aborting." << endl;
        return false;
    }

    SamHeader header = reader.GetHeader();
    if ( !header.HasVersion() )
        header.Version = Constants::SAM_CURRENT_VERSION;

    string header_text = header.ToString();
    RefVector ref = reader.GetReferenceData();

    BamAlignment al;
    vector<BamAlignment> buffer;
    bool bufferFull = false;

	
    int buff_count = 0;
    while ( reader.GetNextAlignment(al)) {

        bufferFull = ( buffer.size() >= SORT_DEFAULT_MAX_BUFFER_COUNT );

        if ( !bufferFull )
            buffer.push_back(al);
        else {
            create_sorted_temp_file(buffer,
                                    out_file_name,
                                    buff_count,
                                    header_text,
                                    ref);
                                    ++buff_count;
            buffer.push_back(al);
        }
    }

    if ( !buffer.empty() ) {
        create_sorted_temp_file(buffer,
                                out_file_name,
                                buff_count,
                                header_text,
                                ref);

        ++buff_count;
    }

    reader.Close();

    return merge_sorted_files(out_file_name, buff_count, header_text, ref);

	for (int i = 0; i < buff_count; ++i) {
    	stringstream temp_name;
    	temp_name << out_file_name << i;
	}
}

bool create_sorted_temp_file(vector<BamAlignment>& buffer,
							 string out_file_name,
							 int num_runs,
							 string header_text,
    						 RefVector &ref)
{
 
    stable_sort(buffer.begin(), buffer.end(), inter_chrom_sort());
  
    stringstream temp_name;
    temp_name << out_file_name << num_runs;

    bool success = write_temp_file(buffer, temp_name.str(), header_text, ref );
    
    buffer.clear();
    
    return success;
}

bool write_temp_file(vector<BamAlignment>& buffer,
					 string temp_file_name,
					 string header_text,
    				 RefVector &ref)
{
    BamWriter temp_writer;
    if ( !temp_writer.Open(temp_file_name, header_text, ref) ) {
        cerr << "sort ERROR: could not open " << temp_file_name
             << " for writing." << endl;
        return false;
    }

    vector<BamAlignment>::const_iterator buffIter = buffer.begin();
    vector<BamAlignment>::const_iterator buffEnd  = buffer.end();
    for ( ; buffIter != buffEnd; ++buffIter )  {
        temp_writer.SaveAlignment(al);
    }

    temp_writer.Close();
    return true;
}

bool merge_sorted_files(string out_file_name,
						int buff_count,
						string header_text,
						RefVector &ref)
{

    priority_queue< BamAlignment, vector<BamAlignment>, inter_chrom_rev_sort > q;

    for (int i = 0; i < buff_count; ++i) {
        stringstream temp_name;
        temp_name << out_file_name << i;


        if ( !reader->Open(temp_name.str()) ) {
            cerr << "sort ERROR: could not open " << 
                    temp_name.str() << " for reading... Aborting." << endl;
            return false;
        }

        bam_readers[temp_name.str()] = reader;
        BamAlignment al;
        if (reader->GetNextAlignment(al))
        q.push(al);
    }

    BamWriter merged_writer;
    if ( !merged_writer.Open(out_file_name, header_text, ref) ) {
        cerr << "sort ERROR: could not open " << out_file_name
                << " for writing." << endl;
        return false;
    }


    while (!q.empty()) {
        BamAlignment al = q.top();
        q.pop();
        merged_writer.SaveAlignment(al);


        BamAlignment new_al;

        if (reader->GetNextAlignment(new_al))
            q.push(new_al);
    }

    merged_writer.Close();

    for (it = bam_readers.begin(); it != bam_readers.end(); ++it) {
        reader->Close();
        delete reader;
        remove(it->first.c_str());
    }

    return true;
}

{
    log_space sum = -INFINITY;

    for (CHR_POS i = 0; i < size; ++ i)
        sum = ls_add(sum, o[i]);

    for (CHR_POS i = 0; i < size; ++ i)
        r[i] = ls_divide(o[i],sum);
}

void parse_exclude_file(string exclude_bed_file,
                        UCSCBins<int> &exclude_regions)
{
    bed_file->Open();
    BED bed;
    while (bed_file->GetNextBed(bed, false)) {
        if (bed_file->_status == BED_VALID) {
            exclude_regions.add(bed.chrom,
                                bed.start,
                                bed.end,
                                '+',
                                1);
        }
    }
    bed_file->Close();
}

uint32_t
count_clipped(vector< CigarOp > cigar_data)
{
    uint32_t match_count = 0;
    vector< CigarOp >::iterator i;
    for (i = cigar_data.begin(); i != cigar_data.end(); ++i) {
        if ( (i->Type == 'M') || (i->Type == 'S') ) {
            uint32_t l = i->Length;
            match_count += l;
        }
    }

    return match_count;
}
