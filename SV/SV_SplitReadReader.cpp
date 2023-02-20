

SV_SplitReadReader::
~SV_SplitReadReader()
{
}

SV_SplitReadReader::
SV_SplitReadReader()
{
	bam_file = "";
	sample_name = "";
	back_distance = 0;
	min_non_overlap = 25;
	weight = -1;
	min_mapping_threshold = 0;
	ev_id = SV_EvidenceReader::counter;
	SV_EvidenceReader::counter = SV_EvidenceReader::counter + 1;
	SV_EvidenceReader::sample_names[ev_id] = sample_name;
	SV_EvidenceReader::ev_types[ev_id] = "SR";
	min_clip = 20;
}

string
SV_SplitReadReader::
check_params()
{
	string msg = "";

	if (bam_file.compare("") == 0)
		msg.append("bam_file ");
	if (sample_name.compare("") == 0)
		msg.append("id ");
	if (back_distance == 0)
		msg.append("back_distance ");
	if (weight == 0)
		msg.append("weight ");

	return msg;
}

bool
SV_SplitReadReader::
{
	if ( strcmp("bam_file", param) == 0 )
		bam_file = val;
	else if ( strcmp("id", param) == 0 ) {
		sample_name = val;
		SV_EvidenceReader::sample_names[ev_id] = sample_name;
		SV_EvidenceReader::ev_types[ev_id] = "SR";
	}
	else if ( strcmp("min_non_overlap", param) == 0 )
		min_non_overlap = atoi(val);
	else if ( strcmp("back_distance", param) == 0 )
		back_distance = atoi(val);
	else if ( strcmp("weight", param) == 0 )
		weight = atoi(val);
	else if ( strcmp("min_mapping_threshold", param) == 0 )
		min_mapping_threshold = atoi(val);
	else if ( strcmp("min_clip", param) == 0 )
		min_clip = atoi(val);
	else if ( strcmp("read_group", param) == 0 )
		read_group.push_back(val);
	else 
		return false;

	return true;
}

void
SV_SplitReadReader::
set_statics()
{
}

void
SV_SplitReadReader::
unset_statics()
{
}

void
SV_SplitReadReader::
initialize()
{
}

string
SV_SplitReadReader::
get_curr_chr()
{
	return refs.at(bam.RefID).RefName;
}

CHR_POS
SV_SplitReadReader::
get_curr_pos()
{
	return bam.Position;
}

void
SV_SplitReadReader::
{
	while (reader.GetNextAlignment(bam)) 
		SV_SplitRead::process_split(bam,
									refs,
									reader.GetHeader().ToString(),
									mapped_splits,
									r_bin,
									weight,
									ev_id,
									this);
}

void
SV_SplitReadReader::
process_input( BamAlignment &_bam,
			   RefVector &_refs,
{
	SV_SplitRead::process_split(_bam,
								_refs,
								mapped_splits,
								r_bin,
								weight,
								ev_id,
								this);
}

void
SV_SplitReadReader::
process_input( BamAlignment &_bam,
			   RefVector &_refs,
			   BamWriter &inter_chrom_reads,
{
	SV_SplitRead::process_intra_chrom_split(_bam,
											_refs,
											inter_chrom_reads,
											mapped_splits,
											r_bin,
											weight,
											ev_id,
											this);
}

void
SV_SplitReadReader::
process_input_chr(string chr,
{
	cerr << "SplitRead:" << chr << endl;
	while ( have_next_alignment &&
			( chr.compare( refs.at(bam.RefID).RefName) == 0 ) ) {

		SV_SplitRead::process_split(bam,
									refs,
									reader.GetHeader().ToString(),
									mapped_splits,
									r_bin,
									weight,
									ev_id,
									this);

		have_next_alignment = reader.GetNextAlignment(bam);
		if ( bam.RefID < 0 )
			have_next_alignment = false;
	}
}

void
SV_SplitReadReader::
process_input_chr_pos(string chr,
					  CHR_POS pos,
{
	cerr << "SplitRead:" << chr << endl;
	while ( have_next_alignment &&
			( chr.compare( refs.at(bam.RefID).RefName) == 0 ) &&
			( bam.Position < pos ) ) {
		SV_SplitRead::process_split(bam,
									refs,
									reader.GetHeader().ToString(),
									mapped_splits,
									r_bin,
									weight,
									ev_id,
									this);

		have_next_alignment = reader.GetNextAlignment(bam);
		if ( bam.RefID < 0 )
			have_next_alignment = false;
	}
}

void 
SV_SplitReadReader::
terminate()
{
	reader.Close();
}

bool
SV_SplitReadReader::
has_next()
{
	return have_next_alignment;
}

string
SV_SplitReadReader::
get_source_file_name()
{
	return bam_file;
}
