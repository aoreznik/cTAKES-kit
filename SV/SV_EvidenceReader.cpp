
using namespace BamTools;

using namespace std;

int  SV_EvidenceReader:: counter = 1;
map<int,string> SV_EvidenceReader:: sample_names;
map<int,string> SV_EvidenceReader:: ev_types;

SV_EvidenceReader::
~SV_EvidenceReader()
{
}

string
SV_EvidenceReader::
check_params()
{
	string msg = "";

	return msg;
}

bool
SV_EvidenceReader::
{
	return false;
}

void
SV_EvidenceReader::
initialize()
{
}

void
SV_EvidenceReader::
set_statics()
{
}

void
SV_EvidenceReader::
unset_statics()
{
}

void
SV_EvidenceReader::
process_input( BamAlignment &_bam,
			   RefVector &_refs,
			   string header,
{
}

void
SV_EvidenceReader::
{
}

void
SV_EvidenceReader::
process_input_chr(string chr,
{
}

void
SV_EvidenceReader::
process_input_chr_pos(string chr,
					  CHR_POS pos,
{
}

void
SV_EvidenceReader::
process_input( BamAlignment &_bam,
			   RefVector &_refs,
			   BamWriter &inter_chrom_reads,
{
	abort();
}

void
SV_EvidenceReader::
process_input( BamAlignment &_bam,
			   RefVector &_refs,
{
	abort();
}


void
SV_EvidenceReader::
terminate()
{
}

string
SV_EvidenceReader::
get_curr_chr()
{
	return "";
}

CHR_POS
SV_EvidenceReader::
get_curr_pos()
{
	CHR_POS x = 0;
	return x;
}

bool
SV_EvidenceReader::
has_next()
{
	return false;
}

string
SV_EvidenceReader::
get_source_file_name()
{
	return "Error";
}

int32_t
SV_EvidenceReader::
get_curr_primary_refid()
{
	cerr << "Error reaching SV_EvidenceReader:: get_curr_primary_refid "
		<<endl;
	abort();
	int32_t a = 1;
	return a;
}

int32_t
SV_EvidenceReader::
get_curr_secondary_refid()
{
	cerr << "Error reaching SV_EvidenceReader:: get_curr_secondary_refid "
		<<endl;
	abort();
	int32_t a = 1;
	return a;
}


string
SV_EvidenceReader::
get_curr_primary_chr()
{
	cerr << "Error reaching SV_EvidenceReader:: get_curr_primary_chr "
		<<endl;
	abort();
	string a = "";
	return a;
}

string
SV_EvidenceReader::
get_curr_secondary_chr()
{
	cerr << "Error reaching SV_EvidenceReader:: get_curr_secondary_chr "
		<<endl;
	abort();
	string a = "";
	return a;
}

CHR_POS
SV_EvidenceReader::
get_curr_primary_pos()
{
	abort();
	CHR_POS a = 0;
	return a;
}

CHR_POS
SV_EvidenceReader::
get_curr_secondary_pos()
{
	abort();
	CHR_POS a = 0;
	return a;
}

void
SV_EvidenceReader::
process_input_chr_pos(string primary_chr,
					  string secondary_chr,
					  CHR_POS pos,
{
	cerr << "Error reaching SV_EvidenceReader:: process_input_chr_pos "
		<<endl;
	abort();
}
