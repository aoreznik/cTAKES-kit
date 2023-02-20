
using namespace BamTools;



using namespace std;

double  SV_Pair:: insert_Z = 0;

SV_Pair::
SV_Pair(const BamAlignment &bam_a,
        const BamAlignment &bam_b,
        const RefVector &refs,
        int _weight,
        int _ev_id,
{
    reader = _reader;

    if ( bam_a.MapQuality < bam_b.MapQuality )
        min_mapping_quality = bam_a.MapQuality;
    else
        min_mapping_quality = bam_b.MapQuality;

    read_id = bam_a.Name;

    struct interval tmp_a, tmp_b;
    tmp_a.start = bam_a.Position;
    tmp_a.end = bam_a.GetEndPosition(false, false) - 1;
    tmp_a.chr = refs.at(bam_a.RefID).RefName;

    if ( bam_a.IsReverseStrand() == true )
        tmp_a.strand = '-';
    else
        tmp_a.strand = '+';

    tmp_b.start = bam_b.Position;
    tmp_b.end = bam_b.GetEndPosition(false, false) - 1;
    tmp_b.chr = refs.at(bam_b.RefID).RefName;

    if ( bam_b.IsReverseStrand() == true )
        tmp_b.strand = '-';
    else
        tmp_b.strand = '+';

    if ( bam_a.RefID < bam_b.RefID ) {
        read_l = tmp_a;
        read_r = tmp_b;
    } else if ( bam_a.RefID > bam_b.RefID) {
        read_l = tmp_b;
        read_r = tmp_a;
        if (tmp_a.start > tmp_b.start) {
            read_l = tmp_b;
            read_r = tmp_a;
        } else {
            read_l = tmp_a;
            read_r = tmp_b;
        }
    }

    weight = _weight;
    ev_id = _ev_id;
}

SV_Pair::
get_bp_interval_probability(char strand,
                            int distro_size,
{
    int size = distro_size;
    int j;
    for (j = 0; j < size; ++j) {
        if (strand == '+')
            tmp_p[j] = get_ls(distro[j]);
        else
            tmp_p[(size - 1) - j] = get_ls(distro[j]);
    }

    return tmp_p;
}

void
SV_Pair::
                          unsigned int back_distance,
                          unsigned int distro_size)
{
    i->i.chr = target_interval->chr;
    i->i.strand = target_interval->strand;
    if ( i->i.strand == '+' ) {

        if (back_distance > target_interval->end) {
            i->i.start = 0;
            i->i.start_clip = target_interval->end;
        }
        else {
            i->i.start = target_interval->end - back_distance;
            i->i.start_clip = 0;
        }

        i->i.end = i->i.start + distro_size - 1;

    } else {
        i->i.end = target_interval->start + back_distance;

        if (distro_size > i->i.end) {
            i->i.start = 0;
            i->i.start_clip = i->i.end;
        }
        else {
            i->i.start = i->i.end - distro_size + 1;
            i->i.start_clip = 0;
        }
    }
}

SV_Pair::
get_bp()
{

    set_bp_interval_start_end(&(new_bp->interval_l),
                              &read_l,
                              &read_r,
                              reader->back_distance,
                              reader->distro_size);
    set_bp_interval_start_end(&(new_bp->interval_r),
                              &read_r,
                              &read_l,
                              reader->back_distance,
                              reader->distro_size);

    if ( (new_bp->interval_l.i.chr.compare(read_r.chr) == 0 ) &&
            (new_bp->interval_l.i.end >= read_r.start) ) {
        new_bp->interval_l.i.end_clip =
            new_bp->interval_l.i.end - read_r.start + 1;
        new_bp->interval_l.i.end = read_r.start - 1;
    }

    if ( (new_bp->interval_r.i.chr.compare(read_l.chr) == 0 ) &&
            (new_bp->interval_r.i.start <= read_l.end) ) {
        new_bp->interval_r.i.start_clip =
            read_l.end - new_bp->interval_r.i.start + 1;
        new_bp->interval_r.i.start = read_l.end + 1;
    }

    new_bp->interval_r.p = NULL;
    new_bp->interval_l.p = NULL;

            new_bp->type = SV_BreakPoint::INVERSION;
            new_bp->type = SV_BreakPoint::DELETION;
            new_bp->type = SV_BreakPoint::DUPLICATION;
            new_bp->type = SV_BreakPoint::INVERSION;

    new_bp->weight = weight;

    return new_bp;
}

bool
SV_Pair::
is_aberrant()
{
    if ( read_l.strand == read_r.strand )
        return true;

    if ( read_l.strand == '-')
        return true;

    if ( (read_r.end - read_l.start) >=
        return true;

    if ( (read_r.end - read_l.start) <=
        return true;

    return false;
}

bool
SV_Pair::
is_sane()
{
    vector<UCSCElement<int> > v = exclude_regions.get(read_l.chr,
                                  read_l.start,
                                  read_l.end,
                                  '+',
                                  false);
    if ( v.size() > 0 )
        return false;

    v = exclude_regions.get(read_r.chr,
                            read_r.start,
                            read_r.end,
                            '+',
                            false);

    if ( v.size() > 0 )
        return false;


    if ( min_mapping_quality < reader->min_mapping_threshold )
        return false;

    int read_len_a = read_l.end - read_l.start;
    int read_len_b = read_r.end - read_r.start;


    cerr << "READ LEN\t" <<
            read_len_a << "\t" <<
            read_len_b << endl;

    if ( (read_len_a <= 0) || (read_len_b <= 0) )
        return false;


    int overlap = min(read_l.end, read_r.end) - max(read_l.start, read_r.start);
    cerr << "OVERLAP\t" <<
            read_l.chr << " " << 
            read_l.start  << " " << 
            read_l.end  << "\t" << 
            read_r.chr << " " << 
            read_r.start  << " " << 
            read_r.end  << "\t" <<
            overlap << endl;

    if (overlap > 0)
        return false;
    else
        return true;
    int non_overlap = min(read_len_a, read_len_b) - overlap;

    if ( (overlap > 0) && (abs(non_overlap) < reader->min_non_overlap) )
        return false;
    else
        return true;
}

ostream& operator << (ostream& out, const SV_Pair& p)
{

    out << p.read_l.chr << "," <<
        p.read_l.start << "," <<
        p.read_l.end << "," <<
        p.read_l.strand <<
        "\t" <<
        p.read_r.chr << "," <<
        p.read_r.start << "," <<
        p.read_r.end << "," <<
        p.read_r.strand;

    return out;
}

void
SV_Pair::
print_evidence()
{
    cout << read_id << "\t";
    print_bedpe(0);
}

void
SV_Pair::
print_bedpe(int score)
{
    string sep = "\t";
    cout <<
         read_l.chr << sep <<
         read_l.start << sep <<
         read_l.end << sep <<
         read_r.chr << sep <<
         read_r.start << sep<<
         read_r.end << sep<<
         this << sep <<
         score << sep <<
         read_l.strand << sep <<
         read_r.strand << sep <<
	 "id:" << ev_id << sep <<
         "weight:" << weight <<
         endl;
}

int
SV_Pair::
set_distro_from_histo (int back_distance,
                       int histo_start,
                       int histo_end,
{
    for (int i = 0; i < histo_end - histo_start + 1; ++i) {
        if (histo[i] == 0) {
            histo_end = i-1;
            break;
        }
    }

    int distro_size = back_distance + histo_end;


    for (int i = 0; i < back_distance; ++i)

    for (int i = back_distance; i < histo_start; ++i)

    double last = 0;
    for (int i = histo_end - 1; i >= histo_start; --i) {
    }


    int zero_i = -1;
    for (int i = 0; i < distro_size; ++i) {
            zero_i = i;
            break;
        }
    }

    if (zero_i != -1) {
        distro_size = zero_i;
        for (int i = 0; i < distro_size; ++i) {
        }
    }


    return distro_size;
}

void
SV_Pair::
process_pair(const BamAlignment &curr,
             const RefVector refs,
             map<string, BamAlignment> &mapped_pairs,
             int weight,
             int ev_id,
{
    if (mapped_pairs.find(curr.Name) == mapped_pairs.end())
        mapped_pairs[curr.Name] = curr;
    else {
                                        curr,
                                        refs,
                                        weight,
                                        ev_id,
                                        reader);
                
        if ( new_pair->is_sane() &&  
             new_pair->is_aberrant() &&
             (count_clipped(curr.CigarData) > 0) &&
             (count_clipped(mapped_pairs[curr.Name].CigarData) > 0) ) {
            cerr << "READ\t" << 
                    refs.at(mapped_pairs[curr.Name].RefID).RefName << "," <<
                    mapped_pairs[curr.Name].Position << "," <<
                    (mapped_pairs[curr.Name].GetEndPosition(false, false) - 1)
                        << "\t" <<
                    refs.at(curr.RefID).RefName << "," <<
                    curr.Position << "," <<
                    (curr.GetEndPosition(false, false) - 1)
                        <<
                    endl;

            new_bp->cluster(r_bin);
        } else {
            delete(new_pair);
        }
        mapped_pairs.erase(curr.Name);
    }
}

void
SV_Pair::
process_intra_chrom_pair(const BamAlignment &curr,
                         const RefVector refs,
                         BamWriter &inter_chrom_reads,
                         map<string, BamAlignment> &mapped_pairs,
                         int weight,
                         int ev_id,
{
    if (curr.RefID == curr.MateRefID) {

        process_pair(curr,
                     refs,
                     mapped_pairs,
                     r_bin,
                     weight,
                     ev_id,
                     reader);

    } else if (curr.IsMapped() &&
               curr.IsMateMapped() &&
               (curr.RefID >= 0) &&
               (curr.MateRefID >= 0) ) {
        BamAlignment al = curr;
        string x = reader->get_source_file_name();
        al.AddTag("LS","Z",x);
        inter_chrom_reads.SaveAlignment(al);
    }
}

string
SV_Pair::
evidence_type()
{
    return "Pair";
}

void
SV_Pair::
{
    for (unsigned int a = 0; a < m->size1 (); ++ a) {
        for (unsigned int b = 0; b < m->size2 (); ++ b) {
            unsigned int test_break_a, test_break_b, a_d, b_d;

            test_break_a = read_l.start + a;
            test_break_b = read_r.start + b;

            if (read_l.strand == '+')
                a_d = test_break_a - read_l.start;
            else
                a_d = read_l.end - test_break_a;


            if (read_r.strand == '+')
                b_d = test_break_b - read_r.start;
            else
                b_d = read_r.end - test_break_b;

            double dist = (a_d + b_d);
            double z_dist = ((a_d + b_d) - insert_mean);

            double p = gsl_ran_gaussian_pdf(z_dist, insert_stdev);
            log_space lp = get_ls(p);

        }
    }
}
void
SV_Pair::
set_distro_from_histo ()
{
    SV_Pair::distro_size = SV_Pair::back_distance + SV_Pair::histo_end;

    for (int i = 0; i < SV_Pair::back_distance; ++i)

    for (int i = SV_Pair::back_distance; i < SV_Pair::histo_start; ++i)
        SV_Pair::distro[i] = 1;

    double last = 0;
    for (int i = SV_Pair::histo_end - 1; i >= SV_Pair::histo_start; --i) {
        SV_Pair::distro[i + SV_Pair::back_distance] =
            SV_Pair::histo[i - SV_Pair::histo_start] + last;
        last = SV_Pair::distro[i + SV_Pair::back_distance];
    }
}


