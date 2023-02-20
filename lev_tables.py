import statistics

import Levenshtein
import pandas as pd
import numpy as np
import time
from math import sqrt
from misc import input_dialog_tables
import count_sam
import logging

log = logging.getLogger(__name__)


def make_table(sequence_json, alignment_json, output_dir):
    # Calculating the Levenshtein distance between sequence/alignment k-mers
    log.info('Calculating the Levenshtein distance between sequence/alignment k-mers')
    log.info('sequence_json: ' + sequence_json)
    log.info('alignment_json: ' + alignment_json)
    lev_index = pd.read_json(sequence_json, typ='series', orient='index')
    lev_columns = pd.read_json(alignment_json, typ='series', orient='records')
    # lev_table = pd.DataFrame(index=lev_index.keys(), columns=lev_columns.keys())
    lev_table = pd.DataFrame(index=lev_index.keys().str.upper(), columns=lev_columns.keys().str.upper())

    for i in range(lev_table.shape[0]):
        cutoff = None
        time_start = time.time()
        index = lev_table.iloc[i:i + 1, :].index[0]
        for c in range(lev_table.shape[1]):
            column = lev_table.iloc[:, c:c + 1].columns[0]
            extract = Levenshtein.distance(index[:cutoff], column[:cutoff])
            lev_table.iloc[i, c] = extract
            if lev_table.iloc[i, c] == 0:
                break
        time_end = time.time() - time_start
        log.debug('Row {} ({}-character cutoff) took {:.3f} s '.format(index[:cutoff], cutoff, time_end) +
                  'or {:.3f} m '.format(time_end / 60.0) +
                  'or {:.3f} h '.format(time_end / 3600.0))
    # lev_table = lev_table.apply(pd.to_numeric)
    lev_table = lev_table.astype(float)
    log.debug('\nlev_table\n' + str(lev_table))
    lev_table_fn = output_dir + '/tmp_aaf.' + 'levenshtein' + alignment_json.split('reads')[-1].split('json')[0] + \
                   sequence_json.split('.')[-3] + '.csv'
    lev_table.to_csv(lev_table_fn, index=True)
    log.debug('Levenshtein distances for all k-mers are written to ' + lev_table_fn)

    return lev_table


def calculating_lev_wagner_fischer(str1, str2):
    m = len(str1)
    n = len(str2)
    d = np.zeros(shape=(m + 1, n + 1), dtype=np.int)

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if str1[i - 1] == str2[j - 1]:
                sub = 0
            else:
                sub = 1
            d[i][j] = min(d[i - 1][j] + 1,  # deletion
                          d[i][j - 1] + 1,  # insertion
                          d[i - 1][j - 1] + sub  # substitution
                          )

    # log.debug('d[m][n] ' + d[m][n])
    return d[m][n]


def run_table(sequence_json=None, alignment_json=None, output_dir=None):
    # Preparing files
    log.info('Preparing files')
    if sequence_json is None or alignment_json is None or output_dir is None:
        sequence_json, alignment_json, output_dir = input_dialog_tables()

    time_start = time.time()
    lev_table = make_table(sequence_json, alignment_json, output_dir)
    time_end = time.time() - time_start
    log.debug('Calculation of the table took {:.3f} s'.format(time_end) +
              'or {:.3f} m'.format(time_end / 60.0) +
              'or {:.3f} h'.format(time_end / 3600.0))

    return output_dir, sequence_json, alignment_json, lev_table


def calc_metrics(ref_mean_start, ref_mean_end, ref_mean_mapq, alt_mean_start, alt_mean_end, alt_mean_mapq):
    euclidean_dist_nr = (ref_mean_start - alt_mean_start) ** 2 + (ref_mean_end - alt_mean_end) ** 2 + \
                        (ref_mean_mapq - alt_mean_mapq) ** 2
    euclidean_dist = sqrt((ref_mean_start - alt_mean_start) ** 2 + (ref_mean_end - alt_mean_end) ** 2 + \
                          (ref_mean_mapq - alt_mean_mapq) ** 2)
    manhattan_dist = abs(ref_mean_start - alt_mean_start) + abs(ref_mean_end - alt_mean_end) + \
                     abs(ref_mean_mapq - alt_mean_mapq)
    chebyshev_dist = max(abs(ref_mean_start - alt_mean_start), abs(ref_mean_end - alt_mean_end),
                         abs(ref_mean_mapq - alt_mean_mapq))
    return euclidean_dist_nr, euclidean_dist, manhattan_dist, chebyshev_dist


def calc_metrics_4(ref_mean_start, ref_mean_end, ref_mean_mapq, alt_mean_start, alt_mean_end, alt_mean_mapq,
                   ref_mode_e2e, alt_mode_e2e):
    euclidean_dist_nr = (ref_mean_start - alt_mean_start) ** 2 + (ref_mean_end - alt_mean_end) ** 2 + \
                        (ref_mean_mapq - alt_mean_mapq) ** 2 + (ref_mode_e2e - alt_mode_e2e) ** 2
    euclidean_dist = sqrt((ref_mean_start - alt_mean_start) ** 2 + (ref_mean_end - alt_mean_end) ** 2 + \
                          (ref_mean_mapq - alt_mean_mapq) ** 2 + (ref_mode_e2e - alt_mode_e2e) ** 2)
    manhattan_dist = abs(ref_mean_start - alt_mean_start) + abs(ref_mean_end - alt_mean_end) + \
                     abs(ref_mean_mapq - alt_mean_mapq) + abs(ref_mode_e2e - alt_mode_e2e)
    chebyshev_dist = max(abs(ref_mean_start - alt_mean_start), abs(ref_mean_end - alt_mean_end),
                         abs(ref_mean_mapq - alt_mean_mapq)), abs(ref_mode_e2e - alt_mode_e2e)
    return euclidean_dist_nr, euclidean_dist, manhattan_dist, chebyshev_dist


def append_with_metrics_for_all(df, variant,
                                euclidean_dist_nr, euclidean_dist, manhattan_dist, chebyshev_dist,
                                exp_euclidean_dist_nr, exp_euclidean_dist, exp_manhattan_dist, exp_chebyshev_dist,
                                index_ref, ref_seq_ext, closest_ref_kmer, closest_ref_kmer_count,
                                index_alt, alt_seq_ext, closest_alt_kmer, closest_alt_kmer_count,
                                ref_mean_start, ref_mean_end, ref_mean_mapq,
                                alt_mean_start, alt_mean_end, alt_mean_mapq,
                                ref_mode_e2e, alt_mode_e2e,
                                alternate_allele_frequency,
                                exp_index_ref, exp_ref_seq_ext, exp_closest_ref_kmer, exp_closest_ref_kmer_count,
                                exp_index_alt, exp_alt_seq_ext, exp_closest_alt_kmer, exp_closest_alt_kmer_count,
                                exp_ref_mean_start, exp_ref_mean_end, exp_ref_mean_mapq,
                                exp_alt_mean_start, exp_alt_mean_end, exp_alt_mean_mapq,
                                exp_ref_mode_e2e, exp_alt_mode_e2e,
                                exp_alternate_allele_frequency
                                ):
    return df.append(
        {'variant': variant,
         'index_ref': index_ref,
         'ref_seq_ext': ref_seq_ext,
         'closest_ref_kmer': closest_ref_kmer,
         'closest_ref_kmer_count': closest_ref_kmer_count,
         'ref_mean_start': ref_mean_start,
         'ref_mean_end': ref_mean_end,
         'ref_mean_mapq': ref_mean_mapq,
         'ref_mode_e2e': ref_mode_e2e,
         'index_alt': index_alt,
         'alt_seq_ext': alt_seq_ext,
         'closest_alt_kmer': closest_alt_kmer,
         'closest_alt_kmer_count': closest_alt_kmer_count,
         'alt_mean_start': alt_mean_start,
         'alt_mean_end': alt_mean_end,
         'alt_mean_mapq': alt_mean_mapq,
         'alt_mode_e2e': alt_mode_e2e,
         'exp_index_ref': exp_index_ref,
         'exp_ref_seq_ext': exp_ref_seq_ext,
         'exp_closest_ref_kmer': exp_closest_ref_kmer,
         'exp_closest_ref_kmer_count': exp_closest_ref_kmer_count,
         'exp_ref_mean_start': exp_ref_mean_start,
         'exp_ref_mean_end': exp_ref_mean_end,
         'exp_ref_mean_mapq': exp_ref_mean_mapq,
         'exp_ref_mode_e2e': exp_ref_mode_e2e,
         'exp_index_alt': exp_index_alt,
         'exp_alt_seq_ext': exp_alt_seq_ext,
         'exp_closest_alt_kmer': exp_closest_alt_kmer,
         'exp_closest_alt_kmer_count': exp_closest_alt_kmer_count,
         'exp_alt_mean_start': exp_alt_mean_start,
         'exp_alt_mean_end': exp_alt_mean_end,
         'exp_alt_mean_mapq': exp_alt_mean_mapq,
         'exp_alt_mode_e2e': exp_alt_mode_e2e,
         'AAF': alternate_allele_frequency,
         'exp_AAF': exp_alternate_allele_frequency,
         'euclidean_dist_nr': euclidean_dist_nr,
         'euclidean_dist': euclidean_dist,
         'manhattan_dist': manhattan_dist,
         'chebyshev_dist': chebyshev_dist,
         'exp_euclidean_dist_nr': exp_euclidean_dist_nr,
         'exp_euclidean_dist': exp_euclidean_dist,
         'exp_manhattan_dist': exp_manhattan_dist,
         'exp_chebyshev_dist': exp_chebyshev_dist
         })


def append_with_metrics_for_exp(df, variant,
                                index_ref, ref_seq_ext, closest_ref_kmer, closest_ref_kmer_count,
                                index_alt, alt_seq_ext, closest_alt_kmer, closest_alt_kmer_count,
                                exp_euclidean_dist_nr, exp_euclidean_dist, exp_manhattan_dist, exp_chebyshev_dist,
                                exp_index_ref, exp_ref_seq_ext, exp_closest_ref_kmer, exp_closest_ref_kmer_count,
                                exp_index_alt, exp_alt_seq_ext, exp_closest_alt_kmer, exp_closest_alt_kmer_count,
                                exp_ref_mean_start, exp_ref_mean_end, exp_ref_mean_mapq,
                                exp_alt_mean_start, exp_alt_mean_end, exp_alt_mean_mapq,
                                exp_ref_mode_e2e, exp_alt_mode_e2e,
                                exp_alternate_allele_frequency, alternate_allele_frequency):
    return df.append(
        {'variant': variant,
         'index_ref': index_ref,
         'ref_seq_ext': ref_seq_ext,
         'closest_ref_kmer': closest_ref_kmer,
         'closest_ref_kmer_count': closest_ref_kmer_count,
         'index_alt': index_alt,
         'alt_seq_ext': alt_seq_ext,
         'closest_alt_kmer': closest_alt_kmer,
         'closest_alt_kmer_count': closest_alt_kmer_count,
         'exp_index_ref': exp_index_ref,
         'exp_ref_seq_ext': exp_ref_seq_ext,
         'exp_closest_ref_kmer': exp_closest_ref_kmer,
         'exp_closest_ref_kmer_count': exp_closest_ref_kmer_count,
         'exp_ref_mean_start': exp_ref_mean_start,
         'exp_ref_mean_end': exp_ref_mean_end,
         'exp_ref_mean_mapq': exp_ref_mean_mapq,
         'exp_ref_mode_e2e': exp_ref_mode_e2e,
         'exp_index_alt': exp_index_alt,
         'exp_alt_seq_ext': exp_alt_seq_ext,
         'exp_closest_alt_kmer': exp_closest_alt_kmer,
         'exp_closest_alt_kmer_count': exp_closest_alt_kmer_count,
         'exp_alt_mean_start': exp_alt_mean_start,
         'exp_alt_mean_end': exp_alt_mean_end,
         'exp_alt_mean_mapq': exp_alt_mean_mapq,
         'exp_alt_mode_e2e': exp_alt_mode_e2e,
         'exp_AAF': exp_alternate_allele_frequency,
         'AAF': alternate_allele_frequency,
         'exp_euclidean_dist_nr': exp_euclidean_dist_nr,
         'exp_euclidean_dist': exp_euclidean_dist,
         'exp_manhattan_dist': exp_manhattan_dist,
         'exp_chebyshev_dist': exp_chebyshev_dist
         })


def append_with_metrics_for_nom(df, variant,
                                exp_index_ref, exp_ref_seq_ext, exp_closest_ref_kmer, exp_closest_ref_kmer_count,
                                exp_index_alt, exp_alt_seq_ext, exp_closest_alt_kmer, exp_closest_alt_kmer_count,
                                euclidean_dist_nr, euclidean_dist, manhattan_dist, chebyshev_dist,
                                index_ref, ref_seq_ext, closest_ref_kmer, closest_ref_kmer_count,
                                index_alt, alt_seq_ext, closest_alt_kmer, closest_alt_kmer_count,
                                ref_mean_start, ref_mean_end, ref_mean_mapq,
                                alt_mean_start, alt_mean_end, alt_mean_mapq,
                                ref_mode_e2e, alt_mode_e2e,
                                alternate_allele_frequency, exp_alternate_allele_frequency):
    return df.append(
        {'variant': variant,
         'index_ref': index_ref,
         'ref_seq_ext': ref_seq_ext,
         'closest_ref_kmer': closest_ref_kmer,
         'closest_ref_kmer_count': closest_ref_kmer_count,
         'ref_mean_start': ref_mean_start,
         'ref_mean_end': ref_mean_end,
         'ref_mean_mapq': ref_mean_mapq,
         'ref_mode_e2e': ref_mode_e2e,
         'index_alt': index_alt,
         'alt_seq_ext': alt_seq_ext,
         'closest_alt_kmer': closest_alt_kmer,
         'closest_alt_kmer_count': closest_alt_kmer_count,
         'alt_mean_start': alt_mean_start,
         'alt_mean_end': alt_mean_end,
         'alt_mean_mapq': alt_mean_mapq,
         'alt_mode_e2e': alt_mode_e2e,
         'exp_index_ref': exp_index_ref,
         'exp_ref_seq_ext': exp_ref_seq_ext,
         'exp_closest_ref_kmer': exp_closest_ref_kmer,
         'exp_closest_ref_kmer_count': exp_closest_ref_kmer_count,
         'exp_index_alt': exp_index_alt,
         'exp_alt_seq_ext': exp_alt_seq_ext,
         'exp_closest_alt_kmer': exp_closest_alt_kmer,
         'exp_closest_alt_kmer_count': exp_closest_alt_kmer_count,
         'AAF': alternate_allele_frequency,
         'exp_AAF': exp_alternate_allele_frequency,
         'euclidean_dist_nr': euclidean_dist_nr,
         'euclidean_dist': euclidean_dist,
         'manhattan_dist': manhattan_dist,
         'chebyshev_dist': chebyshev_dist
         })


def append_without_metrics(df, variant,
                           index_ref, ref_seq_ext, closest_ref_kmer, closest_ref_kmer_count,
                           index_alt, alt_seq_ext, closest_alt_kmer, closest_alt_kmer_count,
                           alternate_allele_frequency,
                           exp_index_ref, exp_ref_seq_ext, exp_closest_ref_kmer, exp_closest_ref_kmer_count,
                           exp_index_alt, exp_alt_seq_ext, exp_closest_alt_kmer, exp_closest_alt_kmer_count,
                           exp_alternate_allele_frequency):
    return df.append({'variant': variant,
                      'index_ref': index_ref,
                      'ref_seq_ext': ref_seq_ext,
                      'closest_ref_kmer': closest_ref_kmer,
                      'closest_ref_kmer_count': closest_ref_kmer_count,
                      'index_alt': index_alt,
                      'alt_seq_ext': alt_seq_ext,
                      'closest_alt_kmer': closest_alt_kmer,
                      'closest_alt_kmer_count': closest_alt_kmer_count,
                      'exp_index_ref': exp_index_ref,
                      'exp_ref_seq_ext': exp_ref_seq_ext,
                      'exp_closest_ref_kmer': exp_closest_ref_kmer,
                      'exp_closest_ref_kmer_count': exp_closest_ref_kmer_count,
                      'exp_index_alt': exp_index_alt,
                      'exp_alt_seq_ext': exp_alt_seq_ext,
                      'exp_closest_alt_kmer': exp_closest_alt_kmer,
                      'exp_closest_alt_kmer_count': exp_closest_alt_kmer_count,
                      'AAF': alternate_allele_frequency,
                      'exp_AAF': exp_alternate_allele_frequency
                      })


def variant_in_kmer(variant, kmer):
    # print(variant, ':', kmer, ':', type(variant), type(kmer))
    # return [str(i) for i in range(len(kmer)) if variant in kmer[i:i + len(variant)]]
    returned = [str(i) for i in range(len(str(kmer))) if str(variant) in str(kmer)[i:i + len(str(variant))]]
    # print('returned', returned, '\n')
    return returned


def intersect_lists(list1, list2):
    intersected = list(set(list1).intersection(list2))
    return intersected if intersected else None


def find_aaf(output_dir, alignment_json, lev_table_ref, lev_table_alt, mini_sam_fn=None, variant=None, k=None):
    # Selecting ref/alt sequence k-mers with min distance to their alignment counterparts
    # AAF = counts of alt k-mer in the alignment divided by the sum of alt+ref k-mer counts in the alignment
    fa_1_ts = time.time()
    log.info('Selecting ref/alt sequence k-mers with min distance to their alignment counterparts')

    aaf = []
    # aj_fn = alignment_json[:-5].split('alignment')
    alignment_kmer_counts = pd.read_json(alignment_json, typ='series', orient='records')
    alignment_kmer_counts.keys = alignment_kmer_counts.keys().str.upper()
    dp = count_sam.counting_dp(mini_sam_fn, variant)

    lev_exp_table_ref = lev_table_ref.copy()
    lev_exp_table_alt = lev_table_alt.copy()

    fa_1_te = time.time() - fa_1_ts
    fa_2_ts = time.time()
    for column in lev_exp_table_ref.columns:
        lev_exp_table_ref[column] = (2 ** (k - lev_exp_table_ref[column])) * (alignment_kmer_counts[column] / dp)
        lev_exp_table_alt[column] = (2 ** (k - lev_exp_table_alt[column])) * (alignment_kmer_counts[column] / dp)
    fa_2_te = time.time() - fa_2_ts

    # print('k', k, 'dp', dp)
    # print('\nlev_exp_table_alt\n', lev_exp_table_alt)
    # print('\nlev_exp_table_alt\n', lev_exp_table_alt)

    fa_3_ts = time.time()
    fa_4_mean = []
    fa_5_mean = []
    for index in range(len(lev_table_ref.index)):
        fa_4_loop_ts = time.time()
        index_ref = lev_table_ref.index[index]
        index_alt = lev_table_alt.index[index]
        exp_index_ref = lev_exp_table_ref.index[index]
        exp_index_alt = lev_exp_table_alt.index[index]
        log.debug('index_ref ' + index_ref + ' index_alt ' + index_alt)

        ref_seq_ext = lev_table_ref.iloc[index, :].min()  # loc[index_ref, :]
        alt_seq_ext = lev_table_alt.iloc[index, :].min()  # loc[index_alt, :]
        closest_ref_kmer = lev_table_ref.iloc[index, :].idxmin()
        closest_alt_kmer = lev_table_alt.iloc[index, :].idxmin()
        exp_ref_seq_ext = lev_exp_table_ref.iloc[index, :].max()  # loc[exp_index_ref, :]
        exp_alt_seq_ext = lev_exp_table_alt.iloc[index, :].max()  # loc[exp_index_alt, :]
        exp_closest_ref_kmer = lev_exp_table_ref.iloc[index, :].idxmax()
        exp_closest_alt_kmer = lev_exp_table_alt.iloc[index, :].idxmax()

        # print('Levenshtein.distance between closest_ref_kmer and closest_alt_kmer',
        #       Levenshtein.distance(closest_ref_kmer, closest_alt_kmer))

        # index_ref_kmer_count = alignment_kmer_counts[index_ref]
        # index_alt_kmer_count = alignment_kmer_counts[index_alt]
        closest_ref_kmer_count = alignment_kmer_counts[closest_ref_kmer]
        closest_alt_kmer_count = alignment_kmer_counts[closest_alt_kmer]
        exp_closest_ref_kmer_count = alignment_kmer_counts[exp_closest_ref_kmer]
        exp_closest_alt_kmer_count = alignment_kmer_counts[exp_closest_alt_kmer]
        # print('closest_ref_kmer_count', closest_ref_kmer_count, 'closest_alt_kmer_count', closest_alt_kmer_count)
        # log.debug('ref_seq_ext ' + str(ref_seq_ext) + ' closest_ref_kmer ' + closest_ref_kmer + ' index_ref==closest_ref_kmer ' + str(index_ref==closest_ref_kmer))
        # log.debug('alt_seq_ext ' + str(alt_seq_ext) + ' closest_alt_kmer ' + closest_alt_kmer + ' index_alt==closest_alt_kmer ' + str(index_alt==closest_alt_kmer))
        # Additional conditions such as below could be applied:
        # ref_seq_ext > alt_seq_ext, closest_ref_kmer_count == closest_alt_kmer_count, index_ref==closest_ref_kmer
        # log.info('alignment_kmer_counts.index == alignment_kmer_counts.keys' if alignment_kmer_counts.index.to_list == alignment_kmer_counts.keys.to_list else 'alignment_kmer_counts.index != alignment_kmer_counts.keys')
        # log.info('alignment_kmer_counts.index.to_list ' + str(alignment_kmer_counts.index.to_list))
        # log.info('alignment_kmer_counts.keys.to_list ' + str(alignment_kmer_counts.keys.to_list))

        index_log_1 = 'index_ref == index_alt' if index_ref == index_alt else 'index_ref != index_alt'
        index_log_2 = ' | index_ref not in alignment_kmer_counts.index' if index_ref not in alignment_kmer_counts.index else ' | index_ref in alignment_kmer_counts.index'
        index_log_3 = ' | index_alt not in alignment_kmer_counts.index' if index_alt not in alignment_kmer_counts.index else ' | index_alt in alignment_kmer_counts.index'
        closest_log_1 = 'closest_ref_kmer == closest_alt_kmer' if closest_ref_kmer == closest_alt_kmer else 'closest_ref_kmer != closest_alt_kmer'
        closest_log_2 = ' | closest_ref_kmer not in alignment_kmer_counts.index' if closest_ref_kmer not in alignment_kmer_counts.index else ' | closest_ref_kmer in alignment_kmer_counts.index'
        closest_log_3 = ' | closest_alt_kmer not in alignment_kmer_counts.index' if closest_alt_kmer not in alignment_kmer_counts.index else ' | closest_alt_kmer in alignment_kmer_counts.index'
        exp_closest_log_1 = 'exp_closest_ref_kmer == exp_closest_alt_kmer' if exp_closest_ref_kmer == exp_closest_alt_kmer else 'exp_closest_ref_kmer != exp_closest_alt_kmer'
        exp_closest_log_2 = ' | exp_closest_ref_kmer not in alignment_kmer_counts.index' if exp_closest_ref_kmer not in alignment_kmer_counts.index else ' | exp_closest_ref_kmer in alignment_kmer_counts.index'
        exp_closest_log_3 = ' | exp_closest_alt_kmer not in alignment_kmer_counts.index' if exp_closest_alt_kmer not in alignment_kmer_counts.index else ' | exp_closest_alt_kmer in alignment_kmer_counts.index'

        log.debug(index_log_1 + index_log_2 + index_log_3)
        log.debug(closest_log_1 + closest_log_2 + closest_log_3)
        log.debug(exp_closest_log_1 + exp_closest_log_2 + exp_closest_log_3)
        log.info(index_log_1 + ' || ' + closest_log_1)

        if closest_alt_kmer == closest_ref_kmer:  # index_ref == index_alt:
            alternate_allele_frequency = closest_alt_kmer_count / closest_alt_kmer_count
            log.info('[Nominal] AAF = closest_alt_kmer_count / closest_alt_kmer_count'
                     ' = ' + str(closest_alt_kmer_count) + ' / ' + str(closest_alt_kmer_count))
            log.info('[Nominal] Apparently, closest_alt_kmer == closest_ref_kmer: may have to reduce the k-mer length')
            # return False, True, [None, None], [None, None]
        else:
            alternate_allele_frequency = closest_alt_kmer_count / (closest_alt_kmer_count + closest_ref_kmer_count)
            log.info('[Nominal] AAF = closest_alt_kmer_count / (closest_alt_kmer_count + closest_ref_kmer_count)'
                     ' = ' + str(closest_alt_kmer_count) + ' / (' +
                     str(closest_alt_kmer_count) + ' + ' + str(closest_ref_kmer_count) + ')')
        log.info('[Nominal] AAF for ' + closest_ref_kmer + ' (ref) and ' + closest_alt_kmer + ' (alt) calculated: '
                 + str(alternate_allele_frequency))

        if exp_closest_alt_kmer == exp_closest_ref_kmer:  # index_ref == index_alt:
            exp_alternate_allele_frequency = exp_closest_alt_kmer_count / exp_closest_alt_kmer_count
            log.info('[Expanded] AAF = exp_closest_alt_kmer_count / exp_closest_alt_kmer_count'
                     ' = ' + str(exp_closest_alt_kmer_count) + ' / ' + str(exp_closest_alt_kmer_count))
            log.info(
                '[Expanded] Apparently, exp_closest_alt_kmer == exp_closest_ref_kmer: may have to reduce the k-mer length')
            # return False, True, [None, None], [None, None]
        else:
            exp_alternate_allele_frequency = exp_closest_alt_kmer_count / (
                    exp_closest_alt_kmer_count + exp_closest_ref_kmer_count)
            log.info(
                '[Expanded] AAF = exp_closest_alt_kmer_count / (exp_closest_alt_kmer_count + exp_closest_ref_kmer_count)'
                ' = ' + str(exp_closest_alt_kmer_count) + ' / (' +
                str(exp_closest_alt_kmer_count) + ' + ' + str(exp_closest_ref_kmer_count) + ')')
        log.info(
            '[Expanded] AAF for ' + exp_closest_ref_kmer + ' (ref) and ' + exp_closest_alt_kmer + ' (alt) calculated: '
            + str(exp_alternate_allele_frequency))
        fa_4_loop_te = time.time() - fa_4_loop_ts

        fa_5_loop_ts = time.time()
        # Using Euclidean distance (nonroot) of read properties to determine fitness of the selected k-mers
        if mini_sam_fn and variant:  # and alternate_allele_frequency < 1.0:
            variant_a = ''.join(char for char in variant.split(':')[-1] if not char.isnumeric()).split('>')
            chunk = 2
            com_ref_kmer = closest_ref_kmer.upper()
            com_alt_kmer = closest_ref_kmer.upper()
            com_ref_var = [variant_a[0].upper()[_:_ + chunk] for _ in range(0, len(variant_a[0]), chunk)]
            com_alt_var = [variant_a[-1].upper()[_:_ + chunk] for _ in range(0, len(variant_a[-1]), chunk)]
            # if any(sub in com_ref_kmer for sub in [v for v in com_ref_var]) or any(sub in com_alt_kmer for sub in [v for v in com_alt_var]):
            # if variant_a[0].upper() in closest_ref_kmer.upper() and variant_a[-1].upper() in closest_alt_kmer.upper():

            ref_mean_start, ref_mean_end, ref_mean_mapq, ref_mode_e2e = count_sam.read_properties(mini_sam_fn,
                                                                                                  closest_ref_kmer)
            alt_mean_start, alt_mean_end, alt_mean_mapq, alt_mode_e2e = count_sam.read_properties(mini_sam_fn,
                                                                                                  closest_alt_kmer)
            exp_ref_mean_start, exp_ref_mean_end, exp_ref_mean_mapq, exp_ref_mode_e2e = count_sam.read_properties(mini_sam_fn,
                                                                                                                  exp_closest_ref_kmer)
            exp_alt_mean_start, exp_alt_mean_end, exp_alt_mean_mapq, exp_alt_mode_e2e = count_sam.read_properties(mini_sam_fn,
                                                                                                                  exp_closest_alt_kmer)

            nom_condition_1 = ref_mean_start and ref_mean_end and ref_mean_mapq and alt_mean_start and alt_mean_end and alt_mean_mapq
            nom_condition_2 = ref_mean_start != alt_mean_start and ref_mean_end != alt_mean_end and ref_mean_mapq != alt_mean_mapq
            exp_condition_1 = exp_ref_mean_start and exp_ref_mean_end and exp_ref_mean_mapq and exp_alt_mean_start and exp_alt_mean_end and exp_alt_mean_mapq
            exp_condition_2 = exp_ref_mean_start != exp_alt_mean_start and exp_ref_mean_end != exp_alt_mean_end and exp_ref_mean_mapq != exp_alt_mean_mapq
            if nom_condition_1 and nom_condition_2 and exp_condition_1 and exp_condition_2:
                euclidean_dist_nr, euclidean_dist, manhattan_dist, chebyshev_dist = calc_metrics(ref_mean_start,
                                                                                                 ref_mean_end,
                                                                                                 ref_mean_mapq,
                                                                                                 alt_mean_start,
                                                                                                 alt_mean_end,
                                                                                                 alt_mean_mapq)
                # euclidean_dist_nr, euclidean_dist, manhattan_dist, chebyshev_dist = calc_metrics_4(ref_mean_start,
                #                                                                                  ref_mean_end,
                #                                                                                  ref_mean_mapq,
                #                                                                                  alt_mean_start,
                #                                                                                  alt_mean_end,
                #                                                                                  alt_mean_mapq,
                #                                                                                  ref_mode_e2e,
                #                                                                                  alt_mode_e2e)
                exp_euclidean_dist_nr, exp_euclidean_dist, exp_manhattan_dist, exp_chebyshev_dist = calc_metrics(
                    exp_ref_mean_start,
                    exp_ref_mean_end,
                    exp_ref_mean_mapq,
                    exp_alt_mean_start,
                    exp_alt_mean_end,
                    exp_alt_mean_mapq)
                # exp_euclidean_dist_nr, exp_euclidean_dist, exp_manhattan_dist, exp_chebyshev_dist = calc_metrics_4(
                #     exp_ref_mean_start,
                #     exp_ref_mean_end,
                #     exp_ref_mean_mapq,
                #     exp_alt_mean_start,
                #     exp_alt_mean_end,
                #     exp_alt_mean_mapq,
                #     exp_ref_mode_e2e, exp_alt_mode_e2e)
                append_with_metrics_for_all(aaf, variant,
                                            euclidean_dist_nr, euclidean_dist, manhattan_dist, chebyshev_dist,
                                            exp_euclidean_dist_nr, exp_euclidean_dist, exp_manhattan_dist,
                                            exp_chebyshev_dist,
                                            index_ref, ref_seq_ext, closest_ref_kmer, closest_ref_kmer_count,
                                            index_alt, alt_seq_ext, closest_alt_kmer, closest_alt_kmer_count,
                                            ref_mean_start, ref_mean_end, ref_mean_mapq,
                                            alt_mean_start, alt_mean_end, alt_mean_mapq,
                                            ref_mode_e2e, alt_mode_e2e,
                                            alternate_allele_frequency,
                                            exp_index_ref, exp_ref_seq_ext, exp_closest_ref_kmer,
                                            exp_closest_ref_kmer_count,
                                            exp_index_alt, exp_alt_seq_ext, exp_closest_alt_kmer,
                                            exp_closest_alt_kmer_count,
                                            exp_ref_mean_start, exp_ref_mean_end, exp_ref_mean_mapq,
                                            exp_alt_mean_start, exp_alt_mean_end, exp_alt_mean_mapq,
                                            exp_ref_mode_e2e, exp_alt_mode_e2e,
                                            exp_alternate_allele_frequency)
            elif nom_condition_1 and nom_condition_2 and (
                    not exp_condition_1 or not exp_condition_2):  # and was here, no brackets
                # log.debug('ref_mean_start ' + str(ref_mean_start) + ' ref_mean_end ' + str(ref_mean_end) + ' ref_mean_mapq ' + str(ref_mean_mapq))
                # log.debug('alt_mean_start ' + str(alt_mean_start) + ' alt_mean_end ' + str(alt_mean_end) + ' alt_mean_mapq ' + str(alt_mean_mapq))
                euclidean_dist_nr, euclidean_dist, manhattan_dist, chebyshev_dist = calc_metrics(ref_mean_start,
                                                                                                 ref_mean_end,
                                                                                                 ref_mean_mapq,
                                                                                                 alt_mean_start,
                                                                                                 alt_mean_end,
                                                                                                 alt_mean_mapq)
                # euclidean_dist_nr, euclidean_dist, manhattan_dist, chebyshev_dist = calc_metrics_4(ref_mean_start,
                #                                                                                  ref_mean_end,
                #                                                                                  ref_mean_mapq,
                #                                                                                  alt_mean_start,
                #                                                                                  alt_mean_end,
                #                                                                                  alt_mean_mapq,
                #                                                                                    ref_mode_e2e,
                #                                                                                    alt_mode_e2e)
                append_with_metrics_for_nom(aaf, variant,
                                            exp_index_ref, exp_ref_seq_ext, exp_closest_ref_kmer, exp_closest_ref_kmer_count,
                                            exp_index_alt, exp_alt_seq_ext, exp_closest_alt_kmer, exp_closest_alt_kmer_count,
                                            euclidean_dist_nr, euclidean_dist, manhattan_dist, chebyshev_dist,
                                            index_ref, ref_seq_ext, closest_ref_kmer, closest_ref_kmer_count,
                                            index_alt, alt_seq_ext, closest_alt_kmer, closest_alt_kmer_count,
                                            ref_mean_start, ref_mean_end, ref_mean_mapq,
                                            alt_mean_start, alt_mean_end, alt_mean_mapq,
                                            ref_mode_e2e, alt_mode_e2e,
                                            alternate_allele_frequency, exp_alternate_allele_frequency)

                log_start = 'The metric was not calculated because'
                log_same = ' the read properties for reference and alternative k-mers are the same'
                log_exp_ref = ' the reference sequence k-mer doesn\'t have supporting reads: ' + exp_closest_ref_kmer
                log_exp_alt = ' the alternative sequence k-mer doesn\'t have supporting reads: ' + exp_closest_alt_kmer
                if not exp_ref_mean_start or not exp_ref_mean_end or not exp_ref_mean_mapq:
                    log.info('[Expanded] ' + log_start + log_exp_ref)
                if not exp_alt_mean_start or not exp_alt_mean_end or not exp_alt_mean_mapq:
                    log.info('[Expanded] ' + log_start + log_exp_alt)
                if exp_ref_mean_start == exp_alt_mean_start and exp_ref_mean_end == exp_alt_mean_end and \
                        exp_ref_mean_mapq == exp_alt_mean_mapq:
                    log.info('[Expanded] ' + log_start + log_same)

            elif (
                    not nom_condition_1 or not nom_condition_2) and exp_condition_1 and exp_condition_2:  # and was here, no brackets
                exp_euclidean_dist_nr, exp_euclidean_dist, exp_manhattan_dist, exp_chebyshev_dist = calc_metrics(
                    exp_ref_mean_start,
                    exp_ref_mean_end,
                    exp_ref_mean_mapq,
                    exp_alt_mean_start,
                    exp_alt_mean_end,
                    exp_alt_mean_mapq)
                # exp_euclidean_dist_nr, exp_euclidean_dist, exp_manhattan_dist, exp_chebyshev_dist = calc_metrics_4(
                #     exp_ref_mean_start,
                #     exp_ref_mean_end,
                #     exp_ref_mean_mapq,
                #     exp_alt_mean_start,
                #     exp_alt_mean_end,
                #     exp_alt_mean_mapq,
                #     exp_ref_mode_e2e, exp_alt_mode_e2e)
                append_with_metrics_for_exp(aaf, variant,
                                            index_ref, ref_seq_ext, closest_ref_kmer, closest_ref_kmer_count,
                                            index_alt, alt_seq_ext, closest_alt_kmer, closest_alt_kmer_count,
                                            exp_euclidean_dist_nr, exp_euclidean_dist,
                                            exp_manhattan_dist, exp_chebyshev_dist,
                                            exp_index_ref, exp_ref_seq_ext, exp_closest_ref_kmer,
                                            exp_closest_ref_kmer_count,
                                            exp_index_alt, exp_alt_seq_ext, exp_closest_alt_kmer,
                                            exp_closest_alt_kmer_count,
                                            exp_ref_mode_e2e, exp_alt_mode_e2e,
                                            exp_ref_mean_start, exp_ref_mean_end, exp_ref_mean_mapq,
                                            exp_alt_mean_start, exp_alt_mean_end, exp_alt_mean_mapq,
                                            exp_alternate_allele_frequency, alternate_allele_frequency)

                log_start = 'The metric was not calculated because'
                log_same = ' the read properties for reference and alternative k-mers are the same'
                log_ref = ' the reference sequence k-mer doesn\'t have supporting reads: ' + closest_ref_kmer
                log_alt = ' the alternative sequence k-mer doesn\'t have supporting reads: ' + closest_alt_kmer
                if not ref_mean_start or not ref_mean_end or not ref_mean_mapq:
                    log.info('[Nominal] ' + log_start + log_ref)
                if not alt_mean_start or not alt_mean_end or not alt_mean_mapq:
                    log.info('[Nominal] ' + log_start + log_alt)
                if ref_mean_start == alt_mean_start and ref_mean_end == alt_mean_end and ref_mean_mapq == alt_mean_mapq:
                    log.info('[Nominal] ' + log_start + log_same)

            else:
                log_start = 'The metric was not calculated because'
                log_same = ' the read properties for reference and alternative k-mers are the same'

                log_ref = ' the reference sequence k-mer doesn\'t have supporting reads: ' + closest_ref_kmer
                log_alt = ' the alternative sequence k-mer doesn\'t have supporting reads: ' + closest_alt_kmer
                if not ref_mean_start or not ref_mean_end or not ref_mean_mapq:
                    log.info('[Nominal] ' + log_start + log_ref)
                if not alt_mean_start or not alt_mean_end or not alt_mean_mapq:
                    log.info('[Nominal] ' + log_start + log_alt)
                if ref_mean_start == alt_mean_start and ref_mean_end == alt_mean_end and ref_mean_mapq == alt_mean_mapq:
                    log.info('[Nominal] ' + log_start + log_same)

                log_exp_ref = ' the reference sequence k-mer doesn\'t have supporting reads: ' + exp_closest_ref_kmer
                log_exp_alt = ' the alternative sequence k-mer doesn\'t have supporting reads: ' + exp_closest_alt_kmer
                if not exp_ref_mean_start or not exp_ref_mean_end or not exp_ref_mean_mapq:
                    log.info('[Expanded] ' + log_start + log_exp_ref)
                if not exp_alt_mean_start or not exp_alt_mean_end or not exp_alt_mean_mapq:
                    log.info('[Expanded] ' + log_start + log_exp_alt)
                if exp_ref_mean_start == exp_alt_mean_start and exp_ref_mean_end == exp_alt_mean_end and \
                        exp_ref_mean_mapq == exp_alt_mean_mapq:
                    log.info('[Expanded] ' + log_start + log_same)

                append_without_metrics(aaf, variant,
                                       index_ref, ref_seq_ext, closest_ref_kmer, closest_ref_kmer_count,
                                       index_alt, alt_seq_ext, closest_alt_kmer, closest_alt_kmer_count,
                                       alternate_allele_frequency,
                                       exp_index_ref, exp_ref_seq_ext, exp_closest_ref_kmer, exp_closest_ref_kmer_count,
                                       exp_index_alt, exp_alt_seq_ext, exp_closest_alt_kmer, exp_closest_alt_kmer_count,
                                       exp_alternate_allele_frequency)
            # else:
            #     log_start = 'Discarding the pair of k-mers because '
            #     log_ref = 'the closest ref k-mer cannot accommodate the ref allele: ' + variant_a[0].upper() + ' not in ' + closest_ref_kmer.upper()
            #     log_alt = 'the closest alt k-mer cannot accommodate the alt allele: ' + variant_a[-1].upper() + ' not in ' + closest_alt_kmer.upper()
            #     if variant_a[0].upper() not in closest_ref_kmer.upper():
            #         log.info(log_start + log_ref)
            #     if variant_a[-1].upper() not in closest_alt_kmer.upper():
            #         log.info(log_start + log_alt)
            #     # return False, False, 0
        else:
            if not mini_sam_fn:
                log.info('Since the alignment file wasn\'t provided, the metric will not be calculated')
            if not variant:
                log.info('Since the variant wasn\'t provided, the metric will not be calculated')

            append_without_metrics(aaf, variant,
                                   index_ref, ref_seq_ext, closest_ref_kmer, closest_ref_kmer_count,
                                   index_alt, alt_seq_ext, closest_alt_kmer, closest_alt_kmer_count,
                                   alternate_allele_frequency,
                                   exp_index_ref, exp_ref_seq_ext, exp_closest_ref_kmer, exp_closest_ref_kmer_count,
                                   exp_index_alt, exp_alt_seq_ext, exp_closest_alt_kmer, exp_closest_alt_kmer_count,
                                   exp_alternate_allele_frequency)
        fa_5_loop_te = time.time() - fa_5_loop_ts
        fa_4_mean.append(fa_4_loop_te)
        fa_5_mean.append(fa_5_loop_te)
    fa_3_te = time.time() - fa_3_ts
    fa_6_ts = time.time()
    aaf = pd.DataFrame(aaf)
    log.info('\nAAF\n' + str(aaf))
    if aaf.empty:
        log.info('Alternate Allele Frequency (AAF) values do not exist, which isn\'t acceptable')
        # log.debug('\nlev_table_ref\n' + str(lev_table_ref))
        # log.debug('\nlev_table_alt\n' + str(lev_table_alt))
        return False, True, [None, None], [None, None], [None, None], [None, None], \
                       [fa_1_te, fa_2_te, fa_3_te,
                        statistics.mean(fa_4_mean), statistics.mean(fa_5_mean),
                        0, 0, 0]
    else:
        # print('k', k)
        if 'closest_alt_kmer' in aaf.columns and all(aaf['closest_alt_kmer'] == aaf['closest_ref_kmer']) and k > 15:
            log.info('[Nominal] Apparently, all of closest_alt_kmer == closest_ref_kmer: reducing the k-mer length')
            return False, True, [None, None], [None, None], [None, None], [None, None], \
                       [fa_1_te, fa_2_te, fa_3_te,
                        statistics.mean(fa_4_mean), statistics.mean(fa_5_mean),
                        0, 0, 0]
        if 'exp_closest_alt_kmer' in aaf.columns and all(aaf['exp_closest_alt_kmer'] == aaf['exp_closest_ref_kmer']) and k > 15:
            log.info(
                '[Expanded] Apparently, all of exp_closest_alt_kmer == exp_closest_ref_kmer: reducing the k-mer length')
            return False, True, [None, None], [None, None], [None, None], [None, None], \
                       [fa_1_te, fa_2_te, fa_3_te,
                        statistics.mean(fa_4_mean), statistics.mean(fa_5_mean),
                        0, 0, 0]

        fa_6_te = time.time() - fa_6_ts
        fa_7_ts = time.time()
        # log.debug('\nlev_table_ref\n' + str(lev_table_ref))
        # log.debug('\nlev_table_alt\n' + str(lev_table_alt))
        aj_fn = alignment_json[:-5].split('alignment')
        aaf_fn = output_dir + '/' + 'alternate_allele_frequency' + '.k' + str(k) + aj_fn[-1] + '.csv'
        aaf.to_csv(aaf_fn, index=False)
        log.info('Alternate Allele Frequency (AAF) for pairs of k-mers is written to ' + aaf_fn)

        final_aaf_fn = output_dir + '/' + 'final_alternate_allele_frequency' + '.k' + str(k) + aj_fn[-1] + '.csv'

        # Selecting candidates with minimum Levenshtein distances
        # Seeking minimum first, then filtering by kmer counts
        aaf_nom = aaf[(aaf.loc[:, 'ref_seq_ext'] == aaf.loc[:, 'ref_seq_ext'].min()) &
                      (aaf.loc[:, 'alt_seq_ext'] == aaf.loc[:, 'alt_seq_ext'].min())]
        if aaf_nom.empty:
            aaf_nom = aaf[(aaf.loc[:, 'ref_seq_ext'] == aaf.loc[:, 'ref_seq_ext'].min()) |
                          (aaf.loc[:, 'alt_seq_ext'] == aaf.loc[:, 'alt_seq_ext'].min())]
        aaf_hc = aaf_nom[(aaf_nom.loc[:, 'closest_ref_kmer_count'] >= 20) &
                         (aaf_nom.loc[:, 'closest_alt_kmer_count'] >= 20)]
        if not aaf_hc.empty:
            aaf_nom = aaf_hc

        # Marking kmer pairs that include the variant
        # if variant_a:
        #     aaf_nom.loc[:, ['variant_pos_ref']] = aaf_nom.apply(lambda x: variant_in_kmer(variant_a[0], x['closest_ref_kmer']), axis=1)
        #     aaf_nom.loc[:, ['variant_pos_alt']] = aaf_nom.apply(lambda x: variant_in_kmer(variant_a[-1], x['closest_alt_kmer']), axis=1)
        #     aaf_nom.loc[:, ['variant_pos_intersection']] = aaf_nom.apply(lambda x: intersect_lists(x['variant_pos_ref'], x['variant_pos_alt']), axis=1)
        #     # print('aaf_nom', aaf_nom[['variant', 'closest_ref_kmer', 'closest_alt_kmer', 'variant_pos_ref', 'variant_pos_alt', 'variant_pos_intersection']])
        #     # print('aaf_nom.shape before dropna', aaf_nom.shape)
        #     print('aaf_nom none count', aaf_nom[aaf_nom['variant_pos_intersection'] == None]['variant_pos_intersection'].count())  # [aaf_nom.loc[:, 'variant_pos_intersection'] is None, ])
        #     if aaf_nom[aaf_nom['variant_pos_intersection'] == None]['variant_pos_intersection'].count() != aaf_nom.shape[0]:
        #         aaf_nom = aaf_nom.dropna(subset=['variant_pos_intersection'])
        #         print('removed nones in aaf_nom')
        #     # print('aaf_nom.shape before dropna', aaf_nom.shape)
        #
        #     aaf.loc[:, ['exp_variant_pos_ref']] = aaf.apply(lambda x: variant_in_kmer(variant_a[0], x['exp_closest_ref_kmer']), axis=1)
        #     aaf.loc[:, ['exp_variant_pos_alt']] = aaf.apply(lambda x: variant_in_kmer(variant_a[-1], x['exp_closest_alt_kmer']), axis=1)
        #     aaf.loc[:, ['exp_variant_pos_intersection']] = aaf.apply(lambda x: intersect_lists(x['exp_variant_pos_ref'], x['exp_variant_pos_alt']), axis=1)
        #     # print('aaf', aaf[['variant', 'exp_closest_ref_kmer', 'exp_closest_alt_kmer', 'exp_variant_pos_ref', 'exp_variant_pos_alt', 'exp_variant_pos_intersection']])
        #     # print('aaf.shape before dropna', aaf.shape)
        #     print('aaf none count', aaf[aaf['exp_variant_pos_intersection'] == None]['exp_variant_pos_intersection'].count())  # [aaf_nom.loc[:, 'exp_variant_pos_intersection'] is None, ].count)
        #     if aaf[aaf['exp_variant_pos_intersection'] == None]['exp_variant_pos_intersection'].count() != aaf.shape[0]:
        #         aaf = aaf.dropna(subset=['exp_variant_pos_intersection'])
        #         print('removed nones in aaf_nom')
        #     # print('aaf.shape after dropna', aaf.shape)
        # print('\naaf_nom after prioritizing kmer pairs that include the variant\n', aaf_nom)

        # Prioritizing kmers located closer to the end of the read (nom: 1 for snps, 3 for indels, exp: 3 for all)
        if 'ref_mode_e2e' in aaf_nom.columns:

            # print("aaf_nom['ref_mode_e2e'].median()", aaf_nom['ref_mode_e2e'].median())
            # print('ref_mode_e2e', aaf_nom['ref_mode_e2e'],
            #       'alt_mode_e2e', aaf_nom['alt_mode_e2e'])
            if variant_a and any((len(variant_a[0]) > 1, len(variant_a[-1]) > 1)):
                aaf_nom = aaf_nom.nsmallest(3, 'ref_mode_e2e')
                # print("aaf_nom['ref_mode_e2e'] == aaf_nom['alt_mode_e2e']", aaf_nom['ref_mode_e2e'] == aaf_nom['alt_mode_e2e'])
            else:
                aaf_nom = aaf_nom.nsmallest(1, 'ref_mode_e2e')
                # print("aaf_nom['ref_mode_e2e'] == aaf_nom['alt_mode_e2e']", aaf_nom['ref_mode_e2e'] == aaf_nom['alt_mode_e2e'])
        if 'exp_ref_mode_e2e' in aaf.columns:

            # print("aaf['exp_ref_mode_e2e'].median()", aaf['exp_ref_mode_e2e'].median())
            # print('exp_ref_mode_e2e', aaf['exp_ref_mode_e2e'],
            #       'exp_alt_mode_e2e', aaf['exp_alt_mode_e2e'])
            aaf = aaf.nsmallest(3, 'exp_ref_mode_e2e')
            # print("aaf['exp_ref_mode_e2e'] == aaf['exp_alt_mode_e2e']", aaf['exp_ref_mode_e2e'] == aaf['exp_alt_mode_e2e'])

        fa_7_te = time.time() - fa_7_ts
        fa_8_ts = time.time()
        # aaf_nom = aaf
        aaf_nom = aaf_nom.dropna(axis='columns', how='all')  #, inplace=True)
        aaf_nom.reset_index(inplace=True, drop=True)
        if 'euclidean_dist_nr' in aaf_nom.columns and 'exp_euclidean_dist_nr' in aaf.columns:  # aaf in prev state
            # final_aaf = aaf[aaf['euclidean_dist_nr'] == aaf['euclidean_dist_nr'].min()]
            # final_aaf.reset_index(inplace=True, drop=True)
            # final_aaf.to_csv(final_aaf_fn, index=False)

            final_aaf_nom = aaf_nom[aaf_nom['euclidean_dist_nr'] == aaf_nom['euclidean_dist_nr'].min()]
            final_aaf_exp = aaf[aaf['exp_euclidean_dist_nr'] == aaf['exp_euclidean_dist_nr'].min()]
            final_aaf_nom.reset_index(inplace=True, drop=True)
            final_aaf_exp.reset_index(inplace=True, drop=True)

            log.info('final_aaf_nom.empty ' + str(final_aaf_nom.empty))
            log.info('final_aaf_exp.empty ' + str(final_aaf_exp.empty))

            # print("final_aaf_nom['ref_mode_e2e'] == final_aaf_nom['alt_mode_e2e']",
            #       final_aaf_nom['ref_mode_e2e'] == final_aaf_nom['alt_mode_e2e'],
            #       final_aaf_nom['ref_mode_e2e'].tolist(), final_aaf_nom['alt_mode_e2e'].tolist())
            # print("final_aaf_exp['exp_ref_mode_e2e'] == final_aaf_exp['exp_alt_mode_e2e']",
            #       final_aaf_exp['exp_ref_mode_e2e'] == final_aaf_exp['exp_alt_mode_e2e'],
            #       final_aaf_exp['exp_ref_mode_e2e'].tolist(), final_aaf_exp['exp_alt_mode_e2e'].tolist())

            pd.concat([final_aaf_nom, final_aaf_exp]).to_csv(final_aaf_fn, index=False)

            if not final_aaf_nom.empty:
                return_nom = final_aaf_nom.at[0, 'AAF']
                log_final_nom = str(return_nom) + ' for the nominal formula with the min metric value ' + \
                                str(final_aaf_nom.at[0, 'euclidean_dist_nr'])
            else:
                return_nom = 0
                log_final_nom = str(return_nom) + ' for the nominal formula'
            if not final_aaf_exp.empty:
                return_exp = final_aaf_exp.at[0, 'exp_AAF']
                log_final_exp = str(return_exp) + ' for the expanded formula with the max metric value ' + \
                                str(final_aaf_exp.at[0, 'exp_euclidean_dist_nr'])
            else:
                return_exp = 0
                log_final_exp = str(return_exp) + ' for the expanded formula'
            log.info('Final AAF for one pair of k-mers is ' + log_final_nom + ' and ' + log_final_exp)
            log.info('The supporting data is written to ' + final_aaf_fn)

            fa_8_te = time.time() - fa_8_ts
            return True, False, [return_nom, return_exp], \
                   [final_aaf_nom.at[0, 'ref_seq_ext'], final_aaf_nom.at[0, 'alt_seq_ext']], \
                   [final_aaf_nom.at[0, 'euclidean_dist_nr'], final_aaf_exp.at[0, 'exp_euclidean_dist_nr']], \
                   [[final_aaf_nom.at[0, 'ref_mode_e2e'], final_aaf_nom.at[0, 'alt_mode_e2e']],
                    [final_aaf_exp.at[0, 'exp_ref_mode_e2e'], final_aaf_exp.at[0, 'exp_alt_mode_e2e']]], \
                       [fa_1_te, fa_2_te, fa_3_te,
                        statistics.mean(fa_4_mean), statistics.mean(fa_5_mean),
                        fa_6_te, fa_7_te, fa_8_te]
        elif 'euclidean_dist_nr' in aaf_nom.columns and 'exp_euclidean_dist_nr' not in aaf.columns:
            final_aaf = aaf_nom[aaf_nom['euclidean_dist_nr'] == aaf_nom['euclidean_dist_nr'].min()]
            final_aaf.reset_index(inplace=True, drop=True)
            final_aaf.to_csv(final_aaf_fn, index=False)

            print('final_aaf.empty', final_aaf.empty)
            print('aaf.empty', aaf.empty)

            # print("final_aaf['ref_mode_e2e'] == final_aaf['alt_mode_e2e']",
            #       final_aaf['ref_mode_e2e'] == final_aaf['alt_mode_e2e'],
            #       final_aaf['ref_mode_e2e'].tolist(), final_aaf['alt_mode_e2e'].tolist())

            log_final_nom = str(final_aaf.at[0, 'AAF']) + ' for the nominal formula with the min metric value ' + \
                            str(final_aaf.at[0, 'euclidean_dist_nr'])
            log.info('Final AAF for one pair of k-mers is ' + log_final_nom)
            log.info('The supporting data is written to ' + final_aaf_fn)
            if not final_aaf.empty:
                return_nom = final_aaf.at[0, 'AAF']
            else:
                return_nom = 0
            if not aaf.empty:
                return_exp = aaf.at[0, 'exp_AAF']
            else:
                return_exp = 0

            fa_8_te = time.time() - fa_8_ts
            return True, False, [return_nom, return_exp], \
                   [final_aaf.at[0, 'ref_seq_ext'], final_aaf.at[0, 'alt_seq_ext']], \
                   [final_aaf.at[0, 'euclidean_dist_nr'], None], \
                   [[final_aaf.at[0, 'ref_mode_e2e'], final_aaf.at[0, 'alt_mode_e2e']], None], \
                       [fa_1_te, fa_2_te, fa_3_te,
                        statistics.mean(fa_4_mean), statistics.mean(fa_5_mean),
                        fa_6_te, fa_7_te, fa_8_te]
        elif 'euclidean_dist_nr' not in aaf_nom.columns and 'exp_euclidean_dist_nr' in aaf.columns:
            final_aaf = aaf[aaf['exp_euclidean_dist_nr'] == aaf['exp_euclidean_dist_nr'].min()]
            final_aaf.reset_index(inplace=True, drop=True)
            final_aaf.to_csv(final_aaf_fn, index=False)

            print('final_aaf.empty', final_aaf.empty)

            # print("final_aaf['exp_ref_mode_e2e'] == final_aaf['exp_alt_mode_e2e']",
            #       final_aaf['exp_ref_mode_e2e'] == final_aaf['exp_alt_mode_e2e'],
            #       final_aaf['exp_ref_mode_e2e'].tolist(), final_aaf['exp_alt_mode_e2e'].tolist())

            log_final_exp = str(final_aaf.at[0, 'exp_AAF']) + ' for the expanded formula with the max metric value ' + \
                            str(final_aaf.at[0, 'exp_euclidean_dist_nr'])
            log.info('Final AAF for one pair of k-mers is ' + log_final_exp)
            log.info('The supporting data is written to ' + final_aaf_fn)
            if not aaf_nom.empty:
                return_nom = aaf_nom.at[0, 'AAF']
            else:
                return_nom = 0
            if not final_aaf.empty:
                return_exp = final_aaf.at[0, 'exp_AAF']
            else:
                return_exp = 0

            fa_8_te = time.time() - fa_8_ts
            return True, False, [return_nom, return_exp], \
                   [final_aaf.at[0, 'ref_seq_ext'], final_aaf.at[0, 'alt_seq_ext']], \
                   [None, final_aaf.at[0, 'exp_euclidean_dist_nr']],                    \
                   [None, [final_aaf.at[0, 'exp_ref_mode_e2e'], final_aaf.at[0, 'exp_alt_mode_e2e']]], \
                       [fa_1_te, fa_2_te, fa_3_te,
                        statistics.mean(fa_4_mean), statistics.mean(fa_5_mean),
                        fa_6_te, fa_7_te, fa_8_te]
            # not zero but first or by min_dist
        else:
            # TODO check if correct
            if not aaf_nom.empty:
                log.info('No metric exists, selecting the first calculated AAF')
                log.info('Final AAF for one pair of k-mers is ' + str(aaf_nom.at[0, 'AAF']))
                aaf_nom.iloc[0:1, :].to_csv(final_aaf_fn, index=False)
                fa_8_te = time.time() - fa_8_ts
                return True, False, [aaf_nom.at[0, 'AAF'], aaf_nom.at[0, 'exp_AAF']], [-1, -1], [None, None], [None, None], \
                       [fa_1_te, fa_2_te, fa_3_te,
                        statistics.mean(fa_4_mean), statistics.mean(fa_5_mean),
                        fa_6_te, fa_7_te, fa_8_te]
                       # [aaf_nom.at[0, 'ref_seq_ext'], aaf_nom.at[0, 'alt_seq_ext']], \
                       # [aaf_nom.at[0, 'euclidean_dist_nr'], aaf_nom.at[0, 'exp_euclidean_dist_nr']], \

            else:
                log.info('No metric exists, aaf_nom.empty = ' + str(aaf_nom.empty))
                log.info('Final AAF for one pair of k-mers is ' + str(0))
                fa_8_te = time.time() - fa_8_ts
                return True, False, [0, 0], [-1, -1], [None, None], [None, None], \
                       [fa_1_te, fa_2_te, fa_3_te,
                        statistics.mean(fa_4_mean), statistics.mean(fa_5_mean),
                        fa_6_te, fa_7_te, fa_8_te]


if __name__ == '__main__':
    # run()

    output_dir, sequence_json, alignment_json, lev_table_ref = run_table()
    lev_table_alt = run_table()[-1]
    find_aaf(output_dir, alignment_json, lev_table_ref, lev_table_alt)
