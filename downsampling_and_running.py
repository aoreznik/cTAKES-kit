import glob
import os
import datetime
import time
import pandas as pd
from misc import input_dialog_downsampling_and_running
from aaf_pipeline_singular import run
import logging

log_fn = 'downsampling_and_running.' + datetime.datetime.now().strftime("%d-%m-%Y.%H-%M-%S") + '.log'
logging.basicConfig(filename=log_fn, format='%(asctime)s %(levelname)s : %(message)s', level=logging.INFO)
# logging.basicConfig(stream=sys.stdout, format='%(asctime)s %(levelname)s : %(message)s', level=logging.INFO)
log = logging.getLogger()


if __name__ == '__main__':
    pivot, run_samtools, input_dir, output_dir, reference, dd = input_dialog_downsampling_and_running()
    pivot = pd.read_table(pivot)
    pivot['k'] = None
    pivot['AF_AOD_DP1000'] = None
    pivot['time_DP1000'] = None


    print('Working table:')
    print(pivot)
    log.info('Working table:')
    log.info(pivot)
    for i in range(len(pivot.index)):
        print('\nRunning:')
        log.info('Running:')
        row = pivot.iloc[i, :]
        samtools_fraction = str(dd / row.loc['DP'])
        samtools_input = glob.glob(input_dir + row.loc['SampleOnly'] + '/raw/*bam')
        if samtools_input:
            samtools_output = output_dir + '/' + row.loc['Panel'] + '.' + row.loc['SampleOnly'] + '.' + \
                              samtools_input[0].split('/')[-1].split('.')[-2] + '.' + \
                              row.loc['variant_name'].replace(':', '-').replace('>', '-') + '.1000.bam'

            command = run_samtools + ' view -s ' + samtools_fraction + ' -@ 6 ' + samtools_input[0] + ' -o ' + samtools_output
            print(command)
            log.info(command)
            command_run = os.system(command)
            print('Result: ' + 'OK' if command_run == 0 else 'FAIL')
            log.info('Result: ' + 'OK' if command_run == 0 else 'FAIL')

            command = run_samtools + ' index ' + samtools_output + ' -@ 6 '
            print(command)
            log.info(command)
            command_run = os.system(command)
            print('Result: ' + 'OK' if command_run == 0 else 'FAIL')
            log.info('Result: ' + 'OK' if command_run == 0 else 'FAIL')

            time_start = time.time()
            # print('\nreference', reference)
            # print('alignment', samtools_output)
            # print('variant', row.loc['variant_name'])
            # print('run_samtools', run_samtools)
            # print('output_dir', output_dir)
            print('Internal call to aaf_pipeline_singular, alternative command:\n', 'python aaf_pipeline_singular.py', reference, samtools_output, '"' + row.loc['variant_name'] + '"', run_samtools, output_dir)
            k, mini_sam_fn, final_value_select, final_value, metric_nom, metric_exp, final_mode_equal = run(reference,
                                                                                                            samtools_output,
                                                                                                            row.loc['variant_name'],
                                                                                                            run_samtools,
                                                                                                            output_dir)
            time_end = time.time() - time_start
            pivot.iat[i, pivot.columns.get_loc('k')] = k
            pivot.iat[i, pivot.columns.get_loc('AF_AOD_DP1000')] = final_value_select
            pivot.iat[i, pivot.columns.get_loc('time_DP1000')] = time_end
            # print(pivot.iloc[i, :])

            print('Removing temporary files...')
            log.info('Removing temporary files...')
            for _ in glob.glob(samtools_output + '*'):
                log.debug('\t' + _)
                os.remove(_)
            print('DONE')
            log.info('DONE')

            if i % 100 == 0:
                pivot.to_csv(output_dir + '/pivot_time.' + str(i) + '.tsv', index=False, sep='\t')
    print(pivot)
    pivot.to_csv(output_dir + '/pivot_time.full.tsv', index=False, sep='\t')
