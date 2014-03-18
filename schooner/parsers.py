import pandas as pd
import numpy as np
import sys
from datetime import datetime

def max_csv(x):
    '''Takes the max of integers separated by commas

    Parameters
    ----------
    x : str
        A string of comma-separated integers

    Returns
    -------
    int
        Maximum of integers
    '''
    return max(map(int, x.split(',')))

def min_csv(x):
    '''Takes the minimum of integers separated by commas

    Parameters
    ----------
    x : str
        A string of comma-separated integers

    Returns
    -------
    int
        Minimum of integers
    '''
    return min(map(int, x.split(',')))


def read_one_miso_summary(filename):
    '''Parse a single MISO summary file and add helpful columns

    Parameters
    ----------
    filename : str
        MISO summary filename

    Returns
    -------
    df : pandas.DataFrame
        A pandas DataFrame of the file, with the addition of these columns:
        1. a copy-paste-able genome location at the end, based on the minimum
           mRNA_starts and maximum mRNA_ends. (df.genome_location)
        2. The difference between df.ci_high and df.ci_low (df.ci_diff)
        3. The left and right halves of the confidence interval, e.g. the right
           half is df.ci_high - df.miso_posterior_mean. (df.ci_left_half and
           df.ci_right_half)
        4. The max of the two left and right confidence interval halves
           (df.ci_halves_max)
    '''
    df = pd.read_table(filename, index_col=0)
    genome_location = pd.DataFrame(
        ['%s:%d-%d' % (chrom, min_csv(starts), max_csv(stops))
         for chrom, starts, stops in zip(df.chrom,
                                         df.mRNA_starts,
                                         df.mRNA_ends)],
        columns=['genome_location'], index=df.index)
    ci_diff = pd.DataFrame(df.ci_high - df.ci_low, columns=['ci_diff'],
                           index=df.index)
    ci_halves = pd.DataFrame(
        {'ci_left_half': (df.ci_high - df.miso_posterior_mean),
         'ci_right_half': (df.miso_posterior_mean - df.ci_low)},
        index=df.index)
    ci_halves_max = pd.DataFrame(ci_halves.max(axis=1),
                                        columns=['ci_halves_max'])
    return pd.concat([df, genome_location, ci_diff, ci_halves,
                      ci_halves_max], axis=1)

def read_sample_info(sample_info_filename):
    """Read a sample info file

    """
    return pd.read_table(sample_info_filename)

def get_miso_summaries(sample_info, ci_halves_max_thresh=0.2,
                       reporting_interval=100):
    """Read a bunch of miso summary files

    Parameters
    ----------
    sample_info : str
        Location of the Schooner-compatible sample info file. Must have a
        column called 'miso-summary_filename'

    ci_halves_max_thresh : float
        Threshold of the maximum confidence interval half to accept for an
        event.

    Returns
    -------
    pandas.DataFrame
        A "tall" dataframe of all the summary files concatenated together,
        on top of each other. All events with ci_halves_max less than the
        threshold have been removed. Also contains all data from the
        sample_info_file.

    """
    dfs = []


    sys.stdout.write(
        'Parsing {} MISO summary files...'.format(sample_info.shape[0]))
    # sample_info_all = pd.read_table(sample_info_filename)

    start_time = datetime.now()
    for i, row in sample_info.iterrows():
        if i % reporting_interval == 0:
            sys.stdout.write('\tparsing {} of {} miso summary files'.format(
                i, sample_info.shape[0]))
        filename = row['miso_summary_filename']
        df = read_one_miso_summary(filename)

        for column_name in row.index:
            if column_name == 'miso_summary_filename':
                continue
            df[column_name] = row[column_name]

        dfs.append(df.reset_index())
    sys.stdout.write('\tDone. Elapsed time: {}'.format(datetime.now() -
                                                    start_time))

    sys.stdout.write('Concatenating and filtering MISO summary files...')
    start_time = datetime.now()
    summary = pd.concat(dfs)
    summary = summary[summary.ci_halves_max <= ci_halves_max_thresh]

    # since we just removed a bunch of rows, reset the index to be simply
    # integers in order, so there's no weird gaps
    summary.index = np.arange(summary.shape[0])

    # Add isoformA and isoformB counts for the two isoforms
    isoform_counts = assigned_counts_to_isoform_counts(summary.assigned_counts)
    summary = summary.join(isoform_counts)
    summary.sort_index(axis=1, inplace=True)
    sys.stdout.write('\tDone. Elapsed time: {}'.format(datetime.now() -
                                                       start_time))
    return summary

def assigned_counts_to_isoform_counts(assigned_counts):
    """Transform an assigned counts column into a dataframe with counts of both isoforms

    """
    int_to_isoform = {0: 'isoformA_counts',
                      1: 'isoformB_counts'}
    assigned_counts = assigned_counts.to_dict()
    isoform_counts = dict((k, dict(
        (int_to_isoform[int(pair.split(':')[0])], int(pair.split(':')[1]))
        for pair in v.split(','))) for k, v in assigned_counts.iteritems())
    isoform_counts = pd.DataFrame.from_dict(isoform_counts, orient='index')
    return isoform_counts