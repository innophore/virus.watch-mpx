# import
import os
import sys
from Bio import SeqIO
import pandas as pd
from multiprocessing import Pool

# define the inputs
path = ''  # path where the MSA file is located
in_msa = 'MSA_full.fasta'  # filename of the MSA, header should begin with ">" and end with a "|" followed by the collection date, e.g. "|2022-05-03"
split_date = '20220401'  # the date where the MSA should be split in "before" and "after" as a string in the following format: "20220401"
consensus_path = path + 'CONSENSUS/'  # path where the MSA file is located
orfscan_path = path + 'ORFScan/'  # path where the output files are located


def sort_by_collection_date():
    """Splits the input MSA into separate files (before and after a split date)"""

    print('Sequences in the MSA file are sorted by date...')

    # check if output files already exist
    msa_before = path + in_msa[:-6] + '_before_{}.fasta'.format(split_date)
    msa_after = path + in_msa[:-6] + '_after_{}.fasta'.format(split_date)
    if os.path.isfile(msa_before):
        sys.exit('ERROR: Outfile {} already exists.'.format(msa_before))
    if os.path.isfile(msa_after):
        sys.exit('ERROR: Outfile {} already exists.'.format(msa_after))

    # open the output file for writing
    msa_before_out = open(msa_before, 'a')
    msa_after_out = open(msa_after, 'a')

    # read the input MSA file
    headers = []
    sequences = []
    collection_dates = []
    with open(path + in_msa, mode='r') as handle1:
        for record in SeqIO.parse(handle1, 'fasta'):
            headers.append(record.description)
            sequences.append(record.seq)
            raw_date = record.description.split('|')[-1]
            if raw_date.count('-') == 2:  # year, month and day given
                the_year, the_month, the_day = raw_date.split('-')
                if len(the_month) == 2 and len(the_day) == 2:  # if date = 2021-01-02
                    collection_dates.append(int(raw_date.replace('-', '')))
                elif len(the_month) == 1 and len(the_day) == 1:  # if date = 2021-1-2
                    collection_dates.append(int(raw_date.replace('-', '0')))
                elif len(the_month) == 2 and len(the_day) == 1:  # if date = 2021-01-2
                    collection_dates.append(int(raw_date.replace('-', '')[:-1] + '0' + raw_date[-1]))
                elif len(the_month) == 1 and len(the_day) == 2:  # if date = 2021-1-02
                    collection_dates.append(int(raw_date.replace('-', '')[:-3] + '0' + raw_date.replace('-', '')[-3:]))

            elif raw_date.count('-') == 1:  # only year and month given
                the_year, the_month = raw_date.split('-')
                if len(the_month) == 2:  # if date = 2021-01
                    collection_dates.append(int(raw_date.replace('-', '') + '00'))
                elif len(the_month) == 1 and len(the_day) == 1:  # if date = 2021-1
                    collection_dates.append(int(raw_date.replace('-', '0') + '00'))

            elif raw_date.count('-') == 0:  # only year given
                if len(raw_date) == 4:
                    collection_dates.append(int(raw_date + '0000'))
                else:
                    print('This is an unknown date format: {}. The sequence "{}" will not be inlcuded in the date-specific output files.'.format(raw_date, record.description))
                    collection_dates.append('else')
            else:
                print(
                    'This is an unknown date format: {}. The sequence "{}" will not be inlcuded in the date-specific output files.'.format(
                        raw_date, record.description))
                collection_dates.append('else')

    # identify the sequences before and after the split date
    headers_before = []
    sequences_before = []
    headers_after = []
    sequences_after = []
    for he, se, co in zip(headers, sequences, collection_dates):
        if co < int(split_date) and co != int(split_date[:-4] + '0000') and co != int(split_date[:-2] + '00'):
            headers_before.append(he)
            sequences_before.append(se)
        if co >= int(split_date):
            headers_after.append(he)
            sequences_after.append(se)

    # write specific sequences to output files
    for x, y in zip(headers_before, sequences_before):
        msa_before_out.write('>{}\n{}\n'.format(x, y))
    for x, y in zip(headers_after, sequences_after):
        msa_after_out.write('>{}\n{}\n'.format(x, y))

    msa_before_out.close()
    msa_after_out.close()

    print('Sorting finished!')


def build_consensus_genomes():
    """Builds a consensus sequence from the MSA files"""

    print('Consensus genome sequences are generated...')

    if not os.path.exists(consensus_path):
        os.mkdir(consensus_path)

    msa_files = [path + in_msa, path + in_msa[:-6] + '_before_{}.fasta'.format(split_date), path + in_msa[:-6] + '_after_{}.fasta'.format(split_date)]

    global process

    def process(item):
        # define output file names
        consensus_outfile_csv = consensus_path + item.split('/')[-1][:-6] + '_consensus.csv'
        consensus_outfile_fasta = consensus_path + item.split('/')[-1][:-6] + '_consensus.fasta'
        # check if output files already exist
        if os.path.isfile(consensus_outfile_csv):
            sys.exit('ERROR: Outfile {} already exists.'.format(consensus_outfile_csv))
        elif os.path.isfile(consensus_outfile_fasta):
            sys.exit('ERROR: Outfile {} already exists.'.format(consensus_outfile_fasta))
        # build the consensus
        else:
            # collect all sequences from MSA
            sequences = []
            with open(item, mode='r') as handle:
                for record1 in SeqIO.parse(handle, 'fasta'):
                    sequences.append(record1)
            # check if all sequences in the msa are of the same length
            it = iter(sequences)
            the_len = len(next(it))
            if not all(len(l) == the_len for l in it):
                sys.exit(
                    'ERROR: The sequences in the MSA file {} are not of the same length.'.format(item.split('/')[-1]))
            # find the distribution of bases per position
            length = len(sequences[0])
            count = 0
            position_A = [0] * (length)
            position_T = [0] * (length)
            position_G = [0] * (length)
            position_C = [0] * (length)
            position_N = [0] * (length)
            position_gap = [0] * (length)
            for record in sequences:
                count = count + 1
                for ind_let, letter in enumerate(list(record)):
                    if letter.upper() == 'A':
                        position_A[ind_let] = position_A[ind_let] + 1
                    if letter.upper() == 'T':
                        position_T[ind_let] = position_T[ind_let] + 1
                    if letter.upper() == 'G':
                        position_G[ind_let] = position_G[ind_let] + 1
                    if letter.upper() == 'C':
                        position_C[ind_let] = position_C[ind_let] + 1
                    if letter.upper() == 'N':
                        position_N[ind_let] = position_N[ind_let] + 1
                    if letter == '-':
                        position_gap[ind_let] = position_gap[ind_let] + 1
            perc_A = []
            perc_T = []
            perc_G = []
            perc_C = []
            perc_N = []
            perc_gap = []
            for a, t, g, c, n, ga in zip(position_A, position_T, position_G, position_C, position_N, position_gap):
                perc_A.append(a / count)
                perc_T.append(t / count)
                perc_G.append(g / count)
                perc_C.append(c / count)
                perc_N.append(n / count)
                perc_gap.append(ga / count)

            # save the fasta and csv files
            df_consensus = pd.DataFrame(list(zip(perc_A, perc_T, perc_G, perc_C, perc_N, perc_gap)),
                              columns=['A', 'T', 'G', 'C', 'N', '-'])
            df_consensus['consensus'] = df_consensus.idxmax(axis=1)
            df_consensus.to_csv(consensus_outfile_csv)

            consensus_sequence = ''.join(df_consensus['consensus'].tolist())
            consensus_out_fasta = open(consensus_outfile_fasta, 'w')
            consensus_out_fasta.write('>{}\n{}\n'.format(item.split('/')[-1][:-6], consensus_sequence))
            consensus_out_fasta.close()

    if __name__ == '__main__':
        with Pool(3) as p:
            p.map(process, msa_files)

    print('Consensus genome sequence generation finished!')


def ORFScan_of_consensus_genomes():
    """Scans all available consensus sequences for potential ORFs. Retrieves position, strand and frame of potential ORFs"""

    print('The consensus genomes are scanned for potential ORFs...')

    if not os.path.exists(orfscan_path):
        os.mkdir(orfscan_path)

    for files in os.listdir(consensus_path):  # perform the ORF scan on every consensus sequence
        if files.endswith('.fasta'):

            # check if output files already exist
            outfile = orfscan_path + '{}_ORFs.csv'.format(files.split('/')[-1][:-6])
            if os.path.isfile(outfile):
                sys.exit('ERROR: Outfile {} already exists.'.format(outfile))

            # Scan for potential ORFs
            with open(consensus_path + files) as handle:
                for genome in SeqIO.parse(handle, 'fasta'):
                    possible_ORFs = []
                    starts = []  # position as in the genome sequence (independent on strand or frame)
                    ends = []
                    strands = []
                    frames = []
                    for strand, nucleotide in [(+1, genome.seq.replace('-', ''))]:  # forward strand
                        for readingframe in range(3):  # in all three reading frames
                            length = 3 * ((len(genome.seq.replace('-', '')) - readingframe) // 3)
                            count = 1 + readingframe
                            translated_frame = nucleotide[readingframe:readingframe + length].translate(1).split("*")
                            if len(translated_frame[0]) >= 10:  # handle the first protein in the genome
                                possible_ORFs.append(translated_frame[0].replace('Z', 'X').replace('J', 'X').replace('B', 'X'))
                                strands.append(strand)
                                frames.append(readingframe + 1)
                                starts.append(count)
                                count = count + len(translated_frame[0] * 3)
                                ends.append(count - 1)
                            count = len(translated_frame[0] * 3) + 1 + readingframe
                            for protein1, protein2 in zip(translated_frame, translated_frame[1:]):
                                if len(protein2) < 10:
                                    count = count + 3 + len(protein2 * 3)  # add 3 because of the stop codon before the ORF
                                elif len(protein2) >= 10:
                                    possible_ORFs.append(protein2.replace('Z', 'X').replace('J', 'X').replace('B', 'X'))
                                    strands.append(strand)
                                    frames.append(readingframe + 1)
                                    starts.append(count + 3)  # add 3 because of the stop codon before the ORF
                                    count = count + 3 + len(protein2 * 3)  # add 3 because of the stop codon before the ORF
                                    ends.append(count - 1)
                    for strand, nucleotide in [(-1, genome.seq.reverse_complement().replace('-', ''))]:  # reverse strand
                        for readingframe in range(3):  # in all three reading frames
                            length = 3 * ((len(genome.seq.replace('-', '')) - readingframe) // 3)
                            count = 0 + readingframe
                            translated_frame = nucleotide[readingframe:readingframe + length].translate(1).split("*")
                            if len(translated_frame[0]) >= 10:  # handle the first protein in the genome
                                possible_ORFs.append(translated_frame[0].replace('Z', 'X').replace('J', 'X').replace('B', 'X'))
                                strands.append(strand)
                                frames.append(readingframe + 1)
                                ends.append(len(genome.seq.replace('-', '')) - count)
                                count = count + len(translated_frame[0] * 3)
                                starts.append(len(genome.seq.replace('-', '')) - count + 1)
                            count = len(translated_frame[0] * 3) + readingframe
                            for protein1, protein2 in zip(translated_frame, translated_frame[1:]):
                                if len(protein2) < 10:
                                    count = count + 3 + len(protein2 * 3)  # add 3 because of the stop codon before the ORF
                                elif len(protein2) >= 10:
                                    possible_ORFs.append(protein2.replace('Z', 'X').replace('J', 'X').replace('B', 'X'))
                                    strands.append(strand)
                                    frames.append(readingframe + 1)
                                    ends.append(len(genome.seq.replace('-',
                                                                       '')) - count - 3)  # add 3 because of the stop codon before the ORF
                                    count = count + 3 + len(protein2 * 3)  # add 3 because of the stop codon before the ORF
                                    starts.append(len(genome.seq.replace('-', '')) - count + 1)

                    # Save the output file
                    df_out = pd.DataFrame(list(zip(possible_ORFs, starts, ends, strands, frames)),
                                          columns=['ORF_Protein_Sequence', 'Start', 'End', 'Strand', 'Frame'])
                    df_out.to_csv(outfile)
    
    print('Scan for potential ORFs finished!')     
              

sort_by_collection_date()

build_consensus_genomes()

ORFScan_of_consensus_genomes()

