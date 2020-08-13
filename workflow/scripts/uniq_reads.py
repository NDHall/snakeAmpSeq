import git
import pysam
try:
    repo = git.Repo(search_parent_directories=True)
    sha = repo.head.object.hexsha
except:
    sha = 'githash unavailable'
__version__ = 1.00
__author__ = 'Nathan D. Hall'
__updated__ = 'December 13, 2019'



def update_coords(read):
    start_coor = read.reference_start
    chrom = read.reference_name
    read_dict = {}
    read_dict = add_new_read(read,read_dict)
    return [start_coor, chrom, read_dict]

def add_new_read(read,read_dict):
    """

    :param read:
    :param read_dict:
    :return:
    """
    read_dict[read.query_sequence] = {'count': 1, 'read': read}
    return read_dict

def move_dict_to_list(read_dict,min,out_seqs):

    for key in read_dict:

        if read_dict[key]['count'] >= min:

            out_seqs.append([read_dict[key]['count'], read_dict[key]['read']])
    return out_seqs


def read_in_sam(samfile,out,min=5,start_coor=0,chrom=None):
    sam = pysam.AlignmentFile(samfile, 'rb')
    read_dict = {}
    out_seqs = []
    counter = 0
    for read in sam.fetch() :
        if read.is_paired :

            if chrom is None:
                """
                Start it all off
                """
                start_coor, chrom, read_dict = update_coords(read)

            elif start_coor != read.reference_start or chrom != read.reference_name:
                """
                Now write out all reads that above the minimum threshold.
                This will limit the output to tractable/sane numbers of reads.
                We could use PCR duplicates marked in Picard, or we don't have too. 
                """

                out_seqs = move_dict_to_list(read_dict,min,out_seqs)


                """
                Reset everything here.
                """
                start_coor,chrom,read_dict = update_coords(read)

            elif start_coor == read.reference_start and chrom == read.reference_name:
                if read.query_sequence not in read_dict:
                    read_dict = add_new_read(read, read_dict)
                else:
                    read_dict[read.query_sequence]['count'] += 1
    out_seqs = move_dict_to_list(read_dict, min, out_seqs)
    fout = open(out,'w')
    for seq in out_seqs:
        count,read = seq
        fout.write('>{ch}_{id} count:{c}\n'.format(
            ch=read.reference_name,
            id=read.query_name,
            c=count
        ))
        fout.write(read.query_sequence+"\n")
    fout.close()


des = """\
In ampseq runs it is likely that pcr duplicates will occur
here we are doing 2 things.

1. assuming that all valid reads will have been marked as
   pcr duplicates by picard tools. This is because we are
   sequening at extreme depths.

2. extracting unique copies of each read that occurs above
   some cutoff.

"""

epi = """\
version:{v}
author: {a}
updated: {u}
written for Python 3.5
githash: {g}

""".format(v=__version__,
           a=__author__,
           u=__updated__,
           g=sha)

if __name__ == '__main__':
    import argparse
    import textwrap

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(des),
        epilog=textwrap.dedent(epi))

    parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-b,--bam_file',
                          help='list of sorted bam files',
                          required=True,
                          dest='bam_file',
                          )
    required.add_argument('-o,--out',
                          help='out file',
                          required=True,
                          dest='out_bam',
                          )
    optional.add_argument('--min',
                          help='minimum number of identical read mapping to be included in reduced file'
                               'default = 5',
                          default=5,
                          required=False,
                          dest='min'
                          )


    argv = parser.parse_args()
    read_in_sam(argv.bam_file,
                argv.out_bam,
                min=argv.min
                )

