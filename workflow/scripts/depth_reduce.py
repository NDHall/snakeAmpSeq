"""
In ampseq runs it is likely that pcr duplicates will occur
here we are doing 2 things.

1. assuming that all valid reads will have been marked as
   pcr duplicates by picard tools. This is because we are
   sequening at extreme depths.

2. limiting number of identical mapped reads, and requiring
   some minimum number of identical reads.

"""
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




def read_in_sam(samfile,samout,limit=20,min=10):
    start_coor = 146
    chrom = 'epsps_inf'
    sam = pysam.AlignmentFile(samfile, 'rb')
    reduced_bam = pysam.AlignmentFile(samout, "wb", template=sam)
    read_dict = {}
    for read in sam.fetch() :
        if read.is_paired :
            if start_coor != read.reference_start or chrom != read.reference_name:
                """
                Now write out all reads that above the minimum threshold.
                This will limit the output to tractable/sane numbers of reads.
                We could use PCR duplicates marked in Picard, or we don't have too. 
                """
                if start_coor is not None and chrom is not None:
                    for key in read_dict:
                        if len(read_dict[key]) >= min:
                            for read in read_dict[key]:
                                reduced_bam.write(read)

                """
                Reset everything here.
                """
                start_coor = read.reference_start
                chrom = read.reference_name
                read_dict = {}
                read_dict[read.query_alignment_sequence] =[read]

            elif start_coor == read.reference_start and chrom == read.reference_name:
                if read.query_alignment_sequence not in read_dict:
                    read_dict[read.query_alignment_sequence] =[]
                if len(read_dict[read.query_alignment_sequence]) < limit:
                    read_dict[read.query_alignment_sequence].append(read)




    reduced_bam.close()
    sam.close()

des = """\
In ampseq runs it is likely that pcr duplicates will occur
here we are doing 2 things.

1. assuming that all valid reads will have been marked as
   pcr duplicates by picard tools. This is because we are
   sequening at extreme depths.

2. limiting number of identical mapped reads, and requiring
   some minimum number of identical reads.
   
   
This is a fairly simple approach to reducing
pcr duplicates. We are trying to make all 
legitimate snps available to bcftools. In
Some ampSeq files we find that the randomness
of pcr and bias causes the signal to be washed
out.
This takes a sorted bam file as input and returns 
a bam with reduced reads for calling variants.

"""


epi="""\
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
                          help='sorted bam file',
                          required=True,
                          dest='bam_file',
                          )
    required.add_argument('-o,--out',
                          help='out_file',
                          required=True,
                          dest='out_bam',
                          )
    optional.add_argument('--min',
                          help='minimum number of identical read mapping to be included in reduced file'
                               'default = 10',
                          default=10,
                          required=False,
                          dest='min'
                          )
    optional.add_argument('--max',
                          help='maximum number of identical read mapping to be included in reduced file'
                               'default = 10',
                          default=20,
                          required=False,
                          dest='max'
                          )

    argv = parser.parse_args()
    read_in_sam(argv.bam_file,
                argv.out_bam,
                min=argv.min,
                limit=argv.max)
