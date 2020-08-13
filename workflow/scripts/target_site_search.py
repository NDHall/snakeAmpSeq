import snp_call_functions
import git
try:
    repo = git.Repo(search_parent_directories=True)
    sha = repo.head.object.hexsha
except:
    sha = 'githash unavailable'
__version__ = 1.00
__author__ = 'Nathan D. Hall'
__updated__ = 'December 10, 2019'



des = """\
This script takes as input a 
    1.) vcf of snps 
    2.) an inframe reference fasta against which all snps are called
    3.) A 5 column bed file of target sites.
        CHROM\tSTART\tSTOP\tname<gene+Ref amino acid+Ref AA position>\t strand
    3.) sample name
    4.) prefix for output files.
It returns
    1.) A file of all SNPs that are translated
        outprefix+.bed
        example:
        #CHROM CODON_START CODON_STOP  POS CODON_POS   REF REF_CODON   REF_AA  ALT ALT_CODON   ALT_AA  SvN QUAL    FILTER  SAMPLE
        epsps_sup	219	222	220	1	A	ACT	Thr	G	GCT	Ala	N	45.1344	None	JP01
        epsps_sup	225	228	228	3	T	GGT	Gly	C	GGC	Gly	S	59	None	JP01

    2.) a list of all snps that overlap target bed sites 
        outprefix_target_snps.txt
        example:
        #sample	named_substitution	ref_chrom	pos_on_ref	codon_pos	ref	alt
        JP01	epspsPro106Ala	epsps_sup	358	1	C	G

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
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = textwrap.dedent(des),
    epilog=textwrap.dedent(epi))

    parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-v,--vcf_file',
                        help='vcf file containing snps',
                        required=True,
                        dest='vcf_file',
                        )
    required.add_argument('-f,--fasta',
                        help='inframe fasta ref used for snp calling',
                        required=True,
                        dest='fasta_file',
                        )
    required.add_argument('-b,--target_bed',
                        required=True,
                        help='5 column bed file with name column consisting of '
                             '<gene name + Ref. AA + Ref. Codon Pos.>',
                        dest='target_bed')
    required.add_argument('-o,--out_prefix',
                        required=True,
                        help='prefix for 2 files written. '
                             'file1: <prefix>.bed contains translation for all snps and refs'
                             'file2: <prefix>_target_snps.txt returns all snps that occur within target_sites',
                        dest='out_prefix')
    optional.add_argument('-s,--sample_name',
                        help='name of sample being proccessed there should be one sample sample per vcf',
                        dest='sample_name',
                        required=False,
                        default=None)
    optional.add_argument('--qual',
                        required=False,
                        help="minimum acceptable vcf quality",
                        dest='qual',
                        default=20)
    optional.add_argument('--filter',
                         required=False,
                         help='filter flag to include',
                         default='PASS',
                         dest='filt')
    argv = parser.parse_args()

    if argv.sample_name is None:
       # strip out path variables and remove file extenstion
       argv.sample_name = '.'.join(argv.vcf_file.split('/')[-1].split('.')[::-1])
    if argv.filt.lower() == 'none':
        #PyVcf4.1 returns a None type for No filter.
        argv.filt = None

    snp_call_functions.target_site_search(
        vcf_path=argv.vcf_file,
        ref_fasta=argv.fasta_file,
        out_prefix=argv.out_prefix,
        sample=argv.sample_name,
        target_sites=argv.target_bed,
        qual=argv.qual,
        filt=argv.filt
   )

