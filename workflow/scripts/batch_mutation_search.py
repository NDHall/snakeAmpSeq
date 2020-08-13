import herbicideMutationFinder as herb
from Bio import SeqIO
sha = 'githash unavailable'
__version__ = 1.00
__author__ = 'Nathan D. Hall'
__updated__ = 'December 13, 2019'

def batch_search(fasta,out,homeolog=True):
    tsr = {}
    stem = fasta.split('/')[-1]
    for seq_record in SeqIO.parse(fasta, 'fasta'):

        """
        Fasta header should be formated as follows:
        <gene>_<homeolog>_<unique id>
        """
        gene,*args = seq_record.description.split('_')
        # args is a list that contains the rest of header

        if homeolog is True:
            homeo = args[0]
        else:
            homeo = None

        res_muts = herb.batchMutationFinder(str(seq_record.seq),gene)
        if len(res_muts) >0:
            detected =[stem,gene,homeo]
            for pos, aa in zip(res_muts[0::2],res_muts[1::2]):
                    if len(aa) > 0 :
                        if 'positions' in pos:
                            detected.append(pos.replace(' positions',''))
                        elif 'position' in pos:
                            detected.append(pos.replace(' position', ''))
                        # Now get total number of seqs with this mutation
                        total_seqs = int(args[-1].split(':')[-1])
                        for a in aa:
                            for ref,mut,codon in zip(a[0::3],a[1::3],a[2::3]):
                                strout = '{r}->{m}({c})'.format(
                                    r=ref,
                                    m=mut,
                                    c=codon
                                )

                                detected.append(strout)
                                key_detected = '\t'.join(detected)
                            if key_detected not in tsr:

                                tsr[key_detected] = [total_seqs,[seq_record.description]]
                            else:
                                tsr[key_detected][0] += total_seqs
                                tsr[key_detected][1].append(seq_record.description)
    with open(out, 'w') as f:
        for snp in tsr:
           f.write(snp+'\t{}\n'.format(tsr[snp][0]))
           for seq in tsr[snp][1]:
               f.write('##\t{id}\t{snp}\n'.format(id=seq,snp=snp))
        f.close()










des = """\
In ampseq runs it is likely that pcr duplicates will occur
here we are doing 2 things.

1. assuming that all valid reads will have been marked as
   pcr duplicates by picard tools. This is because we are
   sequening at extreme depths.

2. extracting unique copies of each read that occurs above
   some cutoff.

Returns a file with fasta name, gene name homeolog designation, AA, AA substition, count  
POA12_uniq.fasta        als     sup     Trp-574 W->L(TTG)       183
##      als_sup_someread3 count:25 POA12_uniq.fasta        als     sup     Trp-574 W->L(TTG)
##      als_sup_someread2 count:66 POA12_uniq.fasta        als     sup     Trp-574 W->L(TTG)
##      als_sup_someread1 count:12 POA12_uniq.fasta        als     sup     Trp-574 W->L(TTG)

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
    required.add_argument('-f,--fasta_file',
                          help='fasta file',
                          required=True,
                          dest='fasta_file',
                          )
    required.add_argument('-o,--out',
                          help='out file',
                          required=True,
                          dest='out_report',
                          )



    argv = parser.parse_args()
    batch_search(argv.fasta_file, argv.out_report, homeolog=True)







