"""
Set of functions for determining if called variants
happen at target sites and if they are silent substitutions.


"""


#import vcf
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import seq3
from pybedtools import BedTool
import textwrap
from operator import itemgetter


def load_vcf(path):
    vcf_rec = vcf.Reader(open(path,'r'))
    return vcf_rec

def get_fasta_dict(fasta_path,seq_type='fasta'):
    """

    :param fasta_path: path to sequence file, fasta is the default format
    :param seq_type: will take any format that Bio.Seq.IO will take just set seq_type accordingly
    :return: dictionary with rec.id from SEQ.IO as the key. This allows us to extract seqs based on seq.id later on.
    """

    _fasta_dict = {}
    for seq_rec in SeqIO.parse(fasta_path,seq_type):
        assert seq_rec.id not in _fasta_dict,"repeated rec.id {id}, please reformat".format(id=seq_rec.id)
        _fasta_dict[seq_rec.id] = seq_rec
    return _fasta_dict


def get_all_codons(vcf_rec,fasta_dict,three_letter_abbr=True):
    """
    Takes a vcf and ref_fasta_dictionary, and returns a list of Record objects with
    codons and aminoacid translations using the standard alphabet.
    CRITICALLY, this functions assumes that all sequences are inframe.
    It will work well for CDS and handcurated genes.

    :param vcf_rec: vcf iterator obj from PyVCF4.1
    :param fasta_dict: Dictionary of Bio.Seq objs ident
    :return: list of modified Record objects that can used write out a file of SNPs

    """
    ret_list = []

    for rec in vcf_rec:
        print(rec.FILTER)
        if rec.POS % 3 == 0 :
            """
            This indicates snp in the 3rd codon position.
            We need to extract the 2 previous nucleotides from 
            the fasta as well. 
            12[3]
            """
            start = rec.POS - 3
            stop = rec.POS
            codon_pos = 3
        elif rec.POS % 3 == 2:
            """
            2nd Codon position
            1[2]3
            """
            start = rec.POS - 2
            stop = rec.POS + 1
            codon_pos = 2
        elif rec.POS % 3 == 1 :
            """
            1st Codon Position
            [1]23
            
            """
            start = rec.POS - 1
            stop = rec.POS + 2
            codon_pos = 1

        """
        get codons and amino acids for mutations
        
        """

        rec.CODON_POS = codon_pos
        rec.CODON_START = start
        rec.CODON_STOP = stop
        rec.CODON_REF = fasta_dict[rec.CHROM].seq[start:stop]
        rec.CODON_ALT = []
        rec.AMINO_REF = rec.CODON_REF.translate()
        if three_letter_abbr is True:
            rec.AMINO_REF = seq3(rec.AMINO_REF)
        rec.AMINO_ALT = []
        for a in rec.ALT :
            """
            ALT alleles come in a list, so alt_codons must follow suit.
            Here we are able to use codon_pos -1 to adjust for python string
            slicing, to replace the ref with the ALT
            """


            alt_codon = list(rec.CODON_REF)
            alt_codon[codon_pos - 1] = a
            alt_codon_str = []
            for nuc in alt_codon:
                alt_codon_str.append(str(nuc))
            del alt_codon # going to create it again as a Seq obj.
            alt_codon = Seq(''.join(alt_codon_str),)


            rec.CODON_ALT.append(alt_codon)
            if three_letter_abbr is True :
                rec.AMINO_ALT.append(seq3(alt_codon.translate()))
            else:
                rec.AMINO_ALT.append(seq3(alt_codon.translate()))
            #print(''.join(alt_codon_str), rec.CODON_REF, alt_codon.translate(),rec.AMINO_REF)
        ret_list.append(rec)
    return ret_list
def write_all_snps(records,out,sample,qual=20.0,filt='PASS'):
    """

    :param records: modified Record object from PyVCF4.1. It contains added categories from
                    func get_all_codons(), written above.
    :param out:     the file to which the results are to be written. It overwrites any existing file
                    without checking.
    :param sample: string designating the name of sample or ampSeq run being processed. As a rule,
                   vcfs should already be sorted by sample, using barcodes and samtools.
    :param qual:    VCF field designating minium quality of SNP to be reported
    :param filt:    VCF field. Default assumes that the user has employed some filtering step. This would be
                    the best practice. However, it can be modified to take any value the Filter Field of a VCF
                    will hold.
    :return:    Nothing this a writing function.

    """
    f = open(out,'w')
    f.write('#CHROM CODON_START CODON_STOP  POS CODON_POS   REF REF_CODON   REF_AA  ALT ALT_CODON'
            '   ALT_AA  SvN QUAL    FILTER  SAMPLE\n')
    for rec in records:
        counter = 0
        for a in rec.ALT:
            if rec.AMINO_REF == rec.AMINO_ALT[counter]:
                syn_vs_nonsyn = 'S'
            else:
                syn_vs_nonsyn = 'N'
            out_string = [str(rec.CHROM),str(rec.CODON_START),str(rec.CODON_STOP),str(rec.POS),str(rec.CODON_POS),
                          str(rec.REF),str(rec.CODON_REF),str(rec.AMINO_REF),
                          str(rec.ALT[counter]), str(rec.CODON_ALT[counter]), str(rec.AMINO_ALT[counter]),
                          str(syn_vs_nonsyn), str(rec.QUAL), str(rec.FILTER), sample+'\n']
            print((rec.QUAL >= qual and( rec.FILTER == filt or len(rec.FILTER) == 0 )),rec.FILTER == filt , rec.FILTER, filt)
            if rec.QUAL >= qual and ( rec.FILTER == filt or len(rec.FILTER) == 0 ):
                """
                Putting in control over filter and quality here.
                defaults assume you have done filtering and the column has passed and that 20 is a sensible
                minimum for quality. 
                """
                f.write('\t'.join(out_string))
            counter += 1
    f.close()

def return_target_snps(target_bed,snp_file):
    """

    :param target_bed: bed file that contains the following fields only
        CHROM   START   STOP    NAME    STRAND
        NAME should be gene+AA+standard AA POS see examples below.
        EXAMPLE:
        tua_inf6	270	273	tuaLeu125	+
        tua_inf6	321	324	tuaLeu136	+
        tua_inf6	519	522	tuaValIlePhe202	+

    :param snp_file: written by function write_all_snps()
    :return: Bedtools intersecttion with flag -wo ( all fields from both intersecting sections, with
    overlapping number of nucleotides supported.
    fields = chrom, codon_start, codon_stop, pos, codon_pos, ref, ref_codon, ref_aa, alt, alt_codon, alt_aa, svn, qual, \
        filter, sample, chrom2, start, stop, name, strand, overlap
    """
    targets = BedTool(target_bed)
    called_snps = BedTool(snp_file)
    intersected_snps = called_snps.intersect(targets,
                                wo=True)
    return intersected_snps
def write_target_snps(target_snps,outpath):
    """

    :param target_snps: takes list of target snps created by function return_target_snps() and returns
                        a report.
    :param outpath: the file to which the report is written. See Sample output below.

    header: #sample	named_substitution	ref_chrom	pos_on_ref	codon_pos	ref	alt
    body:   JP01	epspsPro106Ala	epsps_sup	358	1	C	G



    """

    f = open(outpath,'w')
    f.write('#sample\tnamed_substitution\tref_chrom\tpos_on_ref\tcodon_pos\tref\talt\n')
    for line in target_snps:
        #unpack line into individual variables
        #It is a bit unweildy, but a good quick fix here.
        chrom, codon_start, codon_stop, pos, codon_pos, ref, ref_codon, ref_aa, alt, alt_codon, alt_aa, svn, qual, \
        filter, sample, chrom2, start, stop, name, strand, overlap = line
        outstr = [
            sample,
            name+alt_aa,
            chrom,
            pos,
            codon_pos,
            ref,
            alt+'\n'
        ]
        f.write('\t'.join(outstr))
    f.close()


def target_site_search(vcf_path,ref_fasta,out_prefix,sample,target_sites, qual=0,filt='PASS'):
    vcf_rec = load_vcf(vcf_path)
    seq_dict = get_fasta_dict(ref_fasta)
    records = get_all_codons(vcf_rec,seq_dict)
    write_all_snps(records,out_prefix+'.bed',sample,qual=qual,filt=filt)

    #target_snps() takes file written by write_all_snps()
    target_snps = return_target_snps(target_sites,out_prefix + '.bed')

    write_target_snps(target_snps,out_prefix +'_target_snps.txt' )

"""
We need to also determine if sites are covered so that they can be ruled out.

"""


def determine_coverage(bedfile,bamfile,target_site_mutations):

    #bams and bed through bedtools
    target_sites = BedTool(bedfile)
    bambam = BedTool(bamfile)
    bambed = bambam.bam_to_bed().merge()
    covered_ts = target_sites.intersect(bambed)
    non_covered_ts = target_sites.subtract(bambed)







    site_summary = {}

    """
    Get states for sites associated with each chrom.

    """



    for bed,state in zip([non_covered_ts,covered_ts],['no_coverage','no_mutation']):
        for line in bed:

            # break out parts of bed.
            # bed must have 5 columns.

            chrom, start,stop,name,strand = line
            if chrom not in site_summary:
                site_summary[chrom] = {}
            assert name not in site_summary[chrom],textwrap.dedent("""\
            
            {n} is a repeated site name on {c}.
            All site names must be unique.
            Please modify the name in your
            target site bed file. Sites should
            be named
    
            <gene><3 letter AA symbol><ref. AA position>
    
            Examples:
            alsGly121
            epspsLys201
            
    
            """.format(n=name,c=chrom))

            site_summary[chrom][name]=[int(start),state] # give start position and False for if its covered by reads.

    #read in scored target_site mutations.
    f = open(target_site_mutations)
    tsm = f.read().rstrip().split('\n')
    if len(tsm) >1 : #only these will have TSMs
        for line in tsm[1::]:
            sample, aa_sub, chrom,pos,codon_pos,ref,alt = line.split('\t')
            key_name = aa_sub[:-3]
            assert chrom in site_summary,textwrap.dedent("""\
            
            {c} appears in target site mutations file:
            {tsm}
            but not in the supplied bed file:
            {b}
            
            """.format(c=chrom,tsm=target_site_mutations, b=bedfile))

            assert key_name in site_summary[chrom],textwrap.dedent("""\
            
            {c} appears to be referenced in target site mutations file:
            {tsm}
            but not in the supplied bed file:
            {b}
            
            It may be that amino acid symbol is incorrectly formatted. 
            this program requires a 3 letter symbol, e.g. Pro for proline.
            The error was produced by slicing the last 3 items off this
            string: {s}
            
            
            """.format(c=key_name,tsm=target_site_mutations, b=bedfile,s=aa_sub))
            site_summary[chrom][key_name][-1] =textwrap.dedent("""\
               {aa}
               DNA:{ref}\u2192{alt}
               codon_pos:{c}
               
               """.format(aa=aa_sub[-9::
                             ], ref=ref,alt=alt,c=codon_pos)
           )
            print(site_summary[chrom][key_name][-1])

    """    for chrom in site_summary:
        chrom_sites = []
        for site in site_summary[chrom]:
            chrom_sites.append(site_summary[chrom][site]+[site])
        # sort by starting position using itemgetter

        sorted_chrom = sorted(chrom_sites,key=itemgetter(0))
        print(sorted_chrom)

    """










if __name__ == '__main__':

    determine_coverage('/home/ndh0004/src/ampSeqPipeline/example/sorted_target_site.bed',
                       '/home/ndh0004/src/ampSeqPipeline/example/JP01-Eli-1.sorted.bam',
                       '/home/ndh0004/src/ampSeqPipeline/example/test2_target_snps.txt')