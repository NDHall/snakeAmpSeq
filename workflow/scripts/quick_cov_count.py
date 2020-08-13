




def get_cov_quickcount(file, out,sample):
    qcount = {}
    with open(file=file) as f:
        covf = f.read().rstrip()
        covf = covf.split('\n')
        for line in covf:
            contig, start,stop,site,orientation,count = line.split('\t')
            count = int(count)
            if site not in qcount:
                qcount[site] = count
            else:
                qcount[site] += count
    # sort the keys so we get a consistent output
    qcount_keys = list(qcount.keys())
    qcount_keys.sort()
    with open(out,'w') as fout:
        for qcount_key in qcount_keys:
            fout.write('{sample}\t{site}\t{count}\n'.format(
                sample=sample,
                site=qcount_key,
                count=qcount[qcount_key]
            ))




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
    1.) coverage calculations from bedtools (file)
      sample:
      <contig  start stop site_name orient. count> 
      tua_inf6	519	522	tuaValIlePhe202	+	0
      tua_inf6	630	633	tuaThr239	+	24483
      tua_inf6	641	644	tua243Arg	+	24574
      tua_inf6	717	720	tuaMet268	+	24581
    
    2.) output location
    3.) sample id for matching in later files
It returns summed coverage per site_name. For example you may have multiple contigs
that contain same sites.
    sample out:
    <sample         site_name       aggregated_count>
    JP01-Eli-1      tuaLeu136       0
    JP01-Eli-1      tuaMet268       57736
    JP01-Eli-1      tuaThr239       57635
    JP01-Eli-1      tuaValIlePhe202 0

   

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
    required.add_argument('-f,--file',
                          help='in file a bed coverage file',
                          required=True,
                          dest='file',
                          )
    required.add_argument('-o,--out',
                          required=True,
                          help='output file',
                          dest='out')
    optional.add_argument('-s,--sample_name',
                          help='name of sample being proccessed there should be one sample sample per vcf',
                          dest='sample_name',
                          required=False,
                          default=None)

    argv = parser.parse_args()
    if argv.sample_name is None:
       # strip out path variables and remove file extenstion
       argv.sample_name = '.'.join(argv.file.split('/')[-1].split('.')[:-1])
    get_cov_quickcount(file=argv.file,
                       out=argv.out,
                       sample=argv.sample_name)

