


def combine_reports(cov_report,mut,out):
    cov_dict = {}
    with open(cov_report) as f:
        cov = f.read().rstrip().split('\n')
        for line in cov:
            sample_name, site,coverage = line.split('\t')
            sample_name = convert_sample(sample_name)
            if sample_name not in cov_dict:
                cov_dict[sample_name] = {}
            cov_dict[sample_name][site] = {'coverage':int(coverage),
                                            'mutation': False }
    del f
    with open(mut,'r') as f:
        mut = f.read().rstrip().split('\n')
        for line in mut:
            sample, gene, homeo, mutation = line.split('\t')
            sample_name = convert_sample(sample_name)
            mutation = mutation.split(' ')[0].replace('-','')
            site_name = gene +  mutation
            print(cov_dict.keys())
            assert sample_name in cov_dict, 'sample_name : {}'.format(sample_name)
            assert site_name in cov_dict[sample_name]
            cov_dict[sample_name][site_name]['mutation'] = True
    cov_keys = list(cov_dict.keys())
    cov_keys.sort()
    with open(out, 'w') as fout:
        fout.write('sample\tsite\tcoverage\tmutation\n')
        for cov_key in cov_keys:
            site_keys = list(cov_dict[cov_key].keys())
            site_keys.sort()
            for site_key in site_keys:
                fout.write('{sample}\t{site}\t{count}\t{mutation}\n'.format(
                sample=cov_key,
                site=site_key,
                count=cov_dict[cov_key][site_key]['coverage'],
                mutation=cov_dict[cov_key][site_key]['mutation']
            ))








def convert_sample(sample):
    """

    :param sample: standardizes sample variable for this pipeline only!
    :return: stripped sample that corresponds to snakemake pipeline.
    """
    print(sample)
    if sample.split('_')[-1] == 'uniq.fasta':
        ret = '_'.join(sample.split('_')[:-1])
    elif sample.split('.')[-1] == 'cov':
        ret = '.'.join(sample.split('.')[:-1])
    else:
        ret = sample
    print(ret)
    return sample






if __name__ == '__main__':
    combine_reports(cov_report='/home/ndh0004/Dropbox/ampSeq/7.tally/POA01.tsv',
                    mut='/home/ndh0004/Dropbox/ampSeq/4.tsr_report/POA01_tsr_report.tsv',
                    out='/home/ndh0004/Dropbox/ampSeq/6.cov_per_site/POA01.report'
    )