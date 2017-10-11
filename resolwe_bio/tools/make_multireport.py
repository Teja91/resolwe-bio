#!/usr/bin/env python2
"""Generate amplicon multireport."""

import subprocess
import argparse

DECIMALS=2

parser = argparse.ArgumentParser(description="Fill data into tex template file.")
parser.add_argument('--sample', help="Sample name.", nargs='+')
parser.add_argument('--panel', help="Panel name", nargs='+')
parser.add_argument('--covplot', help="Coverage plot.", nargs='+')
parser.add_argument('--covmetrics', help="Coverge metrics", nargs='+')
parser.add_argument('--cov', help="Amplicon coverage", nargs='+')
parser.add_argument('--metrics', help="CollectTargetedPcrMetrics report file.", nargs='+')
parser.add_argument('--vcfgatkhc', help="File with VCF GATK HaplotypeCaller data.", nargs='+')
parser.add_argument('--vcflf', help="File with VCF Lofreq data.", nargs='+')
parser.add_argument('--template', help="Report template file.")
parser.add_argument('--logo', help="Logo.")

def _escape_latex(string):
    """Format normal string to be LaTeX compliant."""
    string = string.replace('\\', '\\\\')
    string = string.replace('&', '\&')
    string = string.replace('%', '\%')
    string = string.replace('$', '\$')
    string = string.replace('#', '\#')
    string = string.replace('_', '\_')
    string = string.replace('{', '\{')
    string = string.replace('}', '\}')
    string = string.replace('~', '\textasciitilde ')
    string = string.replace('^', '\textasciicircum ')
    return string

def _remove_underscore(string):
    """Remove underscore from the LaTeX compliant string and replaces it with white space."""
    return string.replace('\_', ' ')

def _tsv_to_list(table_file, has_header=False, delimiter='\t', pick_columns=None):
    """Parse *.tsv file into list."""
    table_data = []
    header = None
    with open(table_file, 'r') as tfile:
        if has_header:
            header = next(tfile).strip().split(delimiter)
            common_columns = [x for x in pick_columns if x in header]
            if pick_columns:
                # Find indexes of selected columns
                temp_header = {col: i for i, col in enumerate(header)}
                pick_indexes = [temp_header[col] for col in common_columns]
                header = common_columns
        for line in tfile:
            line_content = line.strip().split(delimiter)
            if pick_columns:
                # Select only specific columns and order them as in pick_columns:
                temp_line_content = {i: col for i, col in enumerate(line_content)}
                line_content = [temp_line_content[i] for i in pick_indexes]
            table_data.append(line_content)
    return table_data, header

def cut_to_pieces(string, piece_size):
    """Cut long string into smaller pieces."""
    pieces = []
    while len(string) > piece_size:
        pieces.append(string[:piece_size])
        string = string[piece_size:]
    pieces.append(string)
    return ' '.join(pieces)

def list_to_tex_table(data, header=None, caption=None, long_columns=False):
    """Make a TeX table from python list."""
    lines = []
    column_alingnment = ['l' for h in header]
    if long_columns:
        for col_index in long_columns:
            column_alingnment[col_index] = 'L'

    if caption:
        lines.append('\\captionof{{table}}{{{}}}'.format(caption))

    lines.append('\\noindent')
    lines.append('\\begin{{longtable}}[l] {{ {} }}'.format('  '.join(column_alingnment)))

    if header:
        lines.append('\\rowcolor{darkblue1}')
        lines.append('\\leavevmode\\color{white}\\textbf{' +
                     '}& \\leavevmode\\color{white}\\textbf{'.join(header) + '} \\\\')

    for line in data:
        if long_columns:
            for col_index in long_columns:
                # If hyperlink line, don't do a thing. Otherwise, insert spaces, so that wrapping can happen:
                new_val = line[col_index] if '\href' in line[col_index] else cut_to_pieces(line[col_index], 8)
                line[col_index] = '\\multicolumn{{1}}{{m{{2.3cm}}}}{{{}}}'.format(new_val)

        lines.append(' & '.join(line) + ' \\\\')

    lines.append('\\end{longtable}')
    lines.append('{\n\\addtocounter{table}{-1}}')
    return '\n'.join(lines)

def snp_href(snpid):
    """Create LaTeX hyperlink for given SNP ID."""
    if snpid.startswith('rs'):
        url = 'http://www.ncbi.nlm.nih.gov/snp/?term={}'.format(snpid)
        pass
    elif snpid.startswith('COSM'):
        url = 'http://grch37-cancer.sanger.ac.uk/cosmic/mutation/overview?id={}'.format(snpid.lstrip('COSM'))
    elif snpid.startswith('COSN'):
        url = 'http://cancer.sanger.ac.uk/cosmic/ncv/overview?id={}'.format(snpid.lstrip('COSN'))
    else:
        return snpid
    return '\\href{{{}}}{{{}}}'.format(url, snpid)


def gene_href(gene_name):
    """Create LaTeX hyperlink for given GENE ID."""
    url = 'http://www.ncbi.nlm.nih.gov/gene/?term={}'.format(gene_name)
    return '\\href{{{}}}{{{}}}'.format(url, gene_name)



def parse_target_pcr_metrics(metrics_report):
    """Parse CollectTargetedPcrMetrics report file."""
    with open(metrics_report) as pcr_metrics:
        labels = []

        for line in pcr_metrics:
            if line.startswith('#'):
                continue
            if len(labels) == 0:
                if line.startswith('CUSTOM_AMPLICON_SET'):
                    labels = line.strip().split('\t')
            else:
                values = line.strip().split('\t')
                break

        return dict(zip(labels, values))

def parse_covmetrics(covmetrics_report):
    """Get coverage stats."""
    with open(covmetrics_report) as covmetrics:
        cov_data = covmetrics.readline().strip().split('\t')
        mean_coverage = float(cov_data[0])
        mean20 = float(cov_data[1])
        cov_unif = float(cov_data[2])
        return cov_data, mean_coverage, mean20, cov_unif

def vcf_table_name(vcf_file_name):
    """Format VCF table caption."""
    if 'gatkhc.finalvars' in vcf_file_name.lower():
        return 'GATK HaplotypeCaller variant calls'
    elif 'lf.finalvars' in vcf_file_name.lower():
        return 'Lofreq variant calls'
    else:
        return os.path.basename(vcf_file_name)

def parse_cov(cov_file_name):
    cov_list, _ = _tsv_to_list(cov_file_name)
    amp_numb = list(range(len(cov_list)))
    covered_amplicons = len([1 for line in cov_list if float(line[8]) >= 1])
    return cov_list, amp_numb, covered_amplicons

def samplelist_to_string(samplelist):
    samples=[]
    for item in samplelist:
        samples.append(item[0]+' ('+item[1]+')')
    return ', '.join(samples)

if __name__ == '__main__':
    args = parser.parse_args()
    
    #Make a index list of sorted sample names (alphabetical order)
    indexlist=[i[0] for i in sorted(enumerate(args.sample), key=lambda x:x[1])]

    # Open template and fill it with data:
    with open(args.template, 'r') as template_in, open('multireport.tex', 'wt') as template_out:
        template = template_in.read()
        template = template.replace('{#LOGO#}', args.logo)

        #create dict with metrics & coverage data
        d={}
        for i in range(len(args.sample)):
            d[args.sample[i]]=parse_target_pcr_metrics(args.metrics[i])
            d[args.sample[i]]['cov_data'], d[args.sample[i]]['mean_coverage'], d[args.sample[i]]['mean20'], d[args.sample[i]]['cov_unif']=parse_covmetrics(args.covmetrics[i])

            #read .cov file
            d[args.sample[i]]['cov_list'], d[args.sample[i]]['amp_numb'], d[args.sample[i]]['covered_amplicons']=parse_cov(args.cov[i]) 
        
        #make QC information table
        QCheader=['Sample name', 'Total reads', 'Aligned reads', 'Aligned bases', 'Mean coverage', 'Threshold coverage', 'Coverage uniformity', 'No. of amplicons', 'Amplicons with 100 coverage']
        lines=[]
        for i in indexlist:
            pct_aligned_reads = str('{0:g}'.format(round(float(d[args.sample[i]]['PCT_PF_UQ_READS_ALIGNED']) * 100, DECIMALS)))
            pct_amplified_bases = str('{0:g}'.format(round(float(d[args.sample[i]]['PCT_AMPLIFIED_BASES']) * 100, DECIMALS)))
        
            lines.append([_escape_latex(args.sample[i]), d[args.sample[i]]['TOTAL_READS'], pct_aligned_reads, pct_amplified_bases, str(round(d[args.sample[i]]['mean_coverage'], DECIMALS)), str(round(d[args.sample[i]]['mean20'], DECIMALS)),str(round(d[args.sample[i]]['cov_unif'], DECIMALS)), str(len(d[args.sample[i]]['cov_list'])), str(d[args.sample[i]]['covered_amplicons'])])

        QCinfo=list_to_tex_table(lines, header=QCheader, long_columns=[1, 2, 3, 4, 5, 6, 7, 8])

        template=template.replace('{#QCTABLE#}', QCinfo)

        #Make Amplicons with coverage < 100% table

        cols=['Sample', 'Amplicon', '\% Covered', 'Less than 20\% of mean']
        not_covered=[]
        for i in indexlist:
            not_covered=not_covered+[(_escape_latex(args.sample[i]), _escape_latex(line[4]), '{:.1f}'.format(float(line[8]) * 100), '0') for line in d[args.sample[i]]['cov_list'] if float(line[8]) < 1]
        if not_covered==[]:
            not_covered=[['/', '/']]

        table_text=list_to_tex_table(not_covered, header=cols)
        template = template.replace('{#BAD_AMPLICON_TABLE#}', table_text)

        #Make VCF GATKHC
        table_text=''
        gatkhc_variants={}
        lf_variants={}
        
        header = ['CHROM', 'POS', 'REF', 'ALT', 'AF', 'DP', 'DP4', 'GEN[0].AD', 'SB', 'FS', 'EFF[*].GENE', 'ID']

        for i in indexlist:
            #GATKHC:
            vcf_table_1, common_columns_1 = _tsv_to_list(args.vcfgatkhc[i], has_header=True, pick_columns=header)

            # Escape user inputs:
            common_columns_1 = [_escape_latex(name) for name in common_columns_1]
            ##########caption = _escape_latex(vcf_table_name(vcf_file))
            vcf_table_1 = [[_escape_latex(value) for value in line] for line in vcf_table_1]

            # Insert space between SNP ID's and create hypelinks:
            vcf_table_1 = [line[:-1] + [' '.join(map(snp_href, line[-1].split(';')))] for line in vcf_table_1]

            #Create dict of shared variants
            for line in vcf_table_1:
                variant=line[-2]+'_chr'+line[0]+'_'+line[1]
                if variant in gatkhc_variants.keys():
                    gatkhc_variants[variant].append([args.sample[i], line[4]])
                else: gatkhc_variants[variant]=[[args.sample[i], line[4]]]

            # Create gene hypelinks:
            vcf_table_1 = [line[:-2] + [gene_href(line[-2])] + [line[-1]] for line in vcf_table_1]

            

            table_text += list_to_tex_table(vcf_table_1, header=common_columns_1, caption='caption', long_columns=[2, 3, -1])
            table_text += '\n\\newpage\n'

            #LF:
            vcf_table_2, common_columns_2 = _tsv_to_list(args.vcflf[i], has_header=True, pick_columns=header)

            # Escape user inputs:
            common_columns_2 = [_escape_latex(name) for name in common_columns_2]
            ########caption = _escape_latex(vcf_table_name(vcf_file))
            vcf_table_2 = [[_escape_latex(value) for value in line] for line in vcf_table_2]

            # Insert space between SNP ID's and create hypelinks:
            vcf_table_2 = [line[:-1] + [' '.join(map(snp_href, line[-1].split(';')))] for line in vcf_table_2]

            #Create dict of shared variants
            for line in vcf_table_2:
                variant=line[-2]+'_chr'+line[0]+'_'+line[1]
                if variant in lf_variants.keys():
                    lf_variants[variant].append([args.sample[i], line[4]])
                else: lf_variants[variant]=[[args.sample[i], line[4]]]

            # Create gene hypelinks:
            vcf_table_2 = [line[:-2] + [gene_href(line[-2])] + [line[-1]] for line in vcf_table_2]

            table_text += list_to_tex_table(vcf_table_2, header=common_columns_2, caption='caption', long_columns=[2, 3, -1])
            table_text += '\n\\newpage\n'

        #Create shared variants tables
        header_shared=['Variant ID', 'Samples sharing this variant (allele frequence)']
        data_gatkhc=[]
        for k, v in gatkhc_variants.items():
            data_gatkhc.append([_escape_latex(k), _escape_latex(samplelist_to_string(v))])
        gatkhc_shared=list_to_tex_table(data_gatkhc, header=header_shared)

        data_lf=[]
        for k, v in lf_variants.items():
            data_lf.append([_escape_latex(k), _escape_latex(samplelist_to_string(v))])
        lf_shared=list_to_tex_table(data_lf, header=header_shared)

        template = template.replace('{#GATKHC_SHARED#}', gatkhc_shared)
        template = template.replace('{#LF_SHARED#}', lf_shared)

        template = template.replace('{#VCF_TABLES#}', table_text)
        
        # Write template to 'report.tex'
        template_out.write(template)

    # Generate PDF:
    args = ['pdflatex', '-interaction=nonstopmode', 'multireport.tex']
    subprocess.call(args)
    subprocess.check_call(args)

    print("Done.")