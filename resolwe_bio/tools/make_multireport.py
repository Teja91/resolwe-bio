#!/usr/bin/env python3
"""Generate amplicon multireport."""

import subprocess
import argparse
import numpy as np
import re

from bokeh import palettes
from bokeh.plotting import figure, save, output_file
from bokeh.models import HoverTool, ColumnDataSource, LinearColorMapper

DECIMALS = 2

parser = argparse.ArgumentParser(description="Fill data into tex template file.")
parser.add_argument('--sample', help="Sample name.", nargs='+')
parser.add_argument('--covmetrics', help="Coverge metrics", nargs='+')
parser.add_argument('--cov', help="Amplicon coverage", nargs='+')
parser.add_argument('--metrics', help="CollectTargetedPcrMetrics report file.", nargs='+')
parser.add_argument('--vcfgatkhc', help="File with VCF GATK HaplotypeCaller data.", nargs='+')
parser.add_argument('--vcflf', help="File with VCF Lofreq data.", nargs='+')
parser.add_argument('--meancov', help="Mean amplicon coverage.", nargs='+')
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
    """Remove underscore from LaTeX compliant string and replace it with white space."""
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
                line_content = [temp_line_content.get(i, '') for i in pick_indexes]
            table_data.append(line_content)
    return table_data, header


def cut_to_pieces(string, piece_size):
    """Cut long string into smaller pieces."""
    assert isinstance(string, str)
    string_size = len(string)
    pieces = [string[i: i + piece_size] for i in range(0, string_size, piece_size)]
    return ' '.join(pieces)


def list_to_tex_table(data, header=None, caption=None, long_columns=False, wide_columns=False, uncut_columns=[]):
    """Make a TeX table from python list."""
    lines = []
    column_alingnment = ['l' for h in header]
    if long_columns:
        for col_index in long_columns:
            column_alingnment[col_index] = 'L'
    if wide_columns:
        for col_index in wide_columns:
            column_alingnment[col_index] = 'W'

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
                # Insert spaces, so that wrapping can happen, for columns
                # without hyperlinks or other specified uncut columns.
                if ('\href' in line[col_index] or col_index in uncut_columns):
                    new_val = line[col_index]
                else:
                    new_val = cut_to_pieces(line[col_index], 8)
                line[col_index] = '\\multicolumn{{1}}{{m{{2.3cm}}}}{{{}}}'.format(new_val)

        if wide_columns:
            for col_index in wide_columns:
                line[col_index] = '\\multicolumn{{1}}{{m{{18cm}}}}{{{}}}'.format(line[col_index])

        lines.append(' & '.join(line) + ' \\\\')

    lines.append('\\end{longtable}')
    lines.append('{\n\\addtocounter{table}{-1}}')
    return '\n'.join(lines)


def snp_href(snpid):
    """Create LaTeX hyperlink for given SNP ID."""
    if snpid.startswith('rs'):
        url = 'http://www.ncbi.nlm.nih.gov/snp/?term={}'.format(snpid)
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


def parse_cov(cov_file_name):
    """Parse *.cov file."""
    cov_list, _ = _tsv_to_list(cov_file_name)
    amp_numb = list(range(len(cov_list)))
    covered_amplicons = len([1 for line in cov_list if float(line[8]) >= 1])
    return cov_list, amp_numb, covered_amplicons


def parse_mean_covd(covd_file_name):
    """Parse *.covd file."""
    covd_list, _ = _tsv_to_list(covd_file_name)
    mean_amp_cov = {}
    for line in covd_list:
        mean_amp_cov[line[0]] = line[1]
    return mean_amp_cov


def samplelist_to_string(samplelist):
    """Converts list of lists to a string, suitable for shared variants tables."""
    samples = []
    for item in samplelist:
        samples.append(item[0] + ' (' + str(round(float(item[1]), 3)) + ')')
    return ', '.join(samples)


def produce_warning(coverage, mean_20):
    """Produces a string with warning if coverage is less than 20% of mean."""
    if coverage < mean_20:
        return 'YES'
    else:
        return ' '


def make_heatmap(samples, variant_dict, fig_name):
    """Creates a heatmap of Samples x Variants."""
    # Prepare data: make sample and variant list and NumPy array of AF.
    x_names = sorted(list(variant_dict.keys()))
    y_names = sorted([x.strip('.bam') for x in samples])
    data = np.zeros((len(x_names), len(y_names),),)
    for k, v in variant_dict.items():
        for item in v:
            variant_index = x_names.index(k)
            sample_index = y_names.index(item[0])
            data[variant_index, sample_index] = round(float(item[1]), 3)
    width = min(len(x_names) * 100, 1100)
    height = min(len(y_names) * 100, 580)

    # Make data into form suitable for bokeh
    x_flat = []
    y_flat = []
    data_flat = []
    for i, var in enumerate(x_names):
        for j, sam in enumerate(y_names):
            x_flat.append(var)
            y_flat.append(sam)
            data_flat.append(data[i, j])
    source = ColumnDataSource(data=dict(xname=x_flat, yname=y_flat, data_flat=data_flat))

    p = figure(
        title="Shared {} across samples".format(fig_name),
        x_axis_location="above",
        tools="hover,save,pan,box_zoom,wheel_zoom,reset",
        plot_width=width,
        plot_height=height,
        x_range=x_names,
        y_range=y_names,
        toolbar_location='below',
    )
    p.grid.grid_line_color = None
    p.axis.axis_line_color = None
    p.axis.major_tick_line_color = None
    p.axis.major_label_text_font_size = "5pt"
    p.axis.major_label_standoff = 0
    p.xaxis.major_label_orientation = None
    p.xaxis.major_label_orientation = np.pi / 2

    palette = palettes.Blues[9]
    mapper = LinearColorMapper(palette=palette, low=0, high=3)
    p.rect('xname', 'yname', 1, 1, source=source, alpha='data_flat',
           fill_color={'field': 'data_flat', 'transform': mapper}, line_color=None, hover_line_color='black')

    p.select_one(HoverTool).tooltips = [
        ('Sample', '@yname'),
        ('Variant', '@xname'),
        ('Allele Frequency', '@data_flat'),
    ]

    output_file("{}.html".format(fig_name.replace(" ", "")), title=fig_name)
    save(p)

def aa_change(aa_list):
    """Creates Amino Acid Change information."""
    if aa_list:
        aa = aa_list[0]
        match_obj = re.match( r'p\.([A-Za-z]*)[0-9]*([A-Za-z]*)', aa)
        if match_obj and match_obj.group(1) == match_obj.group(2):
            return 'Synon'
        else:
            return aa
    else:
        return aa_list

if __name__ == '__main__':
    args = parser.parse_args()

    # Make a index list of sorted sample names (alphabetical order)
    indexlist = [i[0] for i in sorted(enumerate(args.sample), key=lambda x: x[1])]

    # Open template and fill it with data:
    with open(args.template, 'r') as template_in, open('multireport.tex', 'wt') as template_out:
        template = template_in.read()
        template = template.replace('{#LOGO#}', args.logo)

        # create dict with metrics & coverage data
        d = {}
        for i, sample in enumerate(args.sample):
            d[sample] = parse_target_pcr_metrics(args.metrics[i])
            d[sample]['cov_data'], d[sample]['mean_coverage'], d[sample]['mean20'], d[sample]['cov_unif'] = \
                parse_covmetrics(args.covmetrics[i])

            # read .cov file
            d[sample]['cov_list'], d[sample]['amp_numb'], d[sample]['covered_amplicons'] = parse_cov(args.cov[i])

        # make QC information table
        qc_header = ['Sample name', 'Total \\linebreak reads', 'Aligned \\linebreak reads',
                     'Aligned \\linebreak bases\\footnote{on target}', 'Mean \\linebreak coverage',
                     'Threshold \\linebreak coverage\\footnote{20\\% of mean}',
                     'Coverage \\linebreak uniformity\\footnote{\\% bases covered 20\\% above mean}',
                     'No. of \\linebreak amplicons', 'Amplicons with 100\\% coverage']
        lines = []
        for i in indexlist:

            pct_aligned_reads = str('{0:g}'.format(round(float(d[args.sample[i]]['PCT_PF_UQ_READS_ALIGNED']) * 100, DECIMALS)))
            pct_amplified_bases = str('{0:g}'.format(round(float(d[args.sample[i]]['PCT_AMPLIFIED_BASES']) * 100, DECIMALS)))

            lines.append([_escape_latex(args.sample[i]).strip('.bam'), d[args.sample[i]]['TOTAL_READS'], pct_aligned_reads+'\\%',
                          pct_amplified_bases+'\\%', str(round(d[args.sample[i]]['mean_coverage'], DECIMALS)),
                          str(round(d[args.sample[i]]['mean20'], DECIMALS)), str(round(d[args.sample[i]]['cov_unif'], DECIMALS))+'\\%',
                          str(len(d[args.sample[i]]['cov_list'])), str(d[args.sample[i]]['covered_amplicons'])])

        qc_info = list_to_tex_table(lines, header=qc_header, long_columns=[1, 2, 3, 4, 5, 6, 7, 8], caption='QC information')
        template = template.replace('{#QCTABLE#}', qc_info)

        #Make Amplicons with coverage < 100% table

        cols = ['Sample', 'Amplicon', '\% Covered', 'Less than 20\% of mean']
        not_covered = []
        for i in indexlist:
            amp_cov = parse_mean_covd(args.meancov[i])
            not_covered = not_covered+[(_escape_latex(args.sample[i]).strip('.bam'), _escape_latex(line[4]),
                                        '{:.1f}'.format(float(line[8]) * 100), produce_warning(float(amp_cov[line[4]]), d[args.sample[i]]['mean20']))
                                       for line in d[args.sample[i]]['cov_list'] if float(line[8]) < 1]
        if not_covered == []:
            not_covered = [['/', '/']]

        table_text = list_to_tex_table(not_covered, header=cols, caption='Amplicons with coverage < 100\\%')
        template = template.replace('{#BAD_AMPLICON_TABLE#}', table_text)

        # Parse VCF files to get variant tables and dictionaries in the form {variant: samples}
        table_text = ''
        gatkhc_variants = {}
        lf_variants = {}

        header = ['CHROM', 'POS', 'REF', 'ALT', 'AF', 'DP', 'DP4', 'GEN[0].AD', 'SB', 'FS', 'EFF[*].GENE', 'ID',
                  'EFF[*].AA']
        header_glossary = {'GEN[0].AD': 'AD', 'EFF[*].GENE': 'GENE', 'EFF[*].AA': 'AA'}

        for i in indexlist:
            # GATKHC:
            vcf_table_1, common_columns_1 = _tsv_to_list(args.vcfgatkhc[i], has_header=True, pick_columns=header)

            # Escape user inputs and change header:
            common_columns_1 = [header_glossary[x] if (x in header_glossary.keys()) else x for x in common_columns_1]
            caption = _escape_latex('GATK HaplotypeCaller variant calls, sample ' + args.sample[i].strip('.bam'))
            vcf_table_1 = [[_escape_latex(value) for value in line] for line in vcf_table_1]

            # Insert space between SNP ID's and create hypelinks:
            vcf_table_1 = [
                line[:-2] +
                [' '.join(map(snp_href, line[-2].split(';')))] +
                [aa_change(line[-1].split(','))]
                for line in vcf_table_1
            ]

            # Create dict of shared variants
            for line in vcf_table_1:
                variant = line[-3] + '_chr' + line[0] + '_' + line[1]
                if variant in gatkhc_variants.keys():
                    gatkhc_variants[variant].append([args.sample[i].strip('.bam'), line[4]])
                else:
                    gatkhc_variants[variant] = [[args.sample[i].strip('.bam'), line[4]]]

            # Create gene hypelinks:
            vcf_table_1 = [line[:-3] + [gene_href(line[-3])] + line[-2:] for line in vcf_table_1]

            # Set different table counting (Na):
            table_text += '\\renewcommand{\\thetable}{\\arabic{table}a}'
            # Add table text:
            table_text += list_to_tex_table(vcf_table_1, header=common_columns_1, caption=caption,
                                            long_columns=[2, 3, -2, -1], uncut_columns=[-1])
            table_text += '{\n\\addtocounter{table}{-1}}'
            table_text += '\n\\newpage\n'

            # LF:
            vcf_table_2, common_columns_2 = _tsv_to_list(args.vcflf[i], has_header=True, pick_columns=header)

            # Escape user inputs and change header:
            common_columns_2 = [header_glossary[x] if (x in header_glossary.keys()) else x for x in common_columns_2]
            caption = _escape_latex('Lowfreq variant calls, sample ' + args.sample[i].strip('.bam'))
            vcf_table_2 = [[_escape_latex(value) for value in line] for line in vcf_table_2]

            # Insert space between SNP ID's and create hypelinks:
            vcf_table_2 = [
                line[:-2] +
                [' '.join(map(snp_href, line[-2].split(';')))] +
                [aa_change(line[-1].split(','))]
                for line in vcf_table_2
            ]

            # Create dict of shared variants
            for line in vcf_table_2:
                variant = line[-3] + '_chr' + line[0] + '_' + line[1]
                if variant in lf_variants.keys():
                    lf_variants[variant].append([args.sample[i].strip('.bam'), line[4]])
                else:
                    lf_variants[variant] = [[args.sample[i].strip('.bam'), line[4]]]

            # Create gene hypelinks:
            vcf_table_2 = [line[:-3] + [gene_href(line[-3])] + line[-2:] for line in vcf_table_2]

            # Set different table counting (Nb)
            table_text += '\\renewcommand{\\thetable}{\\arabic{table}b}'
            table_text += list_to_tex_table(vcf_table_2, header=common_columns_2, caption=caption,
                                            long_columns=[2, 3, -2, -1], uncut_columns=[-1])
            table_text += '\n\\newpage\n'
        # Set table counter back to normal (N) for further tables
        table_text += '\\renewcommand{\\thetable}{\\arabic{table}}'

        # Create shared variants tables
        header_shared = ['Variant ID', 'Samples sharing this variant (allele frequence)']
        data_gatkhc = []

        # Sort data by number of samples sharing the variant
        gatkhc_sorted = sorted(gatkhc_variants, key=lambda k: (-len(gatkhc_variants[k]),k))
        for k in gatkhc_sorted:
            data_gatkhc.append([_escape_latex(k), _escape_latex(samplelist_to_string(gatkhc_variants[k]))])
            caption = 'Shared GATK HaplotypeCaller variants'
        # Set table counting and add table text
        gatkhc_shared = '\\renewcommand{\\thetable}{\\arabic{table}a}\n'
        gatkhc_shared += list_to_tex_table(data_gatkhc, header=header_shared, caption=caption, wide_columns=[1])
        gatkhc_shared += '{\n\\addtocounter{table}{-1}}'

        data_lf = []

        # ort data by number of samples sharing the variant
        lf_sorted = sorted(lf_variants, key=lambda k: (-len(lf_variants[k]), k))
        for k in lf_sorted:
            data_lf.append([_escape_latex(k), _escape_latex(samplelist_to_string(lf_variants[k]))])
            caption = 'Shared Lowfreq variants'
        # Set table counting and add table text
        lf_shared = '\n\\renewcommand{\\thetable}{\\arabic{table}b}'
        lf_shared += list_to_tex_table(data_lf, header=header_shared, caption=caption, wide_columns=[1])
        gatkhc_shared += '\n\\renewcommand{\\thetable}{\\arabic{table}}'

        template = template.replace('{#GATKHC_SHARED#}', gatkhc_shared)
        template = template.replace('{#LF_SHARED#}', lf_shared)

        template = template.replace('{#VCF_TABLES#}', table_text)

        # Crate heatmaps
        make_heatmap(args.sample, gatkhc_variants, 'GATKHC variants')
        make_heatmap(args.sample, lf_variants, 'Lowfreq variants')

        # Write template to 'report.tex'
        template_out.write(template)

    # Generate PDF. Call teh command two times since the first time
    args = ['pdflatex', '-interaction=nonstopmode', 'multireport.tex']
    subprocess.call(args)
    subprocess.check_call(args)

    print("Done.")
