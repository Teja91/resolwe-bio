#!/usr/bin/env python2
"""Generate amplicon multireport."""

import subprocess
import argparse

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

if __name__ == '__main__':
    args = parser.parse_args()
    
    #Make a index list of sorted sample names (alphabetical order)
    indexlist=[i[0] for i in sorted(enumerate(args.sample), key=lambda x:x[1])]

    # Open template and fill it with data:
    with open(args.template, 'r') as template_in, open('multireport.tex', 'wt') as template_out:
        template = template_in.read()
        template = template.replace('{#LOGO#}', args.logo)

        #make QC information table
        lines=[]
        lines.append('\\begin{longtable}{l l l l l l l l}')
        lines.append('\\rowcolor{darkblue1}')
        lines.append('\\leavevmode\\color{white}\\textbf{' +
                     '}& \\leavevmode\\color{white}\\textbf{'.join(['Sample name', 'Total reads', 'Aligned reads', 'Aligned bases', 'Mean coverage', 'Threshold coverage', 'Coverage uniformity', 'No. of amplicons']) + '} \\\\')
        for i in indexlist:
            lines.append(' & '.join([args.sample[i], '0','0','0','0','0','0','0']) + ' \\\\')

        lines.append('\\end{longtable}')
        QCinfo='\n'.join(lines)

        template=template.replace('{#QCTABLE#}', QCinfo)

        # Write template to 'report.tex'
        template_out.write(template)

    # Generate PDF:
    args = ['pdflatex', '-interaction=nonstopmode', 'multireport.tex']
    subprocess.call(args)
    subprocess.check_call(args)

    print("Done.")
