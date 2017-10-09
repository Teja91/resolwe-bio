# pylint: disable=missing-docstring
from os.path import join

from resolwe_bio.utils.test import skipUnlessLargeFiles, BioProcessTestCase


class ReportProcessorTestCase(BioProcessTestCase):

    def prepare_report_inputs(self, bam_file, mfile, pname, target_pcr, target_cov, annotations, summaries, snpeffs):
        """Prepare report inputs."""
        bam = self.run_process('upload-bam', {'src': join('large', bam_file)})
        master_file = self.prepare_amplicon_master_file(mfile=mfile, pname=pname)
        coverage = self.run_process('coveragebed', {'alignment': bam.id, 'master_file': master_file.id})

        inputs = {'target_pcr_metrics': target_pcr, 'target_coverage': target_cov}
        pcr_metrics = self.run_process('upload-picard-pcrmetrics', inputs)

        annot_vars = []
        assert len(annotations) == len(summaries) == len(snpeffs)
        for ann, summ, snpeff in zip(annotations, summaries, snpeffs):
            inputs = {'annotation': ann, 'summary': summ, 'snpeff_genes': snpeff}
            annot_var = self.run_process('upload-snpeff', inputs)
            annot_vars.append(annot_var.id)

        return coverage, pcr_metrics, master_file, annot_vars

    @skipUnlessLargeFiles('56GSID_10k_mate1_RG.bam')
    def test_amplicon_report(self):
        template = self.run_process('upload-file', {'src': 'report_template.tex'})
        logo = self.run_process('upload-file', {'src': 'genialis_logo.pdf'})

        coverage, pcr_metrics, master_file, annot_vars = self.prepare_report_inputs(
            bam_file='56GSID_10k_mate1_RG.bam',
            mfile='56G_masterfile_test.txt',
            pname='56G panel, v2',
            target_pcr='56gsid_10k.targetPCRmetrics.txt',
            target_cov='56gsid_10k.perTargetCov.txt',
            annotations=['56GSID.lf.finalvars.txt'],
            summaries=['56GSID_1k.gatkHC_snpEff_summary.html'],
            snpeffs=['56GSID_1k.gatkHC_snpEff_genes.txt'],
        )

        report_inputs = {
            'coverage': coverage.id,
            'pcr_metrics': pcr_metrics.id,
            'template': template.id,
            'logo': logo.id,
            'master_file': master_file.id,
            'annot_vars': annot_vars
        }
        self.run_process('amplicon-report', report_inputs)

    def test_multisample_amplicon_report(self):
        self.keep_data()
        template = self.run_process('upload-file', {'src': 'multireport-template.tex'})
        logo = self.run_process('upload-file', {'src': 'genialis_logo.pdf'})

        coverage1, pcr_metrics1, master_file1, annot_vars1 = self.prepare_report_inputs(
            bam_file='56GSID_10k_mate1_RG.bam',
            mfile='56G_masterfile_test.txt',
            pname='56G panel, v2',
            target_pcr='56gsid_10k.targetPCRmetrics.txt',
            target_cov='56gsid_10k.perTargetCov.txt',
            annotations=['56GSID.lf.finalvars.txt', '56GSID.gatkHC.finalvars.txt'],
            summaries=['56GSID_1k.gatkHC_snpEff_summary.html', '56GSID_1k.gatkHC_snpEff_summary.html'],
            snpeffs=['56GSID_1k.gatkHC_snpEff_genes.txt', '56GSID_1k.gatkHC_snpEff_genes.txt'],
        )
        coverage2, pcr_metrics2, master_file2, annot_vars2 = self.prepare_report_inputs(
            bam_file='56GSID_10k_mate1_RG.bam',
            mfile='56G_masterfile_test.txt',
            pname='56G panel, v2',
            target_pcr='56gsid_10k.targetPCRmetrics.txt',
            target_cov='56gsid_10k.perTargetCov.txt',
            annotations=['56GSID.lf.finalvars.txt', '56GSID.gatkHC.finalvars.txt'],
            summaries=['56GSID_1k.gatkHC_snpEff_summary.html', '56GSID_1k.gatkHC_snpEff_summary.html'],
            snpeffs=['56GSID_1k.gatkHC_snpEff_genes.txt', '56GSID_1k.gatkHC_snpEff_genes.txt'],
        )

        report_inputs = {
            'data': [coverage1.id, pcr_metrics1.id, master_file1.id] + annot_vars1 +
                    [coverage2.id, pcr_metrics2.id, master_file2.id] + annot_vars2,
            'fields': ['amplicon_cov'],
            'template': template.id,
            'logo': logo.id,
        }
        self.run_process('amplicon-archive-multi-report', report_inputs)
