#!/usr/bin/env python
# coding: utf-8

import hail as hl
hl.init()

vcf_list = ['gs://dsge-covid19-data/test_SAIGE/chr%s.dose.vcf.gz' % chrom for chrom in range(19, 20)]

mt = hl.import_vcf(vcf_list, force_bgz=True, reference_genome='GRCh38')

t = hl.import_table('gs://dsge-covid19-data/test_SAIGE/chr19.gws.txt',
                    impute=True, delimiter=",", quote='"', no_header=True)

loci_to_extract = t.f0.collect()
mt = mt.filter_rows(hl.literal(loci_to_extract).contains(mt.rsid))

pheno = hl.import_table('gs://dsge-covid19-data/cases_controls_phenotypes_ita_be_brazil_swe_03_09_2020', impute=True)
pheno = pheno.annotate(s="0_"+pheno.ID)
pheno = pheno.key_by('s')

mt = mt.annotate_cols(pheno=pheno[mt.s])

ca = mt.filter_cols(mt.pheno.Analysis_C2 == 1)
co = mt.filter_cols(mt.pheno.Analysis_C2 == 0)

#print(ca.count())
#print(co.count())

ca = hl.variant_qc(ca)
ca_v = ca.rows()
ca_v.select(AC=ca_v.variant_qc.AC,
            AF_ref=ca_v.variant_qc.AF[0],
            AF_alt=ca_v.variant_qc.AF[1],
            n_het=ca_v.variant_qc.n_het,
            het_freq_hwe=ca_v.variant_qc.het_freq_hwe,
            p_value_hwe=ca_v.variant_qc.p_value_hwe).export("gs://dsge-covid19-data/test_SAIGE/chr19_be_cases.tsv")



co = hl.variant_qc(co)
co_v = co.rows()
co_v.select(AC=co_v.variant_qc.AC,
            AF_ref=co_v.variant_qc.AC[0],
            AF_alt=co_v.variant_qc.AF[1],
            n_het=co_v.variant_qc.n_het,
            het_freq_hwe=co_v.variant_qc.het_freq_hwe,
            p_value_hwe=co_v.variant_qc.p_value_hwe).export("gs://dsge-covid19-data/test_SAIGE/chr19_be_ctrls.tsv")



#mt_f = mt_f.key_rows_by('locus', 'alleles', 'rsid')
# mt_f.GT.n_alt_alleles().export('gs://dsgelab-sakari/doc_jukarainen_SNPs_R6_25082020.tsv')
