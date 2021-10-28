import "meta.sub.wdl" as sub

workflow meta_analysis {

    File phenos
    File summary_files
    String conf_template

    Array[Array[String]] pheno = read_tsv(phenos)
    Array[Array[String]] summary_stats = read_tsv(summary_files)

    scatter (i in range(length(pheno))) {
        String conf = sub(conf_template, "\\{PHENO\\}", pheno[i][0])
        call sub.run_meta {
            input: pheno = pheno[i][0], min_n_eff = pheno[i][1], min_n_eff_x = pheno[i][2], conf = conf, summary_stats = summary_stats[i]
        }
    }
}
