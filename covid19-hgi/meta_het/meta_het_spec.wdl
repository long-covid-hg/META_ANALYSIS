task run {

    File variant_shortlist
    File file1
    File file2
    String name1
    String name2
    String meta_name
    String docker

    command <<<

        python3 <<EOF

        import gzip, json

        with open('${variant_shortlist}', 'rt') as f:
            header = f.readline().strip().split('\t')
            variants = {line.strip().split('\t')[header.index('SNP')]: True for line in f}

        files = ['${file1}', '${file2}']
        names = ['${name1}', '${name2}']
        meta = {'meta': []}
        for i,file in enumerate(files):
            out = gzip.open(file + '.filtered.txt.gz', 'wt')
            with gzip.open(file, 'rt') as f:
                headerline = f.readline()
                out.write(headerline)
                header = headerline.strip().split('\t')
                snp_index = header.index('SNP')
                for line in f:
                    s = line.strip().split('\t')
                    if s[snp_index] in variants:
                        out.write(line)
            out.close()
            meta['meta'].append({
                'name': names[i],
                'file': file + '.filtered.txt.gz',
                'n_cases': 0,
                'n_controls': 0,
                'chr':'#CHR',
                'pos':'POS',
                'ref':'REF',
                'alt':'ALT',
                'effect':'all_inv_var_meta_beta',
                'effect_type':'beta',
                'pval':'all_inv_var_meta_p',
                'se':'all_inv_var_meta_sebeta',
                'extra_cols':['all_meta_sample_N', 'all_meta_AF', 'all_inv_var_het_p', 'rsid']
            })

        with open('meta.json', 'wt') as f:
            json.dump(meta, f, indent=4)
        EOF

        /META_ANALYSIS/scripts/meta_analysis.py --is_het_test meta.json "${meta_name}" inv_var

    >>>

    output {
        File out = meta_name + ".gz"
        File out_tbi = meta_name + ".gz.tbi"
    }

    runtime {
        docker: "${docker}"
        cpu: "1"
        memory: "2 GB"
        disks: "local-disk 200 HDD"
        zones: "us-east1-d"
        preemptible: 0
        noAddress: true
    }
}

workflow meta_het_spec {

    File input_file
    Array[Array[String]] inputs = read_tsv(input_file)

    scatter (input_set in inputs) {
        call run {
            input: file1 = input_set[0], file2 = input_set[1], name1 = input_set[2], name2 = input_set[3], variant_shortlist = input_set[4], meta_name = input_set[5]
        }
    }
}
