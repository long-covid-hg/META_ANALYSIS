version 1.0

workflow meta_analysis {

    input {
        File sumstats_loc
        File pheno_confs

        Array[Array[String]] sumstat_files = read_tsv(sumstats_loc)
        Array[String] pheno_conf = read_lines(pheno_confs)
    }

    scatter (i in range(length(pheno_conf))) {

        String pheno = basename(pheno_conf[i], ".json")
        
        scatter (chr in range(23)) {
            call run_range {
                input:
                    pheno = pheno,
                    conf = pheno_conf[i],
                    sumstat_files = sumstat_files[i],
                    chrom = chr+1
            }
        }

        call combine_chrom_metas {
            input:
                pheno = pheno,
                meta_outs = run_range.out
        }

        call add_rsids {
            input:
                meta_file = combine_chrom_metas.meta_out
        }

        call meta_qq {
            input:
                meta_file = combine_chrom_metas.meta_out
        }

        call post_filter {
            input:
                meta_file = add_rsids.meta_out
        }

    }

    output {
        Array[File] metas = combine_chrom_metas.meta_out
        Array[File] metas_with_rsids = add_rsids.meta_out
        Array[File] filtered_metas = post_filter.filtered_meta_out
        Array[Array[File]] meta_pngs = meta_qq.pngs
        Array[Array[File]] meta_lambdas = meta_qq.lambdas
    }
}

# Run meta-analysis for each chromosome separately
task run_range {

    input {
        Array[File] sumstat_files
        String pheno
        String chrom
        File conf

        File script

        String docker
        String method
        String opts
    }

    command <<<

        echo "`date` GWAS meta-analysis"
        echo "docker: ~{docker}"
        echo "pheno: ~{pheno}"
        echo "method: ~{method}"
        echo "options: ~{opts}"
        echo "chromosome: ~{chrom}"
        echo "conf: ~{conf}"

        python3 ~{script} ~{conf} ~{pheno}_chr~{chrom}_meta_out.tsv ~{method} ~{opts} --chrom ~{chrom}

        echo "`date` done"
    >>>

    output {
        File out = pheno + "_chr" + chrom + "_meta_out.tsv.gz"
    }

    runtime {
        docker: "~{docker}"
        cpu: 1
        memory: "10 GB"
        disks: "local-disk " + 4 * ceil(size(sumstat_files, "GB") + 1) + " HDD"
        zones: "us-central1-b"
        preemptible: 2
        noAddress: true
    }
}

# Combine separately run meta-analysis result files
task combine_chrom_metas {

    input {
        Array[File] meta_outs
        String pheno

        String docker
    }

    command <<<

        echo "`date` Combining metadata files"
        echo "docker: ~{docker}"
        echo "pheno: ~{pheno}"

        cat <(zcat ~{meta_outs[0]} | head -1) \
            <(for file in ~{sep=" " meta_outs}; do
                zcat $file | tail -n +2;
            done) \
        | bgzip > ~{pheno}_meta_out.tsv.gz

        echo "`date` tabixing"
        tabix -s 1 -b 2 -e 2 ~{pheno}_meta_out.tsv.gz
        echo "`date` done"
    >>>

    output {
        File meta_out = pheno + "_meta_out.tsv.gz"
        File meta_out_tbi = pheno + "_meta_out.tsv.gz.tbi"
    }
    
    runtime {
        docker: "~{docker}"
        cpu: 1
        memory: "10 GB"
        disks: "local-disk " + 4 * ceil(size(meta_outs, "GB") + 1) + " HDD"
        zones: "us-central1-b"
        preemptible: 2
        noAddress: true
    }
}

# Add rsids
task add_rsids {

    input {
        File meta_file

        File ref_file
        String docker

        String base = basename(meta_file, ".tsv.gz")
    }

    command <<<

        set -euxo pipefail

        echo "`date` Adding rsids"

        python3 <<EOF | bgzip > ~{base}.tsv.gz

        import gzip

        fp_ref = gzip.open('~{ref_file}', 'rt')
        ref_has_lines = True
        ref_chr = 1
        ref_pos = 0
        ref_line = fp_ref.readline()
        while ref_line.startswith("##"):
            ref_line = fp_ref.readline()
        if ref_line.startswith('#'):
            assert ref_line.rstrip('\r\n').split('\t') == '#CHROM POS ID REF ALT QUAL FILTER INFO'.split(), repr(ref_line)
        ref_h_idx = {h:i for i,h in enumerate(ref_line.rstrip('\r\n').split('\t'))}

        with gzip.open('~{meta_file}', 'rt') as f:
            header = f.readline().strip()
            h_idx = {h:i for i,h in enumerate(header.split('\t'))}
            print(header + '\trsid')
            for line in f:
                line = line.strip()
                s = line.split('\t')
                chr = int(s[h_idx['#CHR']])
                pos = int(s[h_idx['POS']])
                ref = s[h_idx['REF']]
                alt = s[h_idx['ALT']]
                ref_vars = []
                while ref_has_lines and int(ref_chr) < chr or (int(ref_chr) == chr and ref_pos < pos):
                    ref_line = fp_ref.readline().rstrip('\r\n').split('\t')
                    try:
                        ref_chr = ref_line[ref_h_idx['#CHROM']]
                        ref_pos = int(ref_line[ref_h_idx['POS']])
                    except ValueError:
                        ref_has_lines = False
                while ref_has_lines and int(ref_chr) == chr and ref_pos == pos:
                    ref_vars.append(ref_line)
                    ref_line = fp_ref.readline().strip().split('\t')
                    try:
                        ref_chr = ref_line[ref_h_idx['#CHROM']]
                        ref_pos = int(ref_line[ref_h_idx['POS']])
                    except ValueError:
                        ref_has_lines = False

                rsid = 'NA'
                for r in ref_vars:
                    if r[ref_h_idx['REF']] == ref and alt in r[ref_h_idx['ALT']].split(','):
                        rsid = r[ref_h_idx['ID']]
                        break

                print(line + '\t' + rsid)

        EOF

        echo "`date` tabixing"
        tabix -s 1 -b 2 -e 2 ~{base}.tsv.gz
        echo "`date` done"

    >>>

    output {
        File meta_out = base + ".tsv.gz"
        File out_tbi = base + ".tsv.gz.tbi"
    }

    runtime {
        docker: "~{docker}"
        cpu: "1"
        memory: "10 GB"
        disks: "local-disk " + 4*ceil(size(meta_file, "G") + size(ref_file, "G")) + " SSD"
        zones: "us-central1-b"
        preemptible: 2
        noAddress: true
    }

}

# Generate qq and manhattan plots from meta-analysis results
task meta_qq {

    input {
        File meta_file

        Int loglog_ylim
        String docker
        String pvals_to_plot

        File script

        String base = basename(meta_file)
    }

    command <<<

        set -euxo pipefail

        mv ~{meta_file} ~{base}

        Rscript ~{script} --file ~{base} --bp_col "POS" --chrcol "#CHR" --pval_col ~{pvals_to_plot} --loglog_ylim ~{loglog_ylim}

    >>>

    output {
        Array[File] pngs = glob("*.png")
        Array[File] lambdas = glob("*.txt")
    }

    runtime {
        docker: "~{docker}"
        cpu: "1"
        memory: "30 GB"
        disks: "local-disk " + 20*ceil(size(meta_file, "G")) + " HDD"
        zones: "us-central1-b"
        preemptible: 2
        noAddress: true
    }
}

# Filter out variants in <2 studies and extract variant details and meta-analysis results columns only
task post_filter {

    input {
        File meta_file

        String docker

        String base = basename(meta_file, ".tsv.gz")
    }

    command <<<

        set -exo pipefail

        # Use the "all_meta_Nstudies" column to remove variants in <2 studies
        zcat ~{meta_file} | awk -v OFS='\t' '
        NR==1 {str="#CHR\tPOS\tREF\tALT\tSNP\tRSID";for(i=1;i<=NF;i++){col[$i]=i;if($i~/^(all_|lmso_)/&&$i~/(_meta_|_het_)/){c[i]++;str=str"\t"$i}};print str}
        (NR>1 && $col["all_meta_Nstudies"] != 1) {str=$col["#CHR"]"\t"$col["POS"]"\t"$col["REF"]"\t"$col["ALT"]"\t"$col["SNP"]"\t"$col["rsid"];for(i=1;i<=NF;i++){if(i in c){str=str"\t"$i}};print str}
        ' | bgzip > ~{base}_filtered.tsv.gz
        tabix -s 1 -b 2 -e 2 ~{base}_filtered.tsv.gz

    >>>

    output {
        File filtered_meta_out = base + "_filtered.tsv.gz"
        File filtered_meta_out_tbi = base + "_filtered.tsv.gz.tbi"
    }

    runtime {
        docker: "~{docker}"
        cpu: "1"
        memory: "2 GB"
        disks: "local-disk " + 3*ceil(size(meta_file, "G")) + " HDD"
        zones: "us-central1-b"
        preemptible: 2
        noAddress: true
    }
}
