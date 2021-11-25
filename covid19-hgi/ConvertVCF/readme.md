### ConvertVCF

1. Check the parameters in the .json file (should change only `convert_bgen.chrom_list` and `convert_bgen.name`, which is the prefix for the generated bgen file (e.g. cov_ita_1.0.bgen))  
1.b. `convert_bgen.chrom_convert.chunk` is the size of each bgen chunk.
2. Connect to Cromwell server using cromwell interact:  
`./CromwellInteract-master/cromwell_interact.py connect cromwell-covid`
3. Submit workflow (use a sensible label so you can retrieve the job id in ./CromwellInteract-master/workflows.log):  
`python3 ./CromwellInteract-master/cromwell_interact.py submit --wdl Projects/covid19-hgi/covid19-hgi-repo/ConvertVCF/wdl/bgen_convert.wdl --inputs Projects/covid19-hgi/covid19-hgi-repo/ConvertVCF/wdl/bgen_convert.json --label convert_ita`
