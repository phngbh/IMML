
process data_partition {
    conda '/home/willy/miniconda3/envs/IMML/'
    input:
    tuple val(target), val(modals_ids), val(target_name), val(out_dir), val(seed), val(k), val(n_iter), val(p), val(package_path)
    output:
    val 0
    
    script:
    """
    Rscript ${package_path}R/data_part.R --target ${target} --modals_ids ${modals_ids} --target_name ${target_name} --out_dir ${out_dir} --seed ${seed} --k ${k} --n_iterations ${n_iter} --percent ${p}
    """
}

process feature_selection {
    conda '/home/willy/miniconda3/envs/IMML/'
    input:
    val config
    val target
    val modals_ids
    val target_name
    val out_dir
    val seed
    val package_path
    val data_part_exit

    script:
    if(config.modality == 'Clinical')
        """
        Rscript ${package_path}R/FS_clinical.R --data ${config.data} --target ${target} --modals_ids ${modals_ids} --target_name ${target_name} --out_dir ${out_dir} --seed ${seed} --n_iterations ${config.n_iterations} --p ${config.train_split} --resampling ${config.resampling} --partitioning ${out_dir}data_partition_${target_name}.rds
        """
    else if(config.modality == 'Genomics')
        """
        Rscript ${package_path}R/FS_genomics.R --target ${target} --modals_ids ${modals_ids} --target_name ${target_name} --out_dir ${out_dir} --seed ${seed} --n_iterations ${config.n_iterations} --p ${config.train_split} --partitioning ${out_dir}data_partition_${target_name}.rds --file_prefix ${config.file_prefix} --geneloc_file ${config.geneloc_file} --geneset_file ${config.geneset_file} 
        """
    else if(config.modality == 'Methylomics')
        """
        Rscript ${package_path}R/FS_methylomics.R --data ${config.data} --target ${target} --modals_ids ${modals_ids} --target_name ${target_name} --out_dir ${out_dir} --seed ${seed} --n_iterations ${config.n_iterations} --p ${config.train_split} --partitioning ${out_dir}data_partition_${target_name}.rds --pathways ${config.pathways} --dualmap_file ${config.dualmap_file} --DM_threshold ${config.DM_threshold}
        """
    else if(config.modality == 'targeted')
        """
        Rscript ${package_path}R/FS_targeted.R --data ${config.data} --target ${target} --modals_ids ${modals_ids} --target_name ${target_name} --out_dir ${out_dir} --seed ${seed} --n_iterations ${config.n_iterations} --p ${config.train_split} --resampling ${config.resampling} --partitioning ${out_dir}data_partition_${target_name}.rds --anno ${config.annotation} --msea_fdr ${config.msea_fdr} --mod_name ${config.name}
        """
    else if(config.modality == 'untargeted')
        """
        Rscript ${package_path}R/FS_untargeted.R --data ${config.data} --target ${target} --modals_ids ${modals_ids} --target_name ${target_name} --out_dir ${out_dir} --seed ${seed} --n_iterations ${config.n_iterations} --p ${config.train_split} --resampling ${config.resampling} --partitioning ${out_dir}data_partition_${target_name}.rds --anno ${config.annotation} --gsea_fdr ${config.gsea_fdr} --mod_name ${config.name}
        """
    else
        error "Invalid modality name: ${config.modality}"

}

workflow {
    data_part_config = tuple( params.targets, params.modals_ids, params.target_name, params.out_dir, params.seed, params.data_partitioning.k, params.data_partitioning.n_iterations, params.data_partitioning.percent, params.package_path )
    data_part_exit = data_partition(data_part_config)

    configs = Channel.from(params.feature_selection)
    feature_selection(configs, params.targets, params.modals_ids, params.target_name, params.out_dir, params.seed, params.package_path, data_part_exit)
}

