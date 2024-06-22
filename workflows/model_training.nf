
process model_training {
    input:
    val iter
    val integration
    val algorithm
    val p_metric
    val feature
    val from_clinical
    val fs_config_file
    val package_path
    output:
    val 0 

    script:
    """
    Rscript ${package_path}R/run_model_training.R --iter ${iter} --integration ${integration} --algorithm ${algorithm} --p_metric ${p_metric} --feature ${feature} --from_clinical ${from_clinical} --config ${fs_config_file}
    """
}

process report {
    input:
    val n_iterations
    val integration
    val algorithm
    val p_metric
    val feature
    val from_clinical
    val fs_config_file
    val package_path
    val annotation_dir
    val target_name
    val out_dir
    val model_train_exit

    script:
    if (from_clinical == 'True')
        """
        Rscript -e "rmarkdown::render('${package_path}report.Rmd', output_dir='${out_dir}', params=list(args=list(fs_config_file='${fs_config_file}', annotation_folder='${annotation_dir}', n_iterations=${n_iterations}, integration='${integration}', algorithm='${algorithm}', p_metric='${p_metric}', feature='${feature}', from_clinical=as.logical('${from_clinical}'))), output_file='${target_name}_${integration}_${algorithm}_${p_metric}_${feature}_fromClinical_model_report.html')"
        """
    else
        """
        Rscript -e "rmarkdown::render('${package_path}report.Rmd', output_dir='${out_dir}', params=list(args=list(fs_config_file='${fs_config_file}', annotation_folder='${annotation_dir}', n_iterations=${n_iterations}, integration='${integration}', algorithm='${algorithm}', p_metric='${p_metric}', feature='${feature}', from_clinical=as.logical('${from_clinical}'))), output_file='${target_name}_${integration}_${algorithm}_${p_metric}_${feature}_model_report.html')"
        """
}

workflow {
    iterations = Channel.from( 1.. params.n_iterations )
    model_train_exit = model_training(iterations, params.integration, params.algorithm, params.p_metric, params.feature, params.from_clinical, params.fs_config_file, params.package_path).collect()
    
    if( params.integration == 'FFS' ) {
        fs_config = new org.yaml.snakeyaml.Yaml().load(file(params.fs_config_file))
        report(params.n_iterations, params.integration, params.algorithm, params.p_metric, params.feature, params.from_clinical, params.fs_config_file, params.package_path, params.annotation_dir, fs_config.target_name, fs_config.out_dir, model_train_exit)
    }
    }