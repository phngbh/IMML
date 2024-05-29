
process model_training {
    input:
    val iter
    val integration
    val algorithm
    val p_metric
    val feature
    val fs_config_file
    val package_path

    script:
    """
    Rscript ${package_path}R/run_model_training.R --iter ${iter} --integration ${integration} --algorithm ${algorithm} --p_metric ${p_metric} --feature ${feature} --config ${fs_config_file}
    """
}

workflow {
    iterations = Channel.from( 1.. params.n_iterations )
    model_training(iterations, params.integration, params.algorithm, params.p_metric, params.feature, params.fs_config_file, params.package_path)
}