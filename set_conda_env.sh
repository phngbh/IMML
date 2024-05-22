#!/usr/bin/env bash

# Path to the YAML file containing the conda environment name
YAML_FILE_PATH="./workflows/FS_config.yml"

# Extract the conda_env value from the YAML file
CONDA_ENV=$(yq e '.conda_env' $YAML_FILE_PATH)

# Path to the target file to modify
TARGET_FILE_PATH="./workflows/feature_selection.nf"

# Replace any characters after 'conda' on lines containing 'conda' with the conda_env value
sed -i "/conda/c\conda '$CONDA_ENV'" "$TARGET_FILE_PATH"

