# /pubhome/xli02/opt/miniconda/envs/reinvent.v3.2/lib/python3.7/site-packages/reinvent_scoring/scoring/enums/component_specific_parameters_enum.py

from dataclasses import dataclass


@dataclass(frozen=True)
class ComponentSpecificParametersEnum:
    SCIKIT = "scikit"
    CLAB_INPUT_FILE = "clab_input_file"
    DESCRIPTOR_TYPE = "descriptor_type"
    TRANSFORMATION = "transformation"

    # structural components
    # ---------
    # AZDOCK
    AZDOCK_CONFPATH = "configuration_path"
    AZDOCK_DOCKERSCRIPTPATH = "docker_script_path"
    AZDOCK_ENVPATH = "environment_path"
    AZDOCK_DEBUG = "debug"

    # DockStream
    DOCKSTREAM_CONFPATH = "configuration_path"
    DOCKSTREAM_DOCKERSCRIPTPATH = "docker_script_path"
    DOCKSTREAM_ENVPATH = "environment_path"
    DOCKSTREAM_DEBUG = "debug"

    # ICOLOS
    ICOLOS_CONFPATH = "configuration_path"
    ICOLOS_EXECUTOR_PATH = "executor_path"
    ICOLOS_VALUES_KEY = "values_key"
    ICOLOS_DEBUG = "debug"
    #######################

    RAT_PK_PROPERTY = "rat_pk_property"
    CLAB_TOP_20_VALUE = "clab_top_20_value"
    ION_CLASS = "Ion class"
    CONTAINER_TYPE = "container_type"

    SMILES = "smiles"
    MODEL_PATH = "model_path"

    #######################
    ARTIFACT = "artifact"

    AIZYNTH_CONFIG_FILE_PATH = "aizynth_config_file_path"

    VALUE_MAPPING = "value_mapping"

    # add for shafts 3D similarity by xli 20230626
    SIMI_CONFIG_PATH = "simi_config_path"
