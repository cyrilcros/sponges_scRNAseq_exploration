#! /g/kreshuk/pape/Work/software/conda/miniconda3/envs/cluster_env37/bin/python

import os
import json
import luigi

from cluster_tools import ThresholdedComponentsWorkflow

CLASS_TO_CHANNEL = {'cells': [2, 3], 'flagella': 4, 'microvillae': 5}


def run_thresh(semantic_class, max_jobs, max_threads, tmp_folder, target):
    input_path = '/g/kreshuk/data/arendt/sponge/data.n5'
    input_key = 'volumes/predictions/semantic_stage1'
    output_key = 'volumes/predictions/classes/%s' % semantic_class
    channel = CLASS_TO_CHANNEL[semantic_class]

    configs = ThresholdedComponentsWorkflow.get_config()

    config_folder = './config_thresh'
    if not os.path.exists(config_folder):
        os.mkdir(config_folder)

    shebang = "#! /g/kreshuk/pape/Work/software/conda/miniconda3/envs/cluster_env37/bin/python"
    global_config = configs['global']
    global_config.update({'shebang': shebang,
                          'block_shape': [256, 256, 256]})
    with open(os.path.join(config_folder, 'global.config'), 'w') as f:
        json.dump(global_config, f)

    task = ThresholdedComponentsWorkflow(tmp_folder=tmp_folder, config_dir=config_folder,
                                         max_jobs=max_jobs, target=target,
                                         input_path=input_path, input_key=input_key,
                                         output_path=input_path, output_key=output_key,
                                         assignment_key='components/%s' % semantic_class,
                                         threshold=.5, threshold_mode='greater',
                                         channel=channel)
    ret = luigi.build([task], local_scheduler=True)
    assert ret, "Failure"


def run_workflow(semantic_class, target):
    if target == 'slurm':
        max_jobs = 200
        max_threads = 8
    else:
        max_jobs = 8
        max_threads = 1

    tmp_folder = './tmp_thresh_%s' % semantic_class
    run_thresh(semantic_class, max_jobs, max_threads, tmp_folder, target=target)


if __name__ == '__main__':
    for class_ in ('flagella', 'microvillae'):
        run_workflow(class_, 'slurm')
