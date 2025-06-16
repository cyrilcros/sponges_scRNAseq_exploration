#! /g/kreshuk/pape/Work/software/conda/miniconda3/envs/cluster_env37/bin/python

import os
import json
import luigi
from cluster_tools.paintera import ConversionWorkflow


def paintera_format(max_jobs, max_threads, tmp_folder, target):
    # input_path = '/g/kreshuk/data/arendt/sponge/data.n5'
    input_path = '/g/emcf/Mizzon/projects/Klaske_FIBSEM/FIB-SEM/paintera/segmentation_v1.n5'

    configs = ConversionWorkflow.get_config()

    config_folder = './config_mc'
    if not os.path.exists(config_folder):
        os.mkdir(config_folder)

    shebang = "#! /g/kreshuk/pape/Work/software/conda/miniconda3/envs/cluster_env37/bin/python"
    global_config = configs['global']
    global_config.update({'shebang': shebang, 'block_shape': [128] * 3})

    with open(os.path.join(config_folder, 'global.config'), 'w') as f:
        json.dump(global_config, f)

    config_names = ["downscaling"]
    for name in config_names:
        config = configs[name]
        config.update({"library_kwargs": {"order": 0}, "mem_limit": 8,
                       "time_limit": 180})
        with open(os.path.join(config_folder, "%s.config" % name), 'w') as f:
            json.dump(config, f)

    block_mapping_conf = configs['label_block_mapping']
    block_mapping_conf.update({'mem_limit': 100, 'time_limit': 360,
                               'threads_per_job': max_threads})
    with open(os.path.join(config_folder, 'label_block_mapping.config'), 'w') as f:
        json.dump(block_mapping_conf, f)

    assignment_path = assignment_key = ''
    label_in_key = 'volumes/segmentation/lifted_multicut'
    label_out_key = 'volumes/segmentation'

    task = ConversionWorkflow(tmp_folder=tmp_folder, config_dir=config_folder,
                              max_jobs=max_jobs, path=input_path,
                              target=target,
                              label_in_key=label_in_key,
                              label_out_key=label_out_key,
                              label_scale=0,
                              raw_key='volumes/raw',
                              resolution=[15, 15, 15],
                              assignment_path=assignment_path,
                              assignment_key=assignment_key)
    ret = luigi.build([task], local_scheduler=True)
    assert ret


if __name__ == '__main__':
    target = 'slurm'
    if target == 'slurm':
        max_jobs = 400
        max_threads = 8
        max_jobs_mc = 8
    else:
        max_jobs = 8
        max_threads = 1
        max_jobs_mc = 1

    tmp_folder = './tmp_export'
    paintera_format(max_jobs, max_threads, tmp_folder, target=target)
