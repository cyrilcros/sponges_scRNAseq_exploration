#! /g/kreshuk/pape/Work/software/conda/miniconda3/envs/cluster_env37/bin/python

import sys
import os
import json
import luigi

from cluster_tools.ilastik import IlastikPredictionWorkflow


def run_pred(ilastik_project, input_key, output_key, n_channels,
             tmp_folder, max_jobs, max_threads, target):

    input_path = '/g/kreshuk/data/arendt/sponge/data.n5'
    output_path = '/g/kreshuk/data/arendt/sponge/data.n5'

    configs = IlastikPredictionWorkflow.get_config()

    config_folder = 'config'
    if not os.path.exists(config_folder):
        os.mkdir(config_folder)

    shebang = "#! /g/kreshuk/pape/Work/software/conda/miniconda3/envs/cluster_env37/bin/python"
    global_config = configs['global']

    # full run
    block_shape = [256, 256, 256]
    roi_begin = roi_end = None

    # debug
    # block_shape = [128, 128, 128]
    # roi_begin = [1592, 1988, 1789]
    # roi_end = [1692, 3012, 2813]

    global_config.update({'shebang': shebang,
                          'block_shape': block_shape,
                          'roi_begin': roi_begin,
                          'roi_end': roi_end})
    with open('./config/global.config', 'w') as f:
        json.dump(global_config, f)

    conf = configs['prediction']
    conf.update({'time_limit': 1440, 'mem_limit': 16, 'threads_per_job': max_threads,
                 'qos': 'normal'})
    with open('./config/prediction.config', 'w') as f:
        json.dump(conf, f)

    ilastik_folder = '/g/kreshuk/pape/Work/software/src/ilastik-meta/ilastik'
    halo = [32, 32, 32]

    ret = luigi.build([IlastikPredictionWorkflow(input_path=input_path, input_key=input_key,
                                                 output_path=output_path, output_key=output_key,
                                                 ilastik_project=ilastik_project,
                                                 ilastik_folder=ilastik_folder,
                                                 halo=halo, n_channels=n_channels,
                                                 config_dir=config_folder,
                                                 tmp_folder=tmp_folder,
                                                 target=target,
                                                 max_jobs=max_jobs)], local_scheduler=True)
    assert ret


def predict_stage0():
    tmp_folder = './tmp'
    target = 'slurm'
    max_jobs = 200
    max_threads = 4
    ilastik_project = '/g/kreshuk/pape/Work/my_projects/spongilla_experiments/ilastik/semantic_pixclass_stage0.ilp'
    in_key = 'volumes/raw/s0'
    out_key = 'volumes/predictions/semantic_stage0'
    n_channels = 6
    run_pred(ilastik_project, in_key, out_key, n_channels,
             tmp_folder, max_jobs, max_threads, target)


def predict_stage1():
    tmp_folder = './tmp_stage1'
    target = 'slurm'
    max_jobs = 225
    max_threads = 4
    ilastik_project = '/g/kreshuk/pape/Work/my_projects/spongilla_experiments/ilastik/semantic_pixclass_stage1.ilp'
    in_key = 'input_stage1'
    out_key = 'volumes/predictions/semantic_stage1'
    n_channels = 6
    run_pred(ilastik_project, in_key, out_key, n_channels,
             tmp_folder, max_jobs, max_threads, target)


def predict_stage2():
    tmp_folder = './tmp_stage2'
    target = 'slurm'
    max_jobs = 225
    max_threads = 4
    ilastik_project = '/g/kreshuk/pape/Work/my_projects/spongilla_experiments/ilastik/binary_pixclass_stage2.ilp'
    in_key = 'input_stage2'
    out_key = 'volumes/predictions/binary_stage2'
    n_channels = 2
    run_pred(ilastik_project, in_key, out_key, n_channels,
             tmp_folder, max_jobs, max_threads, target)


if __name__ == '__main__':
    stage = int(sys.argv[1])
    if stage == 0:
        predict_stage0()
    elif stage == 1:
        predict_stage1()
    elif stage == 2:
        predict_stage2()
    else:
        raise ValueError("Expected stage in [0, 1, 2], got %i" % stage)
