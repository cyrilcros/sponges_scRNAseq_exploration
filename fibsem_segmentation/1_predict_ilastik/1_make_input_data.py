#! /g/kreshuk/pape/Work/software/conda/miniconda3/envs/cluster_env37/bin/python

import sys
import os
import json
import luigi
from cluster_tools.ilastik import StackPredictionsSlurm


def make_stacked_input_n5(pred_key, out_key):
    path = '/g/kreshuk/data/arendt/sponge/data.n5'
    raw_key = 'volumes/raw/s0'
    dtype = 'uint8'

    max_jobs = 200
    tmp_folder = './tmp_stack'
    config_folder = './config_stack'
    os.makedirs(config_folder, exist_ok=True)

    shebang = "#! /g/kreshuk/pape/Work/software/conda/miniconda3/envs/cluster_env37/bin/python"
    global_config = StackPredictionsSlurm.default_global_config()
    block_shape = [128] * 3

    global_config.update({'shebang': shebang, 'block_shape': block_shape})
    with open('./config_stack/global.config', 'w') as f:
        json.dump(global_config, f)

    task = StackPredictionsSlurm(tmp_folder=tmp_folder, max_jobs=max_jobs,
                                 config_dir=config_folder,
                                 raw_path=path, raw_key=raw_key,
                                 pred_path=path, pred_key=pred_key,
                                 output_path=path, output_key=out_key, dtype=dtype)
    ret = luigi.build([task], local_scheduler=True)
    assert ret


def make_input(stage):
    if stage == 1:
        key = 'volumes/predictions/semantic_stage0'
        out_key = 'input_stage1'
    elif stage == 2:
        key = 'volumes/predictions/semantic_stage1'
        out_key = 'input_stage2'
    make_stacked_input_n5(key, out_key)


if __name__ == '__main__':
    stage = int(sys.argv[1])
    if stage not in (1, 2):
        raise ValueError("Expected stage in [1, 2], got %i" % stage)
    make_input(stage)
