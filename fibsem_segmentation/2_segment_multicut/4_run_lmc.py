#! /g/kreshuk/pape/Work/software/conda/miniconda3/envs/cluster_env37/bin/python

import os
import json
import luigi
from cluster_tools import LiftedMulticutSegmentationWorkflow


def run_mc(max_jobs, max_threads, tmp_folder,
           target, max_jobs_mc=10):

    input_path = '/g/kreshuk/data/arendt/sponge/data.n5'
    input_key = 'volumes/predictions/boundaries'
    exp_path = '/g/kreshuk/data/arendt/sponge/multicut_data.n5'
    ws_key = 'volumes/segmentation/watershed_new'
    seg_key = 'volumes/segmentation/lifted_multicut'

    configs = LiftedMulticutSegmentationWorkflow.get_config()

    config_folder = './config_mc'
    if not os.path.exists(config_folder):
        os.mkdir(config_folder)

    shebang = "#! /g/kreshuk/pape/Work/software/conda/miniconda3/envs/cluster_env37/bin/python"
    global_config = configs['global']
    global_config.update({'shebang': shebang,
                          'block_shape': [256, 256, 256]})
    with open(os.path.join(config_folder, 'global.config'), 'w') as f:
        json.dump(global_config, f)

    ws_config = configs['watershed']
    ws_config.update({'threshold': .6, 'alpha': .9,
                      'sigma_weights': 1.6, 'sigma_seeds': 1.,
                      'apply_dt_2d': False, 'apply_ws_2d': False,
                      'non_maximum_suppresion': True, 'two_pass': False,
                      'halo': [25, 25, 25], 'size_filter': 25,
                      'time_limit': 600, 'mem_limit': 6})
    with open(os.path.join(config_folder, 'watershed.config'), 'w') as f:
        json.dump(ws_config, f)

    agglo_config = configs['agglomerate']
    agglo_config.update({'threshold': .5, 'use_mala_agglomeration': False})
    with open(os.path.join(config_folder, 'agglomerate.config'), 'w') as f:
        json.dump(agglo_config, f)

    labeling_config = configs['find_labeling']
    labeling_config.update({'time_limit': 360, 'threads_per_job': max_threads,
                            'qos': 'normal', 'mem_limit': 32})
    with open(os.path.join(config_folder, 'find_labeling.config'), 'w') as f:
        json.dump(labeling_config, f)

    write_config = configs['write']
    write_config.update({'qos': 'normal', 'mem_limit': 6})
    with open(os.path.join(config_folder, 'write.config'), 'w') as f:
        json.dump(write_config, f)

    graph_config = configs['initial_sub_graphs']
    graph_config.update({'qos': 'normal', 'mem_limit': 4})

    subprob_config = configs['solve_lifted_subproblems']
    subprob_config.update({'threads_per_job': max_threads,
                           'time_limit': 720,
                           'mem_limit': 128,
                           'qos': 'normal',
                           'time_limit_solver': 60*60*6})
    with open(os.path.join(config_folder, 'solve_lifted_subproblems.config'), 'w') as f:
        json.dump(subprob_config, f)

    feat_config = configs['block_edge_features']
    feat_config.update({'offsets': None})
    with open(os.path.join(config_folder, 'block_edge_features.config'), 'w') as f:
        json.dump(feat_config, f)

    weight_edges = True
    exponent = 1.
    beta = .5
    costs_config = configs['probs_to_costs']
    costs_config.update({'weight_edges': weight_edges,
                         'weighting_exponent': exponent,
                         'mem_limit': 16, 'qos': 'high',
                         'beta': beta})
    with open(os.path.join(config_folder, 'probs_to_costs.config'), 'w') as f:
        json.dump(costs_config, f)

    # set number of threads for sum jobs
    tasks = ['merge_sub_graphs', 'merge_edge_features', 'map_edge_ids',
             'reduce_lifted_problem', 'solve_lifted_global']

    for tt in tasks:
        config = configs[tt]
        config.update({'threads_per_job': max_threads,
                       'mem_limit': 256,
                       'time_limit': 1440,
                       'qos': 'normal',
                       'time_limit_solver': 60*60*15})
        with open(os.path.join(config_folder, '%s.config' % tt), 'w') as f:
            json.dump(config, f)

    tasks = ['block_node_labels', 'merge_node_labels']
    for tt in tasks:
        config = configs[tt]
        config.update({"time_limit": 160, "mem_limit": 16})
        with open(os.path.join(config_folder, '%s.config' % tt), 'w') as f:
            json.dump(config, f)

    conf = configs['sparse_lifted_neighborhood']
    conf.update({'time_limit': 240, 'mem_limit': 256, 'threads_per_job': max_threads})
    with open(os.path.join(config_folder, 'sparse_lifted_neighborhood.config'), 'w') as f:
        json.dump(conf, f)

    lifted_labels_path = input_path
    lifted_labels_key = 'volumes/predictions/classes/flagella_and_microvilli'

    n_scales = 1
    task = LiftedMulticutSegmentationWorkflow(input_path=input_path, input_key=input_key,
                                              ws_path=input_path, ws_key=ws_key,
                                              problem_path=exp_path,
                                              node_labels_key='node_labels/lmc_s%i' % n_scales,
                                              output_path=input_path,
                                              output_key=seg_key,
                                              n_scales=n_scales,
                                              lifted_labels_path=lifted_labels_path,
                                              lifted_labels_key=lifted_labels_key,
                                              lifted_prefix='flagella_microvilli',
                                              nh_graph_depth=4,
                                              config_dir=config_folder,
                                              tmp_folder=tmp_folder,
                                              target=target,
                                              max_jobs=max_jobs,
                                              max_jobs_multicut=max_jobs_mc,
                                              skip_ws=False,
                                              agglomerate_ws=True,
                                              mode='same')
    ret = luigi.build([task], local_scheduler=True)
    assert ret, "Failure"


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

    tmp_folder = './tmp_mc'
    run_mc(max_jobs, max_threads, tmp_folder, target=target, max_jobs_mc=max_jobs_mc)
