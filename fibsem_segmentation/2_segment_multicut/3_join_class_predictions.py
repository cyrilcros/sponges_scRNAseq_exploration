#! /g/kreshuk/pape/Work/software/conda/miniconda3/envs/cluster_env37/bin/python

import numpy as np
from concurrent import futures

import z5py
import nifty.tools as nt


def join_classes():
    path = '/g/kreshuk/data/arendt/sponge/data.n5'
    key1 = 'volumes/predictions/classes/flagella'
    key2 = 'volumes/predictions/classes/microvillae'
    out_key = 'volumes/predictions/classes/flagella_and_microvilli'

    f = z5py.File(path)
    ds1 = f[key1]
    ds2 = f[key2]

    shape = ds1.shape
    block_shape = list(ds1.chunks)

    ds_out = f.create_dataset(out_key, shape=shape, chunks=tuple(block_shape),
                              compression='gzip', dtype='uint64')

    id_offset = ds1.attrs['maxId']
    blocking = nt.blocking([0, 0, 0], shape, block_shape)
    n_blocks = blocking.numberOfBlocks

    def join_block(block_id):
        print(block_id, "/", n_blocks)
        block = blocking.getBlock(block_id)
        bb = tuple(slice(beg, end) for beg, end in zip(block.begin, block.end))
        flag = ds1[bb]
        mvil = ds2[bb]

        mvil_mask = mvil != 0
        if np.sum(mvil_mask) == 0 and np.sum(flag != 0) == 0:
            return

        mvil[mvil_mask] += id_offset

        out = flag
        mask = flag == 0
        out[mask] = mvil[mask]
        ds_out[bb] = out

    with futures.ThreadPoolExecutor(8) as tp:
        tasks = [tp.submit(join_block, block_id) for block_id in range(n_blocks)]
        [t.result() for t in tasks]


if __name__ == '__main__':
    join_classes()
