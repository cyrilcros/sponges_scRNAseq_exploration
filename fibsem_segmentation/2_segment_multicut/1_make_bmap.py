#! /g/kreshuk/pape/Work/software/conda/miniconda3/envs/cluster_env37/bin/python

import z5py
import nifty.tools as nt
from concurrent import futures


def make_bmap():
    path = '/g/kreshuk/data/arendt/sponge/data.n5'
    in_key = 'volumes/predictions/binary_stage2'
    out_key = 'volumes/predictions/boundaries'
    f = z5py.File(path)
    ds_in = f[in_key]
    chunks = ds_in.chunks[1:]
    ds_out = f.create_dataset(out_key, shape=ds_in.shape[1:], chunks=chunks,
                              dtype=ds_in.dtype, compression='gzip')

    blocking = nt.blocking([0, 0, 0], ds_out.shape, list(chunks))

    def for_chunk(block_id):
        print(block_id, "/", blocking.numberOfBlocks)
        block = blocking.getBlock(block_id)
        bb = tuple(slice(beg, end) for beg, end in zip(block.begin, block.end))
        bb_in = (slice(1, 2),) + bb
        data = ds_in[bb_in].squeeze()
        ds_out[bb] = data

    with futures.ThreadPoolExecutor(8) as tp:
        tasks = [tp.submit(for_chunk, block_id)
                 for block_id in range(blocking.numberOfBlocks)]
        [t.result() for t in tasks]


if __name__ == '__main__':
    make_bmap()
