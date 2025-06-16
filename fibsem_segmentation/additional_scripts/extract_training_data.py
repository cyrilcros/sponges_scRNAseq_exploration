import numpy as np
import h5py
import z5py
from cremi_tools.viewer.volumina import view


#
# Input data for stage 0: semantic pixel classification
#


def extract_raw(bb, out_path, check=False):
    path = '/g/kreshuk/data/arendt/sponge/data.n5'
    with z5py.File(path, 'r') as f:
        raw = f['volumes/raw/s0'][bb]

    if check:
        path = '/g/kreshuk/data/arendt/sponge/training_data/old/train_data_scale0.h5'
        with h5py.File(path) as f:
            print(list(f.keys()))
            raw1 = f['volumes/raw'][:]
        view([raw, raw1])
    else:
        with h5py.File(out_path, 'w') as f:
            f.create_dataset('data', data=raw, compression='gzip',
                             chunks=(32, 32, 32))


def extract_cut1():
    central_micron = [31.0, 28.2, 26.4][::-1]
    resolution = [0.015, 0.015, 0.015]
    central = [int(ce / re) for ce, re in zip(central_micron, resolution)]

    halo = (50, 756, 756)
    shift = (-100, 250, 250)
    bb = tuple(slice(ce - ha + sh,
                     ce + ha + sh)
               for ce, ha, sh in zip(central, halo, shift))

    out_path = '/g/kreshuk/data/arendt/sponge/training_data/train_cut1.h5'
    extract_raw(bb, out_path)


def extract_cut2():
    central_micron = [30.4, 43.9, 39.4][::-1]
    resolution = [0.015, 0.015, 0.015]
    central = [int(ce / re) for ce, re in zip(central_micron, resolution)]

    halo = (50, 756, 756)
    shift = (0, 0, 0)
    bb = tuple(slice(ce - ha + sh,
                     ce + ha + sh)
               for ce, ha, sh in zip(central, halo, shift))

    out_path = '/g/kreshuk/data/arendt/sponge/training_data/train_cut2.h5'
    extract_raw(bb, out_path)


#
# Input data for stage 1: semantic autocontext
#


def extract_ac_sem(bb, out_path):
    path = '/g/kreshuk/data/arendt/sponge/data.n5'
    with z5py.File(path, 'r') as f:
        raw = f['volumes/raw/s0']
        raw.n_threads = 8
        raw = raw[bb]
        print(raw.dtype)
        pred = f['volumes/predictions/semantic_stage0']
        pred.n_threads = 8
        pred = pred[(slice(None),) + bb]
        print(pred.dtype)
    raw = raw.astype('float32')
    raw -= raw.min()
    raw /= raw.max()

    print(raw.shape)
    print(pred.shape)
    data = np.concatenate([raw[None], pred], axis=0)
    print(data.shape)

    with h5py.File(out_path, 'w') as f:
        f.create_dataset('data', data=data, compression='gzip',
                         chunks=(1, 32, 32, 32))


def extract_cut1_ac_sem():
    central_micron = [31.0, 28.2, 26.4][::-1]
    resolution = [0.015, 0.015, 0.015]
    central = [int(ce / re) for ce, re in zip(central_micron, resolution)]

    halo = (50, 756, 756)
    shift = (-100, 250, 250)
    bb = tuple(slice(ce - ha + sh,
                     ce + ha + sh)
               for ce, ha, sh in zip(central, halo, shift))

    out_path = '/g/kreshuk/data/arendt/sponge/training_data/train_cut1_acsem.h5'
    extract_ac_sem(bb, out_path)


def extract_cut2_ac_sem():
    central_micron = [30.4, 43.9, 39.4][::-1]
    resolution = [0.015, 0.015, 0.015]
    central = [int(ce / re) for ce, re in zip(central_micron, resolution)]

    halo = (50, 756, 756)
    shift = (0, 0, 0)
    bb = tuple(slice(ce - ha + sh,
                     ce + ha + sh)
               for ce, ha, sh in zip(central, halo, shift))

    out_path = '/g/kreshuk/data/arendt/sponge/training_data/train_cut2_acsem.h5'
    extract_ac_sem(bb, out_path)

#
# Input data for stage 2: binary autocontext
#


def extract_ac_bin(bb, out_path):
    path = '/g/kreshuk/data/arendt/sponge/data.n5'
    with z5py.File(path, 'r') as f:
        raw = f['volumes/raw/s0']
        raw.n_threads = 8
        raw = raw[bb]
        print(raw.dtype)
        pred = f['volumes/predictions/semantic_stage1']
        pred.n_threads = 8
        pred = pred[(slice(None),) + bb]
        print(pred.dtype)
    pred *= 255
    pred = pred.astype('uint8')
    print(pred.min(), pred.max())

    print(raw.shape)
    print(pred.shape)
    data = np.concatenate([raw[None], pred], axis=0)
    print(data.shape)

    with z5py.File(out_path, 'w') as f:
        f.create_dataset('data', data=data, compression='gzip',
                         chunks=(1, 64, 64, 64))


def extract_cut1_ac_bin():
    central_micron = [31.0, 28.2, 26.4][::-1]
    resolution = [0.015, 0.015, 0.015]
    central = [int(ce / re) for ce, re in zip(central_micron, resolution)]

    halo = (50, 756, 756)
    shift = (-100, 250, 250)
    bb = tuple(slice(ce - ha + sh,
                     ce + ha + sh)
               for ce, ha, sh in zip(central, halo, shift))

    out_path = '/g/kreshuk/data/arendt/sponge/training_data/train_cut1_acbin.n5'
    extract_ac_bin(bb, out_path)


def extract_cut2_ac_bin():
    central_micron = [30.4, 43.9, 39.4][::-1]
    resolution = [0.015, 0.015, 0.015]
    central = [int(ce / re) for ce, re in zip(central_micron, resolution)]

    halo = (50, 756, 756)
    shift = (0, 0, 0)
    bb = tuple(slice(ce - ha + sh,
                     ce + ha + sh)
               for ce, ha, sh in zip(central, halo, shift))

    out_path = '/g/kreshuk/data/arendt/sponge/training_data/train_cut2_acbin.n5'
    extract_ac_bin(bb, out_path)


if __name__ == '__main__':
    # extract_cut1()
    # extract_cut2()

    # extract_cut1_ac_sem()
    # extract_cut2_ac_sem()

    extract_cut1_ac_bin()
    extract_cut2_ac_bin()
