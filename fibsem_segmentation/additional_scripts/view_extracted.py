import h5py
from cremi_tools.viewer.volumina import view


def view_ex():
    out_path = '/g/kreshuk/data/arendt/sponge/training_data/train_cut1_acsem.h5'

    with h5py.File(out_path) as f:
        data = f['data'][:].transpose((1, 2, 3, 0))

    view([data])

view_ex()
