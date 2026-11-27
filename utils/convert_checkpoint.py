#!/usr/bin/env python
'''Convert a checkpoint into a .vtkhdf save/checkpoint pair that DES restarts from.

usage: python utils/convert_checkpoint.py MODEL FRAME -o DIR

MODEL is the run's output path prefix, without an extension; FRAME is the frame of
the checkpoint. The input is a des-binary pair (MODEL.save.NNNNNN and
MODEL.chkpt.NNNNNN, without an extension) or a .vtkhdf pair. A checkpoint without
init_elem_size_n, written by a build without MMG, also gets the records that the
MMG remesher reads:
  init_elem_size_n     the volume-weighted mean size of the elements around each
                       node, as initialize_elem_size_n() in remeshing.cxx computes
                       it at step 0, here on the checkpoint's mesh: the sizes that
                       the last remesh made, or the first mesh's, as deformed since.
                       A size target that an MMG build sets at step 0 (as for a
                       weak zone, in adopt_initial_mmg_output() in mesh.cxx) is not
                       in such a checkpoint;
  refined_elem_size_n  equal to init_elem_size_n, as init() sets it;
  surface dh           zero.
The pair and MODEL.info go to DIR, and DES restarts from them with
  restarting_from_modelname = DIR/<the last part of MODEL>
  restarting_from_frame = FRAME
A des-binary save frame becomes a file of named arrays that only the restart reads;
2vtk-legacy.py converts the frame itself.
'''
import argparse
import os
import shutil
import sys

import h5py
import numpy as np

HEADERLEN = 4096    # the ASCII header of a des-binary file
REVISION = 4        # FILE_REVISION in binaryio.hpp: the layout that both formats share

I4, U4, F8 = np.int32, np.uint32, np.float64


def die(msg):
    sys.exit('Error: ' + msg)


def layout(ndims):
    '''Each des-binary record: its C++ type and its column count (0: one value per row),
    as the .vtkhdf writers give the record of the same name: Output::_write() and
    Output::write_checkpoint() in output.cxx, MarkerSet::write_save_file() and
    write_chkpt_file() in markerset.cxx. Only '<set> size' is des-binary's own (read_binary
    splits it). A key that starts with '.' is a field of a marker set, '<set>.<field>'.'''
    npe, nstr = ndims + 1, ndims * (ndims + 1) // 2
    rec = {'connectivity': (I4, npe), 'segment': (I4, ndims), 'segflag': (I4, 1),
           'bcflag': (U4, 0), '.coord': (F8, ndims), '.eta': (F8, npe)}
    for name in ('coordinate', 'coord0', 'velocity', 'velocity averaged', 'force',
                 'PT residual force', 'displacement'):
        rec[name] = (F8, ndims)
    for name in ('strain', 'strain-rate', 'stress', 'stress averaged'):
        rec[name] = (F8, nstr)
    for name in ('steps', 'nnode', 'nelem', 'nseg', 'info_display_next_step',
                 'dt controls version', 'monitor point count', 'monitor point elements',
                 'PT skip scale update', 'initial equilibrium done', 'output start step',
                 'next regular frame', 'output step interval', 'nremesh',
                 'dv surface acc facets', 'vtopo_len', 'vtopo_int', 'vtopo_dbl_len',
                 '.elem', '.mattype', '.id', '.genesis'):
        rec[name] = (I4, 0)
    for name in ('walltime_sec', 'time_sec', 'dt_sec', 'temperature', 'pore pressure',
                 'radiogenic source', 'plastic strain', 'plastic strain-rate', 'viscosity',
                 'density', 'mesh quality', 'material', 'dynamic friction coefficient',
                 'friction state variable', 'time', 'compensation_pressure',
                 'bottom_temperature', 'dt', 'max_global_vel_mag', 'reference_frame_time',
                 'last_remesh_time', 'dt controls', 'monitor point positions',
                 'monitor point weights', 'pending pore pressure increment',
                 'PT initial residual maximum', 'output start time', 'output time interval',
                 'dv surface acc', 'dhacc', 'volume_old', 'plstrain_remesh',
                 'init_elem_size_n', 'refined_elem_size_n', 'surface dh', 'stressyy',
                 'vtopo_dbl', '.time', '.z', '.distance', '.slope'):
        rec[name] = (F8, 0)
    return rec


def read_binary(path):
    '''The ndims and the named arrays of a des-binary file. Its header lists each record's
    name and byte offset in write order; a record ends where the next one starts.'''
    with open(path, 'rb') as f:
        data = f.read()
    lines = data[:HEADERLEN].split(b'\0')[0].decode().splitlines()
    words = dict(w.split('=') for w in lines[0].split()[2:])
    if int(words['revision']) != REVISION:
        die('%s has layout revision %s; this tool reads revision %d'
            % (path, words['revision'], REVISION))
    ndims = int(words['ndims'])
    entries = [(int(line.rsplit('\t', 1)[1]), line.rsplit('\t', 1)[0]) for line in lines[1:]]
    offsets = [off for off, _ in entries]
    assert offsets == sorted(offsets), path     # an empty record shares the next one's offset
    ends = offsets[1:] + [len(data)]
    table = layout(ndims)
    arrays = {}
    for (off, name), end in zip(entries, ends):
        buf = data[off:end]
        if name == 'provenance':        # build and host text; the restart does not read it
            continue
        if name == 'side_prof_sizes':
            die('%s holds side-wall profiles (remeshing_option = 13), which this tool '
                'does not convert' % path)
        if name.endswith(' size'):      # a marker set: nmarkers, or nmarkers, last_id, reserved_space
            counts = np.frombuffer(buf, I4)
            keys = ('.nmarkers',) if counts.size == 1 else ('.nmarkers', '.last_id', '.reserved_space')
            for key, value in zip(keys, counts):
                arrays[name[:-len(' size')] + key] = np.array([value], I4)
            continue
        key = name[name.index('.'):] if '.' in name else name
        if key not in table:
            die('%s has the record "%s", whose layout this tool does not know' % (path, name))
        dtype, cols = table[key]
        a = np.frombuffer(buf, dtype)
        arrays[name] = a.reshape(-1, cols) if cols else a
    return ndims, arrays


def write_hdf5(path, ndims, arrays):
    '''Write each array under its name at the file root, as HDF5Input in binaryio.cxx
    reads it, and assert that every array reads back unchanged.'''
    with h5py.File(path, 'w') as f:
        f.attrs['ndims'] = np.int32(ndims)
        f.attrs['revision'] = np.int32(REVISION)
        for name, a in arrays.items():
            f.create_dataset(name, data=a)
    with h5py.File(path, 'r') as f:
        for name, a in arrays.items():
            assert np.array_equal(f[name][()], a), name


def element_size_records(coord, conn):
    '''The MMG element-size records of a mesh, by initialize_elem_size_n() in
    remeshing.cxx. 0.433 and 0.118 are its sizefactor (remeshing.cxx), rounded as
    there: the match with DES needs the same values.'''
    nnode, ndims = coord.shape
    edges = coord[conn[:, 1:]] - coord[conn[:, :1]]
    if ndims == 2:
        e0, e1 = edges[:, 0], edges[:, 1]
        vol = 0.5 * np.abs(e0[:, 0] * e1[:, 1] - e0[:, 1] * e1[:, 0])
        size = np.sqrt(vol / 0.433)
    else:
        vol = np.abs(np.linalg.det(edges)) / 6
        size = np.cbrt(vol / 0.118)
    num = np.zeros(nnode)
    den = np.zeros(nnode)
    for k in range(conn.shape[1]):
        np.add.at(num, conn[:, k], size * vol)
        np.add.at(den, conn[:, k], vol)
    init = num / den
    return {'init_elem_size_n': init, 'refined_elem_size_n': init.copy(),
            'surface dh': np.zeros(nnode)}


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument('model', help='output path prefix of the run, without an extension')
    p.add_argument('frame', type=int, help='frame of the checkpoint')
    p.add_argument('-o', dest='outdir', required=True, help='directory for the converted pair')
    args = p.parse_args()

    if os.path.abspath(args.outdir) == os.path.abspath(os.path.dirname(args.model) or '.'):
        die('-o must name another directory than the run\'s: the tool never writes in place')
    base = os.path.basename(args.model)
    src = {k: '%s.%s.%06d' % (args.model, k, args.frame) for k in ('save', 'chkpt')}
    dst = {k: os.path.join(args.outdir, '%s.%s.%06d.vtkhdf' % (base, k, args.frame))
           for k in ('save', 'chkpt')}
    is_vtkhdf = os.path.isfile(src['chkpt'] + '.vtkhdf')
    if not is_vtkhdf and not os.path.isfile(src['chkpt']):
        die('neither %s.vtkhdf nor %s exists' % (src['chkpt'], src['chkpt']))
    if is_vtkhdf:
        with h5py.File(src['chkpt'] + '.vtkhdf', 'r') as chkpt:
            if 'init_elem_size_n' not in chkpt and 'side_prof_sizes' in chkpt:
                die('%s.vtkhdf holds side-wall profiles (remeshing_option = 13) without MMG, '
                    'which this tool does not convert' % src['chkpt'])
    os.makedirs(args.outdir, exist_ok=True)

    if is_vtkhdf:
        for k in src:
            shutil.copyfile(src[k] + '.vtkhdf', dst[k])
        with h5py.File(dst['chkpt'], 'r+') as chkpt:
            if 'init_elem_size_n' in chkpt:
                print('%s: .vtkhdf from a build with MMG, copied unchanged' % src['chkpt'])
            else:
                with h5py.File(dst['save'], 'r') as save:
                    extra = element_size_records(save['coordinate'][()], save['connectivity'][()])
                for name, a in extra.items():
                    chkpt.create_dataset(name, data=a)
                print('%s: .vtkhdf without MMG, element-size records added' % src['chkpt'])
    else:
        ndims, save = read_binary(src['save'])
        ndims_c, chkpt = read_binary(src['chkpt'])
        assert ndims_c == ndims, (ndims, ndims_c)
        note = 'with MMG'
        if 'init_elem_size_n' not in chkpt:
            chkpt.update(element_size_records(save['coordinate'], save['connectivity']))
            note = 'without MMG, element-size records added'
        write_hdf5(dst['save'], ndims, save)
        write_hdf5(dst['chkpt'], ndims, chkpt)
        print('%s: des-binary %s, converted' % (src['chkpt'], note))

    if os.path.isfile(args.model + '.info'):
        shutil.copyfile(args.model + '.info', os.path.join(args.outdir, base + '.info'))
    print('Restart with:\n  restarting_from_modelname = %s\n  restarting_from_frame = %d'
          % (os.path.join(args.outdir, base), args.frame))


if __name__ == '__main__':
    main()
