#!/usr/bin/env python3
"""Check retention cell fields across material interfaces and file containers.

Six-point meshes deliberately have nnode == nelem in both dimensions, so
array length cannot identify field association. Three retention laws share
interface nodes. Independent Python laws and saved marker counts define
the element means; VTK-HDF and converted VTU fields must be cell data.
"""

import argparse
import json
import os
from pathlib import Path
import shutil
import sys
import xml.etree.ElementTree as ET

import numpy as np

import check_bishop_cell as cell
from unsaturated_reference import Retention

FIELDS = ('saturation', 'effective saturation', 'relative permeability', 'bulk density')
LAWS = (Retention('van_genuchten', .001, 1.6, .1),
        Retention('brooks_corey', .001, .8, .05), Retention('gardner', .002, sr=.2))


from des_case import reset_boundaries, run_solver

def polygon(ndims):
    if ndims == 2:
        return ('6 2 0 0\n0 0 0\n1 1 0\n2 1 -1\n3 0 -1\n4 .25 -.5\n5 .75 -.5\n'
                '4 1\n0 0 1 32\n1 1 2 2\n2 2 3 16\n3 3 0 1\n0\n1\n0 .42 -.42 0 2\n')
    return ('6 3 0 0\n0 0 0 0\n1 1 0 0\n2 0 1 0\n3 0 0 -1\n'
            '4 .125 .125 -.5\n5 .25 .25 0\n6 1\n'
            '1 0 32\n3 0 1 5\n1 0 32\n3 1 2 5\n1 0 32\n3 2 0 5\n'
            '1 0 4\n3 0 1 3\n1 0 1\n3 0 2 3\n1 0 2\n3 1 2 3\n'
            '0\n1\n0 .12 .17 -.11 0 2\n')


def config(ndims, pressure, active=True):
    cfg = cell.config(ndims, 1, 'moderate', 1)
    cfg['sim'].update(max_steps='1', max_time_in_yr='1', output_step_interval='1', checkpoint_frame_interval='1',
                      has_initial_checkpoint='true', has_marker_output='true')
    cfg['mesh'].update(meshing_option=str(90 if ndims == 2 else 91), meshing_elem_shape='0', poly_filename='cell.poly',
                       xlength='1', ylength='1', zlength='1', resolution='2',
                       remeshing_option='0', max_steiner_factor='0', min_angle='10')
    cfg['control'].update(fixed_dt='.1', has_moving_mesh='false', has_PT='true', PT_option='1')
    cfg['markers'] = dict(init_marker_option='1', markers_per_element='64')
    reset_boundaries(cfg['bc'], vbc_val='0')
    cfg['ic'].update(pore_pressure_value=str(pressure), weakzone_option='0', mattype_option='1',
                      num_mattype_layers='3', layer_mattypes='[0,1,2]', mattype_layer_depths='[.3,.6]')
    for key, text in list(cfg['mat'].items()):
        if text.startswith('['):
            values = json.loads(text)
            if len(values) == 1: cfg['mat'][key] = json.dumps(values*3)
    cfg['mat'].update(num_materials='3', rho0='[2650,2650,2650]', porosity='[.3,.3,.3]',
                      fluid_rho0='[1000,1000,1000]', hydraulic_perm='[0,0,0]', alpha='[0,0,0]',
                      retention_model='[1,2,3]' if active else '[0,0,0]', retention_alpha='[.001,.001,.002]',
                      retention_n='[1.6,.8,2]', residual_saturation='[.1,.05,.2]',
                      is_plane_strain=str(ndims == 2).lower())
    return cfg


def compare(folder, active, convert):
    model = cell.Dynearthsol(str(folder/'cell')); frame = model.frames[-1]
    coordinate = model.read_field(frame, 'coordinate'); conn = model.read_field(frame, 'connectivity')
    assert len(coordinate) == len(conn) == 6, (len(coordinate), len(conn))
    marker = model.read_markers(frame, 'markerset')
    count = np.zeros((len(conn), 3), dtype=int)
    np.add.at(count, (marker['markerset.elem'], marker['markerset.mattype']), 1)
    assert np.all(count.sum(axis=0) > 0) and np.any((count > 0).sum(axis=1) > 1)
    weights = count/count.sum(axis=1)[:, None]
    if active:
        pressure = model.read_field(frame, 'pore pressure')[conn]
        expected = {}
        for name, method in (('saturation', 'sat'), ('effective saturation', 'se'), ('relative permeability', 'kr')):
            per_material = np.stack([getattr(law, method)(pressure).mean(axis=1) for law in LAWS], axis=1)
            expected[name] = np.sum(per_material*weights, axis=1)
        expected['bulk density'] = 2650*.7+1000*.3*expected['saturation']
        for name, values in expected.items():
            actual = model.read_field(frame, name)
            assert actual.shape == (len(conn),) and np.isfinite(actual).all()
            np.testing.assert_allclose(actual, values, rtol=3e-11, atol=1e-300 if name == 'relative permeability' else 1e-11)
    else:
        for name in FIELDS:
            try: model.read_field(frame, name)
            except KeyError: pass
            else: raise AssertionError(f'Disabled retention emits {name}')
    if model.format == 'hdf5':
        import h5py
        with h5py.File(model.get_fn(frame), 'r') as stream:
            cells, points = stream['VTKHDF/grid/CellData'], stream['VTKHDF/grid/PointData']
            for name in ('stress', 'strain', 'density', 'material', 'viscosity', 'radiogenic source'):
                assert name in cells and name not in points, name
            for name in ('pore pressure', 'temperature', 'velocity', 'force', 'coord0'):
                assert name in points and name not in cells, name
            assert 'fluid balance' not in cells and 'fluid balance' not in points
            np.testing.assert_array_equal(stream['fluid balance'][...], model.read_field(frame, 'fluid balance'))
            for name in FIELDS:
                assert (name in cells) == active and name not in points
                if active: np.testing.assert_array_equal(stream['/'+name][...], cells[name][...])
    elif convert:
        env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
        run_solver([sys.executable, str(cell.ROOT / '2vtk.py'), '-a', str(folder / 'cell'), '0', '1'], log=folder / 'vtk.log', cwd=folder, env=env, check=True, timeout=90)
        vtk = ET.parse(folder/f'cell.{frame:06d}.vtu')
        arrays = {a.get('Name'):a for a in vtk.findall('.//CellData/DataArray')}
        for name in FIELDS:
            assert (name in arrays) == active
            if active: np.testing.assert_allclose(np.fromstring(arrays[name].text, sep=' '), expected[name], rtol=1e-6)
    return dict(nodes=len(coordinate), elements=len(conn), mixed_elements=int(np.sum((count > 0).sum(axis=1) > 1)),
                container=model.format, fields=list(FIELDS) if active else [])


def study(executable, output, ndims):
    output.mkdir(parents=True, exist_ok=False)
    binary_hash = cell.digest(executable)
    report = dict(executable=str(executable), executable_sha256=binary_hash, ndims=ndims, cases={})
    for source in (Path(__file__), Path(cell.__file__), cell.ROOT/'benchmarks/poromechanics/unsaturated_reference.py'):
        shutil.copy2(source, output/source.name)
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    for label, pressure, active in (('wet', 1000., True), ('partial', -2000., True),
                                    ('dry', -1e6, True), ('gradient', 0., True), ('disabled', -2000., False)):
        folder = output/label; folder.mkdir()
        (folder/'cell.poly').write_text(polygon(ndims))
        cfg = config(ndims, pressure, active)
        if label == 'gradient':
            cfg['ic']['pore_pressure_option'] = '1'
            cfg['bc'].update(hbc_x0='2', hbc_x1='2', hbc_pressure='[-2000,1000,0,0,0,0]')
        with (folder/'input.cfg').open('w') as stream: cfg.write(stream)
        run_solver([str(executable), 'input.cfg'], log=folder / 'run.log', cwd=folder, env=env, check=True, timeout=120)
        case = compare(folder, active, convert=True)
        case['files_sha256'] = {p.name:cell.digest(p) for p in sorted(folder.iterdir())
                                if p.suffix in ('.cfg', '.poly', '.info', '.vtu') or '.save.' in p.name or '.chkpt.' in p.name}
        report['cases'][label] = case
        print(label, case['container'], case['mixed_elements'], 'mixed elements: PASS', flush=True)
    assert cell.digest(executable) == binary_hash
    report.update(passed=True, runner_sha256=cell.digest(Path(__file__)))
    (output/'qualification.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path); parser.add_argument('output', type=Path)
    parser.add_argument('--ndims', type=int, choices=(2, 3), default=2)
    args = parser.parse_args()
    study(args.executable.resolve(), args.output.resolve(), args.ndims)
