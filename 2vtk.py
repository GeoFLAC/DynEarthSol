#!/usr/bin/env python
# encoding: utf-8
'''Convert the binary output of DynEarthSol to VTK files.

usage: 2vtk.py [options] modelname [start [end [delta]]]

A model that was not restarted: its .save frames become .vtu/.vtp files beside
them, or with -c in the current directory.
    2vtk.py model               every frame
    2vtk.py model 10 20 2       frames 10, 12, ..., 20 (positions, see below)
    2vtk.py model -1            resume after the last .vtu converted

A restarted model: it holds its frames from its restart frame on, and its
.manifest names the run it restarted from, its parent, which holds the earlier
frames (and may be restarted itself).
    2vtk.py model               this model's own frames only
    2vtk.py -g model            the whole series: the parents' frames are
                                converted beside their own frames (-c: in the
                                current directory) and linked under this
                                model's name
    2vtk.py -g -copy model      the whole series, the parents' frames written
                                under this model's name instead of linked

.vtkhdf frames (an HDF5 build): ParaView reads them as they are. Instead of
writing .vtu files, -u or -hdf adds the calculated fields (requires h5py).
    2vtk.py -u model            into the frames themselves, listed in
                                model.save.vtkhdf.series for ParaView
    2vtk.py -hdf model          into copies, model.NNNNNN.vtkhdf, listed in
                                model.vtkhdf.series; the frames are unchanged
    2vtk.py -u -g model         a restarted model's whole series, the parents'
    2vtk.py -hdf -g model       frames linked under this model's name
                                (-hdf -g -copy: copies under it instead)

options:
    -a          save data in ASCII format (default: binary)
    -c          work in the current directory: the files, links and .series
                this run writes go there (default: the model's own directory)
    -copy       with -g, write the parents' frames under this model's name,
                not linked (not with -u, which updates frames in place)
    -g          also convert (-u: update) the parents' frames, and link them
                under this model's name
    -hdf        write each .vtkhdf frame with calculated fields as a copy
    -m          save marker data
    -p          save principal components (s1 and s3) of deviatoric stress
    -t          save all tensor components (default: only 1st/2nd invariants)
    -u          update the .vtkhdf frames in place with calculated fields,
                also with -c (-hdf writes copies instead)
                WARNING: Do not use this option while the simulation is running
                or accessing the files, as it may corrupt the data.
    -h,--help   show this help

'start', 'end' and 'delta' are positions among the frames converted, from 0,
not frame numbers: a restarted model's own frames begin at its restart frame,
and with -g the whole series is counted. 'end' is inclusive, and -1 or none is
the last; 'start' -1 resumes after the last .vtu converted, from the first if
there is none.

Where each frame goes, W being the work directory (the current one with -c):

               this model's frames     each parent's frames (-g)
  (default)    W/model.NNNNNN.vtu      parent.NNNNNN.vtu beside the parent's
                                       frames (-c: in W), linked as W/model.*
  -copy        as above                W/model.NNNNNN.vtu, not linked
  -hdf         W/model.NNNNNN.vtkhdf   parent.NNNNNN.vtkhdf beside the parent's
                                       frames (-c: in W), linked as
                                       W/model.NNNNNN.vtkhdf
  -hdf -copy   as above                W/model.NNNNNN.vtkhdf, not linked
  -u           updated in place        updated in place, linked as
                                       W/model.save.NNNNNN.vtkhdf

With -u -g, this model's own frames are linked into W as well when W is not
their directory, and a real frame already at a link's name is kept. A parent
with no frame on disk before the restart frame, as when only the restart
frame's .save and .chkpt files were copied for the restart, ends the series
there, with a warning. Caution: a -u -g link has a .save name, and a DES run
writing that frame in W (a fresh run of this model there) truncates the
linked frame through it.
'''

from __future__ import print_function, unicode_literals
import sys, os, shutil
import base64, zlib, glob, itertools, json
import numpy as np

# Disable HDF5 file locking to avoid BlockingIOError on some filesystems
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'

try:
    import h5py
except ImportError:
    h5py = None
from scipy import spatial
from numpy.linalg  import eigh
from fractions import Fraction
from Dynearthsol import Dynearthsol

# Save in ASCII or encoded binary.
# Some old VTK programs cannot read binary VTK files.
output_in_binary = True

# Save the resultant vtu files in current directory?
output_in_cwd = False

# Save indivisual components?
output_tensor_components = True

# Save principle stresses
output_principle_stress = False

# Save markers?
output_markers = True

# Calculate melting?
output_melting = False

# Calculate heat flow?
output_heatflux = False
conductivity = 3.3

# Update the .vtkhdf frames in place (-u)?
update_vtkhdf = False

# Write updated copies as modelname.NNNNNN.vtkhdf instead (-hdf)?
hdf_copies = False

# Link a parent run's output under this model's name (-copy: write it under this name)?
link_parent_frames = True

# Convert or update the parent runs' frames too (-g)?
generate_parent_frames = False

# min mutiprocessing threads
mutiprocessing_threads = 4

# Options the workers read, passed to them: spawned workers (the macOS default) skip __main__.
WORKER_OPTIONS = ('output_in_binary', 'hdf_copies', 'output_tensor_components',
                  'output_principle_stress', 'output_markers', 'output_melting')

########################
# Is numpy version < 1.8?
eigh_vectorized = True
npversion = np.__version__.split('.')
npmajor = int(npversion[0])
npminor = int(npversion[1])
if npmajor < 1 or (npmajor == 1 and npminor < 8):
    eigh_vectorized = False

class Filter():

  def marker(self,par,x, z, m, t):
    ind = (par.xmin <= x) * (x <= par.xmax) * \
          (par.zmin <= z) * (z <= par.zmax)
    x = x[ind]
    z = z[ind]
    m = m[ind]
    t = t[ind]
    return x, z, m, t

  def node(self,x,z,f):
    ind = (f >= 32 ) * (f <= 34)
    x = x[ind]
    z = z[ind]
    return x, z

########################
# Is numpy version < 1.8?
eigh_vectorized = True
npversion = np.__version__.split('.')
npmajor = int(npversion[0])
npminor = int(npversion[1])
if npmajor < 1 or (npmajor == 1 and npminor < 8):
    eigh_vectorized = False


def calculate_derived_data(des, frame):
    point_data = {}
    cell_data = {}
    
    # Node-based
    coord = des.read_field(frame, 'coordinate')
    coord0 = des.read_field(frame, 'coord0')
    nnode = coord.shape[0]
    
    # total displacement
    disp = np.zeros((nnode, 3), dtype=coord.dtype)
    disp[:,0:des.ndims] = coord - coord0
    point_data['total displacement'] = (disp, 3)
    
    # horizon
    horizon = np.zeros((nnode), dtype=coord.dtype)
    horizon[:] = coord0[:,-1]
    point_data['horizon'] = (horizon, 1)

    # location
    map0 = np.zeros((nnode), dtype=coord.dtype)
    map0[:] = coord0[:,0]
    point_data['location_x'] = (map0, 1)
    if des.ndims == 3:
        map1 = np.zeros((nnode), dtype=coord.dtype)
        map1[:] = coord0[:,1]
        point_data['location_y'] = (map1, 1)    
    
    # Element-based
    # Strain Rate
    strain_rate = des.read_field(frame, 'strain-rate')
    srII = second_invariant(strain_rate)
    cell_data['strain-rate II log10'] = (np.log10(srII+1e-45), 1)
    
    if output_tensor_components:
        for d in range(des.nstr):
            cell_data['strain-rate ' + des.component_names[d]] = (strain_rate[:,d], 1)

    # Strain
    strain = des.read_field(frame, 'strain')
    sI = first_invariant(strain)
    sII = second_invariant(strain)
    cell_data['strain I'] = (sI, 1)
    cell_data['strain II'] = (sII, 1)
    
    if output_tensor_components:
        for d in range(des.nstr):
            cell_data['strain ' + des.component_names[d]] = (strain[:,d], 1)

    # Stress
    try:
        stress = des.read_field(frame, 'stress averaged')
    except KeyError:
        stress = des.read_field(frame, 'stress')
    tI = first_invariant(stress)
    tII = second_invariant(stress)
    cell_data['stress I'] = (tI, 1)
    cell_data['stress II'] = (tII, 1)
    
    if output_tensor_components:
        for d in range(des.ndims):
            cell_data['stress ' + des.component_names[d]] = (stress[:,d], 1)
        for d in range(des.ndims, des.nstr):
            cell_data['stress ' + des.component_names[d]] = (stress[:,d], 1)

    # Principal Stress
    if output_principle_stress:
        s1, s3 = compute_principal_stress(stress)
        cell_data['s1'] = (s1, 3)
        cell_data['s3'] = (s3, 3)

    # Effective Viscosity
    effvisc = tII / (srII + 1e-45)
    cell_data['effective viscosity'] = (effvisc, 1)

    # Melting
    if output_melting:
        material = des.read_field(frame, 'material')
        temperature = des.read_field(frame, 'temperature')
        connectivity = des.read_field(frame, 'connectivity')
        nelem = des.nelem_list[des.frames.index(frame)]
        
        # Optimization: Vectorized mean
        ecoord = coord[connectivity].mean(axis=1)
        etemp = temperature[connectivity].mean(axis=1)

        melting = np.zeros(sI.shape)

        # find surface
        bcflag = des.read_field(frame, 'bcflag')
        filter = Filter()
        surfx, surfz = filter.node(coord[:,0],coord[:,1],bcflag)
        orders = np.argsort(surfx)
        surface = np.vstack((surfx[orders],surfz[orders]))
        
        depth = np.interp(ecoord[:,0], surface[0], surface[1]) - ecoord[:,1]
        pressure = depth * 9.8 * 2900.
        
        melting[:] = -1000
        ind = material < 2
        melting[ind] = (etemp[ind]-273. + depth[ind]*3.e-4) - (1120 + (680./7.e9)*pressure[ind])
        cell_data['melting'] = (melting, 1)
        
    return point_data, cell_data

def unlinked(filename):
    '''filename with any link there removed, so writing it replaces a parent's file.'''
    if os.path.islink(filename):
        os.remove(filename)
    return filename


def process_single_frame(args):
    des, output_prefix, i = args

    frame = des.frames[i]
    nnode = des.nnode_list[i]
    nelem = des.nelem_list[i]
    step = des.steps[i]
    time_in_yr = des.time[i] / (365.2425 * 86400)

    des.read_header(frame)
    suffix = '{0:0=6}'.format(frame)

    filename = '{0}.{1}.vtu'.format(output_prefix, suffix)
    fvtu = open(unlinked(filename), 'w')

    try:
        vtu_header(fvtu, nnode, nelem, time_in_yr, step)

        #
        # node-based field
        #
        fvtu.write('  <PointData>\n')

        # averaged velocity is more stable and is preferred
        try:
            convert_field(des, frame, 'velocity averaged', fvtu)
        except KeyError:
            convert_field(des, frame, 'velocity', fvtu)

        convert_field(des, frame, 'force', fvtu)
        
        # Calculate derived data
        point_data, cell_data = calculate_derived_data(des, frame)
        
        # Write Point Data
        for name, (data, comps) in point_data.items():
            vtk_dataarray(fvtu, data, name, comps)

        '''
        # find nearest neighbour marker of nodes
        markersetname = 'markerset'
        marker_data = des.read_markers(frame, markersetname)
        nmarkers = marker_data['size']
        if nmarkers <= 0:
            raise MarkerSizeError()
        marker_coord = marker_data[markersetname + '.coord']
        marker_mattype = marker_data[markersetname + '.mattype']
        kdtree = spatial.KDTree(marker_coord)
        nn = kdtree.query(coord,1)
        nnmattype = np.zeros((nnode), dtype=marker_mattype.dtype)

        try:
            marker_time = marker_data[markersetname + '.time']
            nnchron = np.zeros((nnode), dtype=marker_time.dtype)
            nnchron[:] = marker_time[nn[1]]
        except:
            pass

        # abjust horizon of sediment node
        # Note: horizon is now calculated in calculate_derived_data, but we need to read it back if we want to modify it?
        # The original code here was modifying 'horizon' variable.
        # But since the sediment adjustment block is commented out, I will ignore it.
        
        try:
            vtk_dataarray(fvtu, nnchron, 'chron', 1)
        except:
            pass
        '''

        convert_field(des, frame, 'temperature', fvtu)
        convert_field(des, frame, 'bcflag', fvtu)
        #convert_field(des, frame, 'mass', fvtu)
        #convert_field(des, frame, 'tmass', fvtu)
        #convert_field(des, frame, 'volume_n', fvtu)
        convert_field(des, frame, 'pore pressure', fvtu)

        # node number for debugging
        vtk_dataarray(fvtu, np.arange(nnode, dtype=np.int32), 'node number')

        fvtu.write('  </PointData>\n')
        #
        # element-based field
        #
        fvtu.write('  <CellData>\n')

        #convert_field(des, frame, 'volume', fvtu)
        #convert_field(des, frame, 'edvoldt', fvtu)

        convert_field(des, frame, 'mesh quality', fvtu)
        convert_field(des, frame, 'plastic strain', fvtu)
        convert_field(des, frame, 'plastic strain-rate', fvtu)
        try:
            convert_field(des, frame, 'radiogenic source', fvtu)
        except KeyError:
            # Field not present in this dataset; skip it.
            pass
        except NameError:
            # Optional field or missing symbol; report and continue without radiogenic source.
            print(
                "Warning: 'radiogenic source' field not written for frame {} "
                "because required name is not defined.".format(frame),
                file=sys.stderr,
            )

        # Optional RSF cell fields.
        try:
            convert_field(des, frame, 'dynamic friction coefficient', fvtu)
        except (KeyError, NameError):
            # Optional RSF field not present or not defined in this dataset; skip it.
            pass

        try:
            convert_field(des, frame, 'friction state variable', fvtu)
        except (KeyError, NameError):
            # Optional RSF friction state variable field not present; skip it.
            if i == 0:
                print(
                    "Info: 'friction state variable' field is not available in this dataset"
                    " and will be skipped.",
                    file=sys.stderr,
                )

        # Write Cell Data
        for name, (data, comps) in cell_data.items():
            vtk_dataarray(fvtu, data, name, comps)

        convert_field(des, frame, 'density', fvtu)
        convert_field(des, frame, 'material', fvtu)
        convert_field(des, frame, 'viscosity', fvtu)
        
        # element number for debugging
        vtk_dataarray(fvtu, np.arange(nelem, dtype=np.int32), 'elem number')

        # # heat flux
        # # 3D is not implemented and tested yet
        # if output_heatflux:               
        #     flux, flux_val = des.load_calculation(frame, 'heat flux')

        #     vtk_dataarray(fvtu, flux[0], 'heat flux x')
        #     if des.ndims == 3:
        #         vtk_dataarray(fvtu, flux[1], 'heat flux y')
        #     vtk_dataarray(fvtu, flux[-1], 'heat flux z')
        #     vtk_dataarray(fvtu, flux_val, 'heat flux magnitude')

        fvtu.write('  </CellData>\n')

        #
        # node coordinate
        #
        fvtu.write('  <Points>\n')
        convert_field(des, frame, 'coordinate', fvtu)
        fvtu.write('  </Points>\n')

        #
        # element connectivity & types
        #
        fvtu.write('  <Cells>\n')
        convert_field(des, frame, 'connectivity', fvtu)
        vtk_dataarray(fvtu, (des.ndims+1)*np.array(range(1, nelem+1), dtype=np.int32), 'offsets')
        if des.ndims == 2:
            # VTK_ TRIANGLE == 5
            celltype = 5
        else:
            # VTK_ TETRA == 10
            celltype = 10
        vtk_dataarray(fvtu, celltype*np.ones((nelem,), dtype=np.int32), 'types')
        fvtu.write('  </Cells>\n')

        vtu_footer(fvtu)
        fvtu.close()

    except:
        # delete partial vtu file
        fvtu.close()
        os.remove(filename)
        raise

    #
    # Converting marker
    #
    if output_markers:
        # ordinary markerset
        filename = '{0}.{1}.vtp'.format(output_prefix, suffix)
        output_vtp_file(des, frame, filename, 'markerset', time_in_yr, step)

        # hydrous markerset
        if 'hydrous-markerset size' in des.field_pos:
            filename = '{0}.hyd-ms.{1}.vtp'.format(output_prefix, suffix)
            output_vtp_file(des, frame, filename, 'hydrous-markerset', time_in_yr, step)

    return suffix


def process_vtkhdf_update(args):
    des, output_prefix, i = args
    frame = des.frames[i]
    
    src_filename = des.get_fn(frame)
    assert src_filename.endswith('.vtkhdf'), src_filename   # main() refuses des-binary runs

    suffix = '{0:0=6}'.format(frame)

    if hdf_copies:
        # -hdf: update a copy, named like the .vtu output
        target_filename = '{0}.{1}.vtkhdf'.format(output_prefix, suffix)
        if os.path.abspath(src_filename) != os.path.abspath(target_filename):
            shutil.copy2(src_filename, unlinked(target_filename))
    else:
        # In-place update
        target_filename = src_filename

    # Now open target_filename in r+ mode
    try:
        with h5py.File(target_filename, 'r+') as f:
            # Helper to create dataset and link
            def write_dataset(path, data, name_attr=None):
                # path e.g. /VTKHDF/grid/CellData/stress II
                if path in f:
                    del f[path]
                dset = f.create_dataset(path, data=data, compression="gzip", compression_opts=9, shuffle=True)
                if name_attr:
                    dset.attrs['Name'] = np.bytes_(name_attr) # VTK expects string attributes? or just Name?
                
                # Create soft link at root
                # e.g. /stress II -> /VTKHDF/grid/CellData/stress II
                link_name = '/' + os.path.basename(path)
                if link_name in f:
                    del f[link_name]
                f[link_name] = h5py.SoftLink(path)

            # Calculate derived data
            point_data, cell_data = calculate_derived_data(des, frame)
            
            # Write Point Data
            for name, (data, comps) in point_data.items():
                write_dataset('/VTKHDF/grid/PointData/' + name, data, name)
                
            # Write Cell Data
            for name, (data, comps) in cell_data.items():
                write_dataset('/VTKHDF/grid/CellData/' + name, data, name)

    except Exception as e:
        print(f"Error processing frame {frame}: {e}", file=sys.stderr)
        raise

    return suffix


def read_run_records(modelname):
    '''[runtime.model] records of modelname.manifest, oldest first ([] without one).
    Each restart appends one; one after the "could not start" seam wrote no frames.'''
    records, rec, failed = [], None, False
    if not os.path.isfile(modelname + '.manifest'):
        return records
    with open(modelname + '.manifest') as f:
        for line in f:
            # headers and seams begin at column 0; no line of the embedded code diff begins with '['
            if line.startswith('# ---- a run that could not start'):
                failed = True
            elif line.startswith('['):
                rec = None
                if line.strip() == '[runtime.model]':
                    rec = {}
                    if not failed:
                        records.append(rec)
                    failed = False
            elif rec is not None and '=' in line:
                key, value = line.split('=', 1)
                rec[key.strip()] = value.strip()
    return records


def restart_origin(modelname):
    '''(parent, first frame of its own) of the run whose frames modelname holds, (None, 0)
    if fresh; a same-name resume keeps its earlier frames, so the record before it decides.'''
    for rec in reversed(read_run_records(modelname)):
        if rec['restarting'] != 'yes':
            return None, 0
        # DES resolved both names against its cwd, and '..' after symlinks as realpath does
        cwd = modelname
        for _ in os.path.normpath(rec['modelname']).split(os.sep):
            cwd = os.path.dirname(cwd)
        parent = os.path.realpath(os.path.join(cwd, rec['restart_from_model']))
        if parent != os.path.realpath(modelname):
            return os.path.relpath(parent), int(rec['restart_from_frame'])
    return None, 0


def saved_frames(modelname):
    '''Frame numbers of the modelname.save.* files on disk, whatever its .info lists.'''
    frames = set()
    for fn in glob.glob(modelname + '.save.*'):
        num = fn[len(modelname + '.save.'):].split('.')[0]
        if num.isdigit():
            frames.add(int(num))
    return frames


def restart_series(modelname, with_parents):
    '''(des, frame index) of the frames to convert, oldest first: modelname's own, from its
    restart frame on, and with_parents the earlier ones of every run up its restart chain.'''
    series, stop, seen = [], None, set()
    while True:
        if os.path.realpath(modelname) in seen:   # restarted from its own descendant
            print(f'Warning: the restart chain returns to {modelname}; the frames before '
                  f'{stop} are not converted.', file=sys.stderr)
            return series
        seen.add(os.path.realpath(modelname))
        des = Dynearthsol(modelname)
        parent, first = restart_origin(modelname)
        series[:0] = [(des, i) for i, f in enumerate(des.frames)
                      if f >= first and (stop is None or f < stop)]
        if parent is None:
            return series
        # Frames on disk only: a copied restart frame (or .info) is not the parent run.
        parent_frames = any(f < first for f in saved_frames(parent))
        if not with_parents:
            if parent_frames:
                print(f'Note: {modelname} restarted from frame {first} of {parent}; -g adds its frames.',
                      file=sys.stderr)
            return series
        if not parent_frames:
            print(f'Warning: {modelname} restarted from frame {first} of {parent}, which has no '
                  f'frame before {first} here; those frames are not converted.', file=sys.stderr)
            return series
        print(f'Series: {modelname} restarted from frame {first} of {parent}.', file=sys.stderr)
        modelname, stop = parent, first


def init_worker(options):
    globals().update(options)


def prefix_of(modelname):
    '''Where the files named after modelname go: the current directory with -c, else beside
    its frames.'''
    return os.path.basename(modelname) if output_in_cwd else modelname


def announce_model(batch, ndone, prefix, names):
    '''Name the model whose frames batch converts, below the last progress line.'''
    des, own, i = batch[0]
    if ndone:
        print(file=sys.stderr)   # the progress line ends in '\r'; keep it
    linked = f', linked as {prefix}.*' if names and own != prefix else ''
    print(f'Working on model {des.modelname} (frames {des.frames[i]}-{des.frames[batch[-1][2]]}){linked}.',
          file=sys.stderr)


VTU_NAMES = ('{0}.{1}.vtu', '{0}.{1}.vtp', '{0}.hyd-ms.{1}.vtp')
VTKHDF_NAMES = ('{0}.{1}.vtkhdf',)


def link_file(target, link):
    '''Point link at target by a path relative to the link's directory, replacing what is there.'''
    if os.path.lexists(link):
        os.remove(link)
    os.symlink(os.path.relpath(os.path.realpath(target), os.path.realpath(os.path.dirname(link))), link)


def link_frames(batch, prefix, names):
    '''Give prefix's names to the files a parent's batch wrote under its own prefix.'''
    if not names or batch[0][1] == prefix:
        return
    for des, own, i in batch:
        suffix = '{0:0=6}'.format(des.frames[i])
        for name in names:
            target = name.format(own, suffix)
            if os.path.exists(target):   # a frame without markers has no .vtp
                link_file(target, name.format(prefix, suffix))


def write_series_file(series, file_of, filename):
    '''List the frames of the series whose file_of(des, i) exists in filename, the JSON file ParaView
    opens as one series: names relative to it, times each frame's time_yr (the .info time keeps
    7 digits).'''
    here = os.path.realpath(os.path.dirname(filename))
    files = []
    for des, i in series:
        fn = file_of(des, i)
        if os.path.exists(fn):
            files.append({'name': os.path.relpath(os.path.realpath(fn), here),
                          'time': des.read_field(des.frames[i], 'time_yr')[0]})
    with open(unlinked(filename), 'w') as f:
        json.dump({'file-series-version': '1.0', 'files': files}, f, indent=1)
    nrun = len({des.modelname for des, _ in series})
    print(f'Listed {len(files)} of the {len(series)} frames of {nrun} runs in {filename}.', file=sys.stderr)


def convert_batches(batches, prefix, target_func, imap, names):
    '''Convert the batches through imap, one model at a time, linking each parent's files
    of the given names as ours.'''
    nout = sum(len(b) for b in batches)
    width, ndone = len(str(nout)), 0
    for batch in batches:
        announce_model(batch, ndone, prefix, names)
        for result in imap(target_func, batch):
            ndone += 1
            print(f'Frame #{result} converted ({ndone:{width}d}/{nout}).', end='\r', file=sys.stderr)
        link_frames(batch, prefix, names)


def main(modelname, start, end, delta):
    series = restart_series(modelname, generate_parent_frames)
    binary = sorted({d.modelname for d, _ in series if d.format != 'hdf5'})
    if (update_vtkhdf or hdf_copies) and binary:
        print(f'Error: -u and -hdf update .vtkhdf frames; des-binary frames in {", ".join(binary)}.',
              file=sys.stderr)
        sys.exit(1)
    frames = [d.frames[i] for d, i in series]
    prefix = prefix_of(modelname)

    if start == -1:
        # after the last .vtu of these frames; a -g run may have left the parents' here too
        done = [int(f[len(prefix) + 1:-4]) for f in sorted(glob.glob(prefix + '.*.vtu'))]
        done = [n for n in done if n in frames]
        start = frames.index(done[-1]) + 1 if done else 0
    if end == -1:
        end = len(frames)

    selected = series[start:end:delta]
    # a parent's files go under its own name, to be linked or listed as ours, unless -copy
    # writes them under ours; -u updates every frame in its own file whatever the name
    def own_prefix(d):
        return prefix if d.modelname == modelname or not link_parent_frames else prefix_of(d.modelname)
    args_list = [(d, own_prefix(d), i) for d, i in selected]
    # one batch per model, in series order, so the model being worked on can be named
    batches = [list(b) for _, b in itertools.groupby(args_list, key=lambda args: args[0])]
    # a parent's .vtu/.vtp or -hdf copies are linked as it is converted, -u frames after the
    # update (below); under -copy they carry this model's name already
    if hdf_copies:
        link_names = VTKHDF_NAMES
    elif update_vtkhdf:
        link_names = None
    else:
        link_names = VTU_NAMES

    # frame numbers, not series positions: with -g the series spans several runs
    converted = [d.frames[i] for batch in batches for d, _, i in batch]
    nout = len(converted)

    if nout == 0:
        print(f'No frames to convert (Avail. frames: {frames[0]} to {frames[-1]}).', file=sys.stderr)
        return
    
    target_func = process_vtkhdf_update if update_vtkhdf or hdf_copies else process_single_frame
    try:
        import multiprocessing as mp
        print(f'Using {mutiprocessing_threads} threads (-ncpu {mutiprocessing_threads}) for conversion (system max: {mp.cpu_count()}).', file=sys.stderr)
        print(f'Converting {nout} frames from {converted[0]} to {converted[-1]} with step {delta}.',
              file=sys.stderr)

        options = {name: globals()[name] for name in WORKER_OPTIONS}
        with mp.Pool(processes = mutiprocessing_threads, initializer = init_worker,
                     initargs = (options,)) as pool:
            convert_batches(batches, prefix, target_func, pool.imap_unordered, link_names)

    except ImportError:
        print('Multiprocessing is not available, using single thread instead.')
        convert_batches(batches, prefix, target_func, map, link_names)
            
        
    print(file=sys.stderr)   # ends the last progress line
    # the .series lists the series' frames on disk: -u the .save frames, -hdf the copies
    if update_vtkhdf:
        if generate_parent_frames:   # every frame in the work directory as this model's .save frame
            nlinked = 0
            for d, i in series:
                target = d.get_fn(d.frames[i])
                link = '{0}.save.{1:0=6}.vtkhdf'.format(prefix, d.frames[i])
                if os.path.realpath(link) == os.path.realpath(target):   # there already
                    continue
                if os.path.exists(link) and not os.path.islink(link):   # a restart left its old frame
                    print(f'Warning: {link} is a frame, not a link; left as it is.', file=sys.stderr)
                    continue
                link_file(target, link)
                nlinked += 1
            print(f'Linked {nlinked} frames as {prefix}.save.NNNNNN.vtkhdf.', file=sys.stderr)
        write_series_file(series, lambda d, i: d.get_fn(d.frames[i]), prefix + '.save.vtkhdf.series')
    elif hdf_copies:
        write_series_file(series, lambda d, i: '{0}.{1:0=6}.vtkhdf'.format(own_prefix(d), d.frames[i]),
                          prefix + '.vtkhdf.series')
    return


def output_vtp_file(des, frame, filename, markersetname, time_in_yr, step):
    fvtp = open(unlinked(filename), 'w')

    class MarkerSizeError(RuntimeError):
        pass

    try:
        # read data
        marker_data = des.read_markers(frame, markersetname)
        nmarkers = marker_data['size']

        if nmarkers <= 0:
            raise MarkerSizeError()

        # write vtp header
        vtp_header(fvtp, nmarkers, time_in_yr, step)

        # point-based data
        fvtp.write('  <PointData>\n')
        name = markersetname + '.mattype'
        marker_type = marker_data[name]
        vtk_dataarray(fvtp, marker_type, name)
        name = markersetname + '.elem'
        vtk_dataarray(fvtp, marker_data[name], name)
        name = markersetname + '.id'
        vtk_dataarray(fvtp, marker_data[name], name)
        for name in (markersetname + '.time',markersetname + '.z', markersetname + '.distance',markersetname+'.slope'):    
            try:
                vtk_dataarray(fvtp, marker_data[name], name)
            except:
                pass
        fvtp.write('  </PointData>\n')

        # point coordinates
        fvtp.write('  <Points>\n')
        field = marker_data[markersetname + '.coord']
        if des.ndims == 2:
            # VTK requires vector field (velocity, coordinate) has 3 components.
            # Allocating a 3-vector tmp array for VTK data output.
            tmp = np.zeros((nmarkers, 3), dtype=field.dtype)
            tmp[:,:des.ndims] = field
        else:
            tmp = field

        vtk_dataarray(fvtp, tmp, markersetname + '.coord', 3)
        fvtp.write('  </Points>\n')

        vtp_footer(fvtp)
        fvtp.close()

    except MarkerSizeError:
        # delete partial vtp file
        fvtp.close()
        os.remove(filename)
        # skip this frame

    except:
        # delete partial vtp file
        fvtp.close()
        os.remove(filename)
        raise

    return


def convert_field(des, frame, name, fvtu):
    field = des.read_field(frame, name)
    if name in ('coordinate', 'velocity', 'velocity averaged', 'force'):
        if des.ndims == 2:
            # VTK requires vector field (velocity, coordinate) has 3 components.
            # Allocating a 3-vector tmp array for VTK data output.
            i = des.frames.index(frame)
            tmp = np.zeros((des.nnode_list[i], 3), dtype=field.dtype)
            tmp[:,:des.ndims] = field
        else:
            tmp = field

        # Rename 'velocity averaged' to 'velocity'
        if name == 'velocity averaged': name = 'velocity'

        vtk_dataarray(fvtu, tmp, name, 3)
    else:
        vtk_dataarray(fvtu, field, name)
    return


def vtk_dataarray(f, data, data_name=None, data_comps=None):
    if data.dtype in (np.int32, np.uint32):
        dtype = 'Int32'
    elif data.dtype in (np.single, np.float32):
        dtype = 'Float32'
    elif data.dtype in (np.double, np.float64):
        dtype = 'Float64'
    else:
        raise Error('Unknown data type: ' + name)

    name = ''
    if data_name:
        name = 'Name="{0}"'.format(data_name)

    ncomp = ''
    if data_comps:
        ncomp = 'NumberOfComponents="{0}"'.format(data_comps)

    if output_in_binary:
        fmt = 'binary'
    else:
        fmt = 'ascii'
    header = '<DataArray type="{0}" {1} {2} format="{3}">\n'.format(
        dtype, name, ncomp, fmt)
    f.write(header)
    if output_in_binary:
        header = np.zeros(4, dtype=np.int32)
        header[0] = 1
        a = data.tobytes()
        header[1] = len(a)
        header[2] = len(a)
        b = zlib.compress(a)
        header[3] = len(b)
        f.write(base64.standard_b64encode(header.tobytes()).decode('ascii'))
        f.write(base64.standard_b64encode(b).decode('ascii'))
    else:
        data.tofile(f, sep=' ')
    f.write('\n</DataArray>\n')
    return


def vtu_header(f, nnode, nelem, time, step):
    f.write(
'''<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">
<UnstructuredGrid>
<FieldData>
  <DataArray type="Float32" Name="TIME" NumberOfTuples="1" format="ascii">
    {2}
  </DataArray>
  <DataArray type="Float32" Name="CYCLE" NumberOfTuples="1" format="ascii">
    {3}
  </DataArray>
</FieldData>
<Piece NumberOfPoints="{0}" NumberOfCells="{1}">
'''.format(nnode, nelem, time, step))
    return


def vtu_footer(f):
    f.write(
'''</Piece>
</UnstructuredGrid>
</VTKFile>
''')
    return


def vtp_header(f, nmarkers, time, step):
    f.write(
'''<?xml version="1.0"?>
<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">
<PolyData>
<FieldData>
  <DataArray type="Float32" Name="TIME" NumberOfTuples="1" format="ascii">
    {1}
  </DataArray>
  <DataArray type="Float32" Name="CYCLE" NumberOfTuples="1" format="ascii">
    {2}
  </DataArray>
</FieldData>
<Piece NumberOfPoints="{0}">
'''.format(nmarkers, time, step))
    return


def vtp_footer(f):
    f.write(
'''</Piece>
</PolyData>
</VTKFile>
''')
    return


def first_invariant(t):
    nstr = t.shape[1]
    ndims = 2 if (nstr == 3) else 3
    return np.sum(t[:,:ndims], axis=1) / ndims


def second_invariant(t):
    '''The second invariant of the deviatoric part of a symmetric tensor t,
    where t[:,0:ndims] are the diagonal components;
      and t[:,ndims:] are the off-diagonal components.'''
    nstr = t.shape[1]

    # second invariant: sqrt(0.5 * t_ij**2)
    if nstr == 3:  # 2D
        return np.sqrt(0.25 * (t[:,0] - t[:,1])**2 + t[:,2]**2)
    else:  # 3D
        a = (t[:,0] + t[:,1] + t[:,2]) / 3
        return np.sqrt( 0.5 * ((t[:,0] - a)**2 + (t[:,1] - a)**2 + (t[:,2] - a)**2) +
                        t[:,3]**2 + t[:,4]**2 + t[:,5]**2)


def compute_principal_stress(stress):
    '''The principal stress (s1 and s3) of the deviatoric stress tensor.'''

    nelem = stress.shape[0]
    nstr = stress.shape[1]
    # VTK requires vector field (velocity, coordinate) has 3 components.
    # Allocating a 3-vector tmp array for VTK data output.
    s1 = np.zeros((nelem, 3), dtype=stress.dtype)
    s3 = np.zeros((nelem, 3), dtype=stress.dtype)

    if nstr == 3:  # 2D
        sxx, szz, sxz = stress[:,0], stress[:,1], stress[:,2]
        mag = np.sqrt(0.25*(sxx - szz)**2 + sxz**2)
        theta = 0.5 * np.arctan2(2*sxz, sxx-szz)
        cost = np.cos(theta)
        sint = np.sin(theta)

        s1[:,0] = mag * sint
        s1[:,1] = mag * cost
        s3[:,0] = mag * cost
        s3[:,1] = -mag * sint

    else:  # 3D
        # lower part of symmetric stress tensor
        s = np.zeros((nelem, 3,3), dtype=stress.dtype)
        s[:,0,0] = stress[:,0]
        s[:,1,1] = stress[:,1]
        s[:,2,2] = stress[:,2]
        s[:,1,0] = stress[:,3]
        s[:,2,0] = stress[:,4]
        s[:,2,1] = stress[:,5]

        # eigenvalues and eigenvectors
        if eigh_vectorized:
            # Numpy-1.8 or newer
            w, v = eigh(s)
        else:
            # Numpy-1.7 or older
            w = np.zeros((nelem,3), dtype=stress.dtype)
            v = np.zeros((nelem,3,3), dtype=stress.dtype)
            for e in range(nelem):
                w[e,:], v[e,:,:] = eigh(s[e])

        # isotropic part to be removed
        m = np.sum(w, axis=1) / 3

        p = w.argmin(axis=1)
        t = w.argmax(axis=1)
        #print(w.shape, v.shape, p.shape)

        for e in range(nelem):
            s1[e,:] = (w[e,p[e]] - m[e]) * v[e,:,p[e]]
            s3[e,:] = (w[e,t[e]] - m[e]) * v[e,:,t[e]]

    return s1, s3

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def is_fraction(s):
    try:
        Fraction(s)
        return True
    except (ValueError, ZeroDivisionError):
        return False

def parse_thread_number(value):
    try:
        value = float(value)
        if (value < 1 and value > 0):
            value = int(mp.cpu_count() * value)
        else:
            value = min(mp.cpu_count(), int(value))
        value = int(value)
    except ValueError:
        if value in ('True', 'true', 'all', 'yes'):
            value = mp.cpu_count()
        elif is_fraction(value):
            value = int(Fraction(value) * mp.cpu_count())
        else:
            raise ValueError(f'Invalid value for -ncpu: {value}')

    return min(mp.cpu_count(), value)

def read_and_remove_option_value(option):
    idx = sys.argv.index(option)
    try:
        value = sys.argv[idx+1]
    except IndexError:
        raise ValueError(f'Option {option} requires a value, but no value was provided.')
    if value[0] == '-' and not is_number(value[1:]):
        raise ValueError(f'Option {option} requires a value, but got {value}.')
        
    del sys.argv[idx+1:idx+2]
    return value


if __name__ == '__main__':

    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    else:
        for arg in sys.argv[1:]:
            if arg.lower() in ('-h', '--help'):
                print(__doc__)
                sys.exit(0)

    if '-a' in sys.argv:
        output_in_binary = False
    if '-c' in sys.argv:
        output_in_cwd = True
    if '-p' in sys.argv:
        output_principle_stress = True
    if '-t' in sys.argv:
        output_tensor_components = True
    if '-m' in sys.argv:
        output_markers = True
    if '-u' in sys.argv or '--update-vtkhdf' in sys.argv:
        update_vtkhdf = True
    if '-hdf' in sys.argv:
        hdf_copies = True
    if '-copy' in sys.argv:
        link_parent_frames = False
    if '-g' in sys.argv:
        generate_parent_frames = True
    for bad, why in ((h5py is None and (update_vtkhdf or hdf_copies),
                      '-u and -hdf need h5py, which is not installed'),
                     (update_vtkhdf and hdf_copies,
                      '-u updates the .save frames in place and -hdf writes updated copies: use one'),
                     (not link_parent_frames and not generate_parent_frames,
                      '-copy copies the parent frames that -g adds'),
                     (not link_parent_frames and update_vtkhdf,
                      '-u updates every frame in place, so -copy has nothing to copy')):
        if bad:
            print(f'Error: {why}.', file=sys.stderr)
            sys.exit(1)
    if '-melt' in sys.argv:
        output_melting = True
    if '-heat' in sys.argv:
        output_heatflux = False
    if '-ncpu' in sys.argv:
        try:
            import multiprocessing as mp
            value = read_and_remove_option_value('-ncpu')
            mutiprocessing_threads = parse_thread_number(value)
        except ImportError:
            raise ImportError('Multiprocessing is not available, please install it to use -ncpu option.')

    # delete options
    narvg = len(sys.argv)
    idx = 0
    while (idx < narvg):
        if sys.argv[idx][0] == '-' and not is_number(sys.argv[idx][1:]):
            del sys.argv[idx]
            narvg -= 1
        else:
            idx += 1

    modelname = sys.argv[1]

    if len(sys.argv) < 3:
        start = 0
    else:
        start = int(sys.argv[2])

    if len(sys.argv) < 4 or int(sys.argv[3]) == -1:
        end = -1
    else:
        end = int(sys.argv[3]) + 1

    if len(sys.argv) < 5:
        delta = 1
    else:
        delta = int(sys.argv[4])

    main(modelname, start, end, delta)
