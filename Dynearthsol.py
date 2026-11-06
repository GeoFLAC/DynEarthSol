#!/usr/bin/env python

from __future__ import print_function, unicode_literals
import sys, os
import numpy as np

# must match constants.hpp
YEAR2SEC = 365.2422 * 86400


def _scan_frame(fname):
    import h5py
    with h5py.File(fname, 'r') as f:
        # HDF5Output links every dataset to the root under its plain name, so
        # the .info columns are read the same way as any field.
        try:
            return dict(
                steps=int(f['steps'][0]),
                time=float(f['time_sec'][0]),
                dt=float(f['dt_sec'][0]),
                walltime=float(f['walltime_sec'][0]),
                nnode=int(f['nnode'][0]),
                nelem=int(f['nelem'][0]),
                nseg=int(f['nseg'][0]),
            )
        except KeyError as e:
            raise KeyError('%s missing in %s' % (e, fname))


def scan_frames(modelname):
    '''Read the per-frame .info metadata (frame, steps, time, dt, walltime,
    nnode, nelem, nseg) directly from <modelname>.save.NNNNNN.vtkhdf frame
    files, sorted by frame number.'''
    import glob, re
    pat = re.compile(re.escape(os.path.basename(modelname))
                     + r'\.save\.(\d+)\.vtkhdf$')
    found = []
    for fname in glob.glob(modelname + '.save.*'):
        m = pat.search(os.path.basename(fname))
        if m:
            found.append((int(m.group(1)), fname))
    found.sort()

    rows = []
    for frame, fname in found:
        d = _scan_frame(fname)
        d['frame'] = frame
        rows.append(d)
    return rows


class Dynearthsol:
    '''Read output file of 2D/3D DynEarthSol'''

    def __init__(self, modelname):
        self.suffix = 'save'
        self.modelname = modelname
        self.read_info()
        self.read_header(self.frames[0])
        return


    def read_info(self):
        # The .info index is one small file, so it is read when present; the
        # metadata embedded in every .save frame is scanned only without it.
        info = self.modelname + '.info'
        if os.path.isfile(info):
            tmp = np.fromfile(info, dtype=float, sep=' ')
            tmp.shape = (-1, 8)
            self.frames = list(tmp[:,0].astype(int))
            self.steps = list(tmp[:,1].astype(int))
            self.time = list(tmp[:,2].astype(float))
            self.nnode_list = tmp[:,5].astype(int)
            self.nelem_list = tmp[:,6].astype(int)
            return
        print('Warning: cannot open info file {0}; using metadata embedded in '
              'the {1}.save.* frames.'.format(info, self.modelname), file=sys.stderr)
        try:
            rows = scan_frames(self.modelname)
        except KeyError as e:
            raise IOError('{0}; no {1} fallback found'.format(e, info))
        if not rows:
            raise IOError('neither {0}.save.* frame files nor {0}.info '
                          'found'.format(self.modelname))
        self.frames = [r['frame'] for r in rows]
        self.steps = [r['steps'] for r in rows]
        self.time = [r['time'] for r in rows]
        self.nnode_list = np.array([r['nnode'] for r in rows], dtype=int)
        self.nelem_list = np.array([r['nelem'] for r in rows], dtype=int)


    def get_fn(self, frame):
        return '{0}.{1}.{2:0=6}.vtkhdf'.format(self.modelname, self.suffix, frame)

    def read_header(self, frame):
        self._header_frame = frame
        import h5py
        with h5py.File(self.get_fn(frame), 'r') as f:
            self.ndims = f.attrs['ndims']
            self.revision = f.attrs['revision']

            if self.ndims == 2:
                self.nstr = 3
                self.component_names = ('XX', 'ZZ', 'XZ')
            else:
                self.nstr = 6
                self.component_names = ('XX', 'YY', 'ZZ', 'XY', 'XZ', 'YZ')
        return


    def read_field(self, frame, name):
        if frame != self._header_frame: self.read_header(frame)
        import h5py
        with h5py.File(self.get_fn(frame), 'r') as f:
            field = f[name][:]
        return field
    
    def load_calculation(self, frame, name):
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

        if name == 'heat flux':
            conductivity = 3.3
            coord = self.read_field(frame, 'coordinate')
            temperature = self.read_field(frame, 'temperature')
            connectivity = self.read_field(frame, 'connectivity')
            
            p = np.transpose(coord[connectivity], (1,2,0))
            t = np.transpose(temperature[connectivity], (1,0))

            v_arr = np.zeros((self.ndims,self.ndims+1,connectivity.shape[0]))                
            v_arr[:,:-1] = p[:-1] - p[-1]
            v_arr[:,-1] = t[:-1] - t[-1]

            nv = np.cross(v_arr[0].T,v_arr[1].T)
            v_slope = np.cross(nv, np.cross(nv, [0,0,1])).T
            
            norm = np.linalg.norm(v_slope[:-1], axis=0)
            
            flux_val = -1.e3 * conductivity * v_slope[-1] / norm

            flux = v_slope[:-1] * flux_val / norm
            return flux, flux_val
        
        elif name == 'strain II':
            strain = self.read_field(frame, 'strain')
            sII = second_invariant(strain)
            return sII
        elif name == 'ecoordinate':
            coord = self.read_field(frame, 'coordinate')
            connectivity = self.read_field(frame, 'connectivity')
            ecoord = coord[connectivity]
            return ecoord
        else:
            raise NameError('uknown field name: ' + name)
        


    def marker_sets(self, frame=None):
        '''Names of the marker sets in this frame; "hydrous-markerset" is there
        only with control.has_hydration_processes.'''
        if frame is None: frame = self._header_frame
        if frame != self._header_frame: self.read_header(frame)
        import h5py
        with h5py.File(self.get_fn(frame), 'r') as f:
            return sorted(k[:-len('.nmarkers')] for k in f.keys()
                          if k.endswith('.nmarkers'))


    def read_markers(self, frame, markername):
        'Read and return marker data'
        if frame != self._header_frame: self.read_header(frame)
        import h5py
        with h5py.File(self.get_fn(frame), 'r') as f:
            marker_data = {}
            for var in f.keys():
                if var.startswith(markername):
                    marker_data[var] = f[var][:]
            marker_data['size'] = marker_data[markername+'.elem'].shape[0]
        return marker_data




class DynearthsolCheckpoint(Dynearthsol):
    '''Read chkpt file of 2D/3D DynEarthSol'''

    def __init__(self, modelname, frame):
        self.suffix = 'chkpt'
        self.modelname = modelname
        self.read_info()
        self.read_header(frame)
        return


    def read_markers(self, frame, markername):
        '''Read and return the marker data a checkpoint holds: last_id, reserved_space and
        genesis. Positions, element, mattype and time are in the save frame
        (Dynearthsol.read_markers).'''
        if frame != self._header_frame: self.read_header(frame)
        import h5py
        with h5py.File(self.get_fn(frame), 'r') as f:
            last_id = f[markername+'.last_id'][0]
            reserved_space = f[markername+'.reserved_space'][0]
            genesis = f[markername+'.genesis'][:]
        return {'size': len(genesis), markername+'.last_id': int(last_id),
                markername+'.reserved_space': int(reserved_space), markername+'.genesis': genesis}
