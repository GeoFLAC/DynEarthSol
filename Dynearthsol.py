#!/usr/bin/env python

from __future__ import print_function, unicode_literals
import sys, os
import numpy as np

class Dynearthsol:
    '''Read output file of 2D/3D DynEarthSol'''

    def __init__(self, modelname):
        self.suffix = 'save'
        self.modelname = modelname
        self.format = ''
        self.read_info()
        self.read_header(self.frames[0])
        return


    def read_info(self):
        tmp = np.fromfile(self.modelname + '.info', dtype=float, sep=' ')
        tmp.shape = (-1, 8)
        self.frames = list(tmp[:,0].astype(int))
        self.steps = list(tmp[:,1].astype(int))
        self.time = list(tmp[:,2].astype(float))
        self.nnode_list = tmp[:,5].astype(int)
        self.nelem_list = tmp[:,6].astype(int)
        return


    def get_fn(self, frame):
        filename = '{0}.{1}.{2:0=6}'.format(self.modelname, self.suffix, frame)
        if self.format == '':
            if os.path.isfile(filename):
                self.format = 'binary'
            else:
                filename += '.nc'
                self.format = 'netcdf'
            return filename
        elif self.format == 'binary':
            return filename
        elif self.format == 'netcdf':
            return filename + '.nc'

    def read_header(self, frame):
        self._header_frame = frame
        headerlen = 4096
        fname = self.get_fn(frame)
        self.field_pos = {}

        if self.format == 'netcdf':
            import netCDF4
            with netCDF4.Dataset(fname, 'r') as f:
                self.ndims = f.getncattr('ndims')
                self.revision = f.getncattr('revision')
                
                if self.ndims == 2:
                    self.nstr = 3
                    self.component_names = ('XX', 'ZZ', 'XZ')
                else:
                    self.nstr = 6
                    self.component_names = ('XX', 'YY', 'ZZ', 'XY', 'XZ', 'YZ')
        else:
            with open(fname, 'rb') as f:
                header = f.read(headerlen).splitlines()
                #print(header)

            # parsing 1st line
            first = header[0].split(b' ')
            if (first[0] != b'#' or
                first[1] != b'DynEarthSol' or
                first[2].split(b'=')[0] != b'ndims' or
                first[3].split(b'=')[0] != b'revision'):
                print('Error:', fname, 'is not a valid DynEarthSol output file!')
                sys.exit(1)

            self.ndims = int(first[2].split(b'=')[1])
            self.revision = int(first[3].split(b'=')[1])
            if self.ndims == 2:
                self.nstr = 3
                self.component_names = ('XX', 'ZZ', 'XZ')
            else:
                self.nstr = 6
                self.component_names = ('XX', 'YY', 'ZZ', 'XY', 'XZ', 'YZ')

            # parsing other lines
            for line in header[1:]:
                # test for null in python3 bytes and python2 str
                if line[0] in (0, '\x00'): break  # end of record
                name, pos = line.split(b'\t')
                self.field_pos[name.decode('ascii')] = int(pos)

            #print(self.field_pos)
        return


    def _get_dtype_count_shape(self, frame, name):
        i = self.frames.index(frame)
        nnode = self.nnode_list[i]
        nelem = self.nelem_list[i]

        dtype = np.float64 if name not in ('connectivity', 'bcflag') else np.int32

        if name in set(['strain', 'strain-rate', 'stress', 'stress averaged']):
            count = self.nstr * nelem
            shape = (nelem, self.nstr)
        elif name in set(['density', 'material', 'mesh quality',
                          'plastic strain', 'plastic strain-rate',
                          'viscosity', 'edvoldt', 'volume','radiogenic source']):
            count = nelem
            shape = (nelem, )
        elif name in set(['connectivity']):
            count = (self.ndims + 1) * nelem
            shape = (nelem, self.ndims+1)
        elif name in set(['coordinate', 'velocity', 'velocity averaged', 'force', 'coord0']):
            count = self.ndims * nnode
            shape = (nnode, self.ndims)
        elif name in set(['bcflag', 'temperature', 'pore pressure', 'mass', 'tmass', 'volume_n']):
            count = nnode
            shape = (nnode, )
        else:
            raise NameError('uknown field name: ' + name)
        return dtype, count, shape


    def read_field(self, frame, name):
        if frame != self._header_frame: self.read_header(frame)
        fname = self.get_fn(frame)

        if self.format == 'netcdf':
            import netCDF4
            with netCDF4.Dataset(fname, 'r') as f:
                field = f.variables[name][:]
        else:
            dtype, count, shape = self._get_dtype_count_shape(frame, name)

            pos = self.field_pos[name]
            with open(fname,'r') as f:
                f.seek(pos)
                field = np.fromfile(f, dtype=dtype, count=count).reshape(shape)
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
        


    def overwrite_field(self, frame, name, data):
        if frame != self._header_frame: self.read_header(frame)
        if frame != self._header_frame: read_header(frame)
        if name.startswith(('markerset.', 'hydrous-markerset.')):
            dtype = data.dtype
            count = len(data)
            shape = (count,)
        else:
            dtype, count, shape = self._get_dtype_count_shape(frame, name)
        if data.shape != shape:
            raise Error('Shape of {0} field is changed! Expecting {1}, got {2}.'.format(name, shape, data.shape))

        pos = self.field_pos[name]
        fname = self.get_fn(frame)
        with open(fname, 'r+') as f:
            f.seek(pos)
            f.write(data.tostring())
        return


    def read_markers(self, frame, markername):
        'Read and return marker data'
        if frame != self._header_frame: self.read_header(frame)
        fname = self.get_fn(frame)
        if self.format == 'netcdf':
            import netCDF4
            with netCDF4.Dataset(fname, 'r') as f:
                marker_data = {}
                for var in f.variables:
                    if var.startswith(markername):
                        marker_data[var] = f.variables[var][:]
                marker_data['size'] = marker_data[markername+'.elem'].shape[0]
        else:
            with open(fname) as f:

                pos = self.field_pos[markername+' size']
                f.seek(pos)
                nmarkers = np.fromfile(f, dtype=np.int32, count=1)[0]

                marker_data = {'size': nmarkers}

                # floating point
                for name in (markername+'.coord',):
                    pos = self.field_pos[name]
                    f.seek(pos)
                    tmp = np.fromfile(f, dtype=np.float64, count=nmarkers*self.ndims)
                    marker_data[name] = tmp.reshape(-1, self.ndims)
                    #print(marker_data[name].shape, marker_data[name])

                try:
                    for name in (markername+'.eta',):
                        pos = self.field_pos[name]
                        f.seek(pos)
                        tmp = np.fromfile(f, dtype=np.float64, count=nmarkers*(self.ndims+1))
                        marker_data[name] = tmp.reshape(-1, (self.ndims+1))
                except:
                    pass

                # int
                for name in (markername+'.elem', markername+'.mattype', markername+'.id'):
                    pos = self.field_pos[name]
                    f.seek(pos)
                    marker_data[name] = np.fromfile(f, dtype=np.int32, count=nmarkers)
                    #print(marker_data[name].shape, marker_data[name])

                # float
                for name in (markername+'.time',markername+'.z',markername+'.distance',markername+'.slope'):
                    try:
                        pos = self.field_pos[name]
                        f.seek(pos)
                        marker_data[name] = np.fromfile(f, dtype=np.float64, count=nmarkers)
                    except:
                        pass

        return marker_data




class DynearthsolCheckpoint(Dynearthsol):
    '''Read chkpt file of 2D/3D DynEarthSol'''

    def __init__(self, modelname, frame):
        self.suffix = 'chkpt'
        self.modelname = modelname
        self.read_info()
        self.read_header(frame)
        return


    def _get_dtype_count_shape(self, frame, name):
        i = self.frames.index(frame)
        nnode = self.nnode_list[i]
        nelem = self.nelem_list[i]

        dtype = np.float64 if name != 'connectivity' else np.int32

        if name in set(['volume_old']):
            count = nelem
            shape = (nelem, )
        else:
            raise NameError('uknown field name: ' + name)
        return dtype, count, shape


    def read_markers(self, frame, markername):
        'Read and return marker data'
        if frame != self._header_frame: read_header(frame)
        fname = self.get_fn(frame)
        with open(fname) as f:

            pos = self.field_pos[markername+' size']
            f.seek(pos)
            nmarkers = np.fromfile(f, dtype=np.int32, count=1)[0]

            marker_data = {'size': nmarkers}

            # floating point
            for name in (markername+'.eta',):
                pos = self.field_pos[name]
                f.seek(pos)
                tmp = np.fromfile(f, dtype=np.float64, count=nmarkers*(self.ndims+1))
                marker_data[name] = tmp.reshape(-1, self.ndims+1)
                #print(marker_data[name].shape, marker_data[name])

            # int
            for name in (markername+'.elem', markername+'.mattype', markername+'.id'):
                pos = self.field_pos[name]
                f.seek(pos)
                marker_data[name] = np.fromfile(f, dtype=np.int32, count=nmarkers)
                #print(marker_data[name].shape, marker_data[name])

            # float
            for name in (markername+'.time',markername+'.z',markername+'.distance',markername+'.slope'):
                try:
                    pos = self.field_pos[name]
                    f.seek(pos)
                    marker_data[name] = np.fromfile(f, dtype=np.float64, count=nmarkers)
                except:
                    pass

        return marker_data
