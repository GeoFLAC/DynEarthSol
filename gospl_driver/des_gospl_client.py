"""Root-side adapter for the existing C++ bridge; imports no MPI or GoSPL.

The launcher creates a private local socket and token. Requests contain only
coupling boundary data; the model stays in its independent worker MPI world.
"""
import atexit
import os
from multiprocessing.connection import Client

_connection = None
_models = {}


def _rpc(operation, *args, **kwargs):
    global _connection
    if _connection is None:
        _connection = Client(os.environ['DES_GOSPL_ENDPOINT'], family='AF_UNIX',
                             authkey=bytes.fromhex(os.environ['DES_GOSPL_TOKEN']))
    _connection.send((operation, args, kwargs))
    timeout = float(os.environ.get('DES_GOSPL_TIMEOUT', '300'))
    if not _connection.poll(timeout):
        raise TimeoutError('GoSPL worker did not respond within %s seconds' % timeout)
    result = _connection.recv()
    if not result['ok']:
        raise RuntimeError(result['error'])
    return result['value']


class _Model:
    def __init__(self, handle):
        self.handle = handle

    @property
    def tNow(self):
        return get_current_time(self.handle)

    @property
    def verbose(self):
        return False

    @verbose.setter
    def verbose(self, value):
        _rpc('verbose', self.handle, bool(value))


def create_enhanced_model(config_path):
    handle = _rpc('create_enhanced_model', config_path)
    if handle >= 0:
        _models[handle] = _Model(handle)
    return handle


def destroy_model(handle):
    result = _rpc('destroy_model', handle)
    if result == 0:
        _models.pop(handle, None)
    return result


def _forward(name):
    def call(*args, **kwargs):
        return _rpc(name, *args, **kwargs)
    call.__name__ = name
    return call


for _name in ('run_processes_for_dt', 'run_processes_for_steps', 'run_processes_until_time',
              'apply_velocity_data', 'apply_elevation_data', 'interpolate_elevation_to_points',
              'get_current_time', 'get_time_step', 'set_surface_velocity', 'set_uplift_rate',
              'run_and_get_erosion', 'apply_drift_correction','save_checkpoint','load_checkpoint','prepare_restart'):
    globals()[_name] = _forward(_name)


@atexit.register
def _close():
    if _connection is not None:
        _connection.close()
