#!/usr/bin/env python3
"""Persistent real GoSPL model(s) in a separate MPI world, controlled by DES root."""
import argparse
import importlib.util
import json
import os
import sys
import socket
from pathlib import Path
import time
import traceback
from multiprocessing.connection import Listener
from mpi4py import MPI
import numpy as np

OPERATIONS = {'create_enhanced_model', 'destroy_model', 'run_processes_for_dt',
              'run_processes_for_steps', 'run_processes_until_time', 'apply_velocity_data',
              'apply_elevation_data', 'interpolate_elevation_to_points', 'get_current_time',
              'get_time_step', 'set_surface_velocity', 'set_uplift_rate',
              'run_and_get_erosion', 'apply_drift_correction', 'verbose',
              'save_checkpoint','load_checkpoint','prepare_restart'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--endpoint', required=True)
    parser.add_argument('--interface', type=Path, required=True)
    parser.add_argument('--audit-dir', type=Path)
    args = parser.parse_args()
    comm = MPI.COMM_WORLD
    rank = comm.rank
    hosts=comm.allgather(socket.gethostname())
    # Instrumented interfaces may write per-call arrays; give each rank its own
    # directory so no rank duplicates another rank's files.
    if args.audit_dir:
        folder = args.audit_dir / ('rank-%04d' % rank)
        folder.mkdir(parents=True, exist_ok=True)
        os.environ['GOSPL_AUDIT_DIR'] = str(folder.resolve())
    spec = importlib.util.spec_from_file_location('_des_real_gospl_interface', args.interface.resolve())
    interface = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(interface)
    listener = connection = log = None
    if rank == 0:
        listener = Listener(args.endpoint, family='AF_UNIX',
                            authkey=bytes.fromhex(os.environ['DES_GOSPL_TOKEN']))
        Path(args.endpoint + '.ready').write_text(json.dumps(dict(pid=os.getpid(), mpi_size=comm.size, hosts=hosts)))
        log = Path(args.endpoint + '.jsonl').open('a', buffering=1)
        connection = listener.accept()
    exit_status = 0
    try:
        while True:
            request = None
            if rank == 0:
                try:
                    request = connection.recv()
                except (EOFError, OSError):
                    request = ('stop', (), {})
            operation, values, keywords = comm.bcast(request, root=0)
            if operation == 'stop':
                break
            started = time.perf_counter()
            value = None
            error = None
            try:
                if operation not in OPERATIONS:
                    raise ValueError('Unsupported coupling operation: ' + operation)
                if operation in ('prepare_restart','save_checkpoint','load_checkpoint'):
                    import checkpoint
                    if operation == 'prepare_restart':
                        value=checkpoint.prepare(comm,*values)
                    else:
                        action=checkpoint.save if operation=='save_checkpoint' else checkpoint.load
                        value=action(interface._models[values[0]],comm,*values[1:])
                elif operation == 'verbose':
                    interface._models[values[0]].verbose = values[1]
                    value = 0
                else:
                    value = getattr(interface, operation)(*values, **keywords)
                    if operation == 'create_enhanced_model' and value >= 0:
                        # The installed bridge assumes a sequential hGlobal
                        # array. Adapt only coupling transfers to owned IDs.
                        from parallel_surface import bind_parallel_surface
                        bind_parallel_surface(interface._models[value], comm)
                if value is None or (isinstance(value, (int, float)) and value < 0):
                    raise RuntimeError('GoSPL bridge returned failure for ' + operation)
                if operation == 'run_and_get_erosion':
                    array = np.asarray(value)
                    if array.shape != (int(values[3]),) or not np.isfinite(array).all():
                        raise RuntimeError('Invalid erosion array returned by GoSPL')
            except Exception:
                error = traceback.format_exc()
            errors = comm.allgather(error)
            error = next((e for e in errors if e), None)
            if rank == 0:
                duration = time.perf_counter() - started
                log.write(json.dumps(dict(operation=operation, mpi_size=comm.size,
                                          seconds=duration, ok=error is None)) + '\n')
                connection.send(dict(ok=error is None, value=value if error is None else None,
                                     error=error))
            if error:
                # All ranks have agreed on the error. Destroy collectively and
                # leave through MPI_Finalize instead of concurrent MPI_Abort.
                exit_status = 53
                break
    finally:
        # Every rank destroys exactly its own distributed part of each model.
        for handle in list(interface._models):
            interface.destroy_model(handle)
        if rank == 0:
            connection.close()
            listener.close()
            log.close()
            Path(args.endpoint + '.ready').unlink(missing_ok=True)
    return exit_status


if __name__ == '__main__':
    try:
        status = main()
    except (Exception, KeyboardInterrupt, SystemExit):
        # An uncaught failure or interrupt on one rank must abort its peers;
        # otherwise they can remain blocked in a collective or finalization.
        traceback.print_exc()
        MPI.COMM_WORLD.Abort(53)
    else:
        sys.exit(status)
