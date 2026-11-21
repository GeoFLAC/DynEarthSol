"""Check scheduler guards without launching MPI or contacting a scheduler."""
import hashlib
import json
import os
from pathlib import Path
import socket
import subprocess
import sys
import tempfile
import unittest


LAUNCHER = Path(__file__).with_name('launch_sge.py').resolve()


class SchedulerGuards(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='des-sge-policy-')
        self.addCleanup(self.temp.cleanup)
        self.directory = Path(self.temp.name)
        self.hostfile = self.directory / 'pe-hostfile'
        self.host = socket.gethostname().split('.')[0].lower()
        self.hostfile.write_text(self.host + ' 2 all.q UNDEFINED\n')
        self.env = dict(os.environ, JOB_ID='123', NSLOTS='2',
                        PE_HOSTFILE=str(self.hostfile), OMP_NUM_THREADS='1')

    def reject(self, arguments, expected, env=None):
        result = subprocess.run([sys.executable, str(LAUNCHER)] + arguments,
                                cwd=str(self.directory), env=self.env if env is None else env,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                universal_newlines=True, timeout=10)
        self.assertNotEqual(result.returncode, 0, result.stdout)
        self.assertIn(expected, result.stdout)

    def test_requires_qsub_context(self):
        env = {key:value for key, value in self.env.items()
               if key not in ('JOB_ID', 'NSLOTS', 'PE_HOSTFILE')}
        self.reject(['-n', '2', sys.executable], 'Submit through qsub', env)

    def test_rejects_inconsistent_slot_accounting(self):
        self.env['NSLOTS'] = '8'
        self.reject(['-n', '2', sys.executable], 'slot total does not match NSLOTS')

    def test_rejects_hybrid_oversubscription_before_mpi(self):
        self.env['OMP_NUM_THREADS'] = '2'
        self.reject(['--mpiexec', '/missing-mpiexec', '-n', '2', sys.executable],
                    'exceed the allocated slot capacity')

    def test_rejects_unallocated_node_count(self):
        self.reject(['--min-nodes', '2', '-n', '2', sys.executable],
                    'too few distinct hosts')

    def worker(self, context, expected):
        path = self.directory/'allocation.json'
        path.write_text(json.dumps(context))
        self.env.update(OMPI_COMM_WORLD_RANK='0', OMPI_COMM_WORLD_LOCAL_RANK='0',
                        OMPI_COMM_WORLD_LOCAL_SIZE='1')
        self.reject(['--rank', str(path), sys.executable], expected)

    def test_rejects_worker_on_unallocated_host(self):
        self.worker(dict(job_id='123', host_slots={'unallocated-host.invalid':2}),
                    'outside the scheduler allocation')

    def test_rejects_worker_with_another_job_id(self):
        self.worker(dict(job_id='456', host_slots={self.host:2}), 'job ID differs')

    def test_rejects_changed_worker_executable(self):
        self.worker(dict(job_id='123', host_slots={self.host:2}, threads=1,
                         executable_sha256='invalid'), 'executable differs')

    def test_rejects_per_host_oversubscription(self):
        self.worker(dict(job_id='123', host_slots={self.host:1}, threads=2,
                         executable_sha256=hashlib.sha256(Path(sys.executable).read_bytes()).hexdigest()),
                    'exceed allocated slots')


if __name__ == '__main__':
    unittest.main()
