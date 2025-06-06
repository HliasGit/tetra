import os
import subprocess
import time
import shutil

mol = "1e7p"

results_dir = '/home/fs72740/evaglietti/tetra/marching_tetrahedron/script/results'
mol_dir = os.path.join(results_dir, mol)

os.makedirs(results_dir, exist_ok=True)

if os.path.exists(mol_dir):
    shutil.rmtree(mol_dir)
os.makedirs(mol_dir)

os.chdir(mol_dir)

modality = [{'mod': "cpu_unique",
            'name': "test_cpu/main",
            'args': ['4e-4', mol, 'midpoint']},
            {'mod': "cpu_nonunique",
            'name': "test_cpu/non_unique",
            'args': ['4e-4', mol]},
            {'mod': "gpu_preprocessing",
            'name': "test_gpu/test",
            'args': []},
            {'mod': "gpu_noprep",
            'name': "test_gpu/test_wo_preprocessing",
            'args': []}]

for m in modality:
    with open(m['mod'] + ".txt", 'w') as f:
        for i in range(30):
            start_time = time.time()
            subprocess.run(['../../../build/' + m['name']] + m['args'])
            end_time = time.time()
            f.write(f'{end_time - start_time}\n')