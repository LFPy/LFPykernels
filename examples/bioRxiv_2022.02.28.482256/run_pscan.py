import os
import subprocess as sp
import json
import hashlib
from parameters import ParameterSpace, ParameterSet, ParameterRange
import example_network_parameters as params

for dir in ['jobs', 'parameters', 'logs', 'output']:
    if not os.path.isdir(dir):
        os.mkdir(dir)

# Survey connection weights (max syn. conductance) parameterspace
PS0 = ParameterSpace(
    ParameterSet(
        dict(
            weight_EE=ParameterRange([0.00015]),  # E onto E
            weight_IE=ParameterRange([0.000125]),  # E onto I
            weight_EI=ParameterRange([0.0045]),  # I onto E
            weight_II=ParameterRange([0.0020]),  # I onto I
            # linear scaling of all connection weights
            weight_scaling=ParameterRange([0.25, 0.5, 1., 2., 4.]),
            n_ext=ParameterRange([[465, 160]]),
        )
    )
)

PS0.save('PS0.txt')

# Parameterspace for reconstructed signals generation
PS1 = ParameterSpace(
    ParameterSet(
        dict(
            weight_EE=ParameterRange([0.00015]),
            weight_IE=ParameterRange([0.000125]),
            weight_EI=ParameterRange([0.0045]),
            weight_II=ParameterRange([0.0020]),
            weight_scaling=ParameterRange([1.]),
            biophys=ParameterRange(['lin',
                                    'frozen']),
            i_syn=ParameterRange([True]),
            n_ext=ParameterRange([[465, 160]]),
            g_eff=ParameterRange([False, True])
        )
    )
)

PS1.save('PS1.txt')

# Parameterspace for reconstructed signals generation
PS2 = ParameterSpace(
    ParameterSet(
        dict(
            weight_EE=ParameterRange([0.00015]),
            weight_IE=ParameterRange([0.000125]),
            weight_EI=ParameterRange([0.0045]),
            weight_II=ParameterRange([0.0020]),
            weight_scaling=ParameterRange([1.]),
            biophys=ParameterRange(['lin',
                                    'frozen']),
            n_ext=ParameterRange([[465, 160]]),
            t_E=ParameterRange([200.]),
            t_I=ParameterRange([400.]),
            g_eff=ParameterRange([True])
        )
    )
)

PS2.save('PS2.txt')


#######
# PS0 jobs
#######

# slurm job settings (shared)
if 'HOSTNAME' in os.environ.keys():
    if 'BUDGET_ACCOUNTS' in os.environ.keys():
        ACCOUNT = os.environ['BUDGET_ACCOUNTS']
    if os.environ['HOSTNAME'].rfind('jusuf') >= 0:
        PARTITION = 'batch'
    elif os.environ['HOSTNAME'].rfind('jr') >= 0:
        PARTITION = 'dc-cpu'
else:
    ACCOUNT = None
    PARTITION = None


# estimate run times
# (s) buffer + setup + create + connect + simulation + save
wt = 120 + 2 + 10 + 180 + 240 / 6000 * params.networkParameters['tstop'] + 180
TIME = '%i:%.2i:%.2i' % (wt // 3600,
                         (wt - wt // 3600 * 3600) // 60,
                         (wt - wt // 60 * 60))
LNODES = 8
NTASKS = 1024

job = """#!/bin/bash
##################################################################
#SBATCH --account {}
#SBATCH --partition {}
#SBATCH --job-name {}
#SBATCH --time {}
#SBATCH -o logs/{}_stdout.txt
#SBATCH -e logs/{}_error.txt
#SBATCH -N {}
#SBATCH --ntasks {}
##################################################################
# from here on we can run whatever command we want
unset DISPLAY # DISPLAY somehow problematic with Slurm
srun --mpi=pmi2 python -u example_network.py {}
"""

if 'HOSTNAME' in os.environ.keys():
    if os.environ['HOSTNAME'].rfind('jusuf') >= 0 or \
            os.environ['HOSTNAME'].rfind('jr') >= 0:
        # container for job IDs
        jobIDs = []
        for pset in PS0.iter_inner():
            # sorted json dictionary
            js = json.dumps(pset, sort_keys=True).encode()
            md5 = hashlib.md5(js).hexdigest()

            # save parameter file
            pset.save(url=os.path.join('parameters', '{}.txt'.format(md5)))

            # create job script
            with open(os.path.join('jobs', '{}.job'.format(md5)), 'w') as f:
                f.writelines(job.format(
                    ACCOUNT,
                    PARTITION,
                    md5,
                    TIME,
                    md5,
                    md5,
                    LNODES,
                    NTASKS,
                    md5
                ))
            cmd = ' '.join(['sbatch',
                            '{}'.format(os.path.join('jobs',
                                        '{}.job'.format(md5)))])
            print(cmd)
            output = sp.getoutput(cmd)
            jobid = output.split(' ')[-1]
            jobIDs.append((md5, jobid))
    else:
        raise NotImplementedError(
            f"do not recognize HOST {os.environ['HOSTNAME']}")
else:
    for pset in PS0.iter_inner():
        # sorted json dictionary
        js = json.dumps(pset, sort_keys=True).encode()
        md5 = hashlib.md5(js).hexdigest()

        # save parameter file
        pset.save(url=os.path.join('parameters', '{}.txt'.format(md5)))

        # run model serially
        cmd = 'python example_network.py {}'.format(md5)
        print(cmd)
        sp.run(cmd.split(' '))


#########
# PS1
#########

job = """#!/bin/bash
##################################################################
#SBATCH --account {}
#SBATCH --partition {}
#SBATCH --job-name {}
#SBATCH --time {}
#SBATCH -o logs/{}_stdout.txt
#SBATCH -e logs/{}_error.txt
#SBATCH -N {}
#SBATCH --ntasks {}
##################################################################
unset DISPLAY # DISPLAY somehow problematic with Slurm
srun --mpi=pmi2 python -u example_network_reconstruction.py {}
"""

# estimate run times
# (s) buffer + setup + create + connect + simulation + save
wt = 120 + 2 + 10 + 30 + 240 / 6000 * params.networkParameters['tstop'] + 30
TIME = '%i:%.2i:%.2i' % (wt // 3600,
                         (wt - wt // 3600 * 3600) // 60,
                         (wt - wt // 60 * 60))


if 'HOSTNAME' in os.environ.keys():
    if os.environ['HOSTNAME'].rfind('jusuf') >= 0 or \
            os.environ['HOSTNAME'].rfind('jr') >= 0:
        for pset in PS1.iter_inner():
            # sorted json dictionary
            js = json.dumps(pset, sort_keys=True).encode()
            md5 = hashlib.md5(js).hexdigest()

            # save parameter file
            pset.save(url=os.path.join('parameters', '{}.txt'.format(md5)))

            # create job script
            with open(os.path.join('jobs', '{}.job'.format(md5)), 'w') as f:
                f.writelines(job.format(
                    ACCOUNT,
                    PARTITION,
                    md5,
                    TIME,
                    md5,
                    md5,
                    LNODES,
                    NTASKS,
                    md5
                ))

            # figure out job dependency:
            pset_0 = pset.copy()
            for key in ['biophys', 'i_syn', 'g_eff']:
                del pset_0[key]
            md5_0 = hashlib.md5(
                json.dumps(
                    pset_0,
                    sort_keys=True).encode()).hexdigest()
            for md5_submitted, jobid in jobIDs:
                if md5_submitted == md5_0:
                    break

            # submit job and get job ID:
            cmd = ' '.join(['sbatch',
                            '--dependency=afterok:{}'.format(jobid),
                            os.path.join('jobs', '{}.job'.format(md5))
                            ]
                           )
            print(cmd, f'; depends on {md5_submitted}')
            output = sp.getoutput(cmd)
    else:
        raise NotImplementedError(
            f"do not recognize HOST {os.environ['HOSTNAME']}")
else:
    for pset in PS1.iter_inner():
        # sorted json dictionary
        js = json.dumps(pset, sort_keys=True).encode()
        md5 = hashlib.md5(js).hexdigest()

        # save parameter file
        pset.save(url=os.path.join('parameters', '{}.txt'.format(md5)))

        cmd = 'python -u example_network_reconstruction.py {}'.format(md5)
        print(cmd)
        sp.run(cmd.split(' '))


#########
# PS2
#########

job = """#!/bin/bash
##################################################################
#SBATCH --account {}
#SBATCH --partition {}
#SBATCH --job-name {}
#SBATCH --time {}
#SBATCH -o logs/{}_stdout.txt
#SBATCH -e logs/{}_error.txt
#SBATCH -N {}
#SBATCH --ntasks {}
##################################################################
unset DISPLAY # DISPLAY somehow problematic with Slurm
srun --mpi=pmi2 python -u example_network_kernel.py {}
"""

# estimate run times
# (s) buffer + setup + create + connect + simulation + save
wt = 120 + 2 + 10 + 30 + 60 + 60
TIME = '%i:%.2i:%.2i' % (wt // 3600,
                         (wt - wt // 3600 * 3600) // 60,
                         (wt - wt // 60 * 60))

if 'HOSTNAME' in os.environ.keys():
    if os.environ['HOSTNAME'].rfind('jusuf') >= 0 or \
            os.environ['HOSTNAME'].rfind('jr') >= 0:
        for pset in PS2.iter_inner():
            # sorted json dictionary
            js = json.dumps(pset, sort_keys=True).encode()
            md5 = hashlib.md5(js).hexdigest()

            # save parameter file
            pset.save(url=os.path.join('parameters', '{}.txt'.format(md5)))

            # create job script
            with open(os.path.join('jobs', '{}.job'.format(md5)), 'w') as f:
                f.writelines(job.format(
                    ACCOUNT,
                    PARTITION,
                    md5,
                    TIME,
                    md5,
                    md5,
                    LNODES,
                    NTASKS,
                    md5
                ))

            # figure out job dependency:
            pset_0 = pset.copy()
            for key in ['biophys', 't_E', 't_I', 'g_eff']:
                del pset_0[key]
            md5_0 = hashlib.md5(
                json.dumps(
                    pset_0,
                    sort_keys=True).encode()).hexdigest()
            for md5_submitted, jobid in jobIDs:
                if md5_submitted == md5_0:
                    break

            # submit job and get job ID:
            cmd = ' '.join(['sbatch',
                            '--dependency=afterok:{}'.format(jobid),
                            os.path.join('jobs', '{}.job'.format(md5))])
            print(cmd, f'; depends on {md5_submitted}')
            output = sp.getoutput(cmd)
    else:
        raise NotImplementedError(
            f"do not recognize HOST {os.environ['HOSTNAME']}")
else:
    for pset in PS2.iter_inner():
        # sorted json dictionary
        js = json.dumps(pset, sort_keys=True).encode()
        md5 = hashlib.md5(js).hexdigest()

        # save parameter file
        pset.save(url=os.path.join('parameters', '{}.txt'.format(md5)))

        cmd = 'python -u example_network_kernel.py {}'.format(md5)
        print(cmd)
        sp.run(cmd.split(' '))