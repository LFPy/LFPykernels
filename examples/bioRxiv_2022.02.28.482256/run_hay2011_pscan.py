import os
import subprocess as sp
import json
import hashlib
from parameters import ParameterSpace, ParameterSet, ParameterRange
import hay2011_network_parameters as params

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
            weight_scaling=ParameterRange([0.975, 1., 1.025, 1.050, 1.075]),
            n_ext=ParameterRange([[920, 160]]),
        )
    )
)

PS0.save('hay2011_PS0.txt')

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
            n_ext=ParameterRange([[920, 160]]),
            g_eff=ParameterRange([False, True]),
            perseg_Vrest=ParameterRange([False])
        )
    )
)

PS1.save('hay2011_PS1.txt')

# Parameterspace for kernel signals generation
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
            n_ext=ParameterRange([[920, 160]]),
            t_E=ParameterRange([200.]),
            t_I=ParameterRange([400.]),
            g_eff=ParameterRange([True]),
            perseg_Vrest=ParameterRange([False])
        )
    )
)

PS2.save('hay2011_PS2.txt')


##############
# Singularity prompt
###############
i = 0
while i < 2:
    answer = input(
        'Run jobs using Singularity container "singularity.sif" built from Dockerfile (see README)? y/N: ')
    if any(answer.lower() == f for f in ["y", 'Y', '1']):
        print("Yes")
        singularity = True
        break
    elif any(answer.lower() == f for f in ['N', 'n', '0']):
        print("No")
        singularity = False
        break
    else:
        i += 1
        if i < 2:
            print('Please type "Y" or "n"')
        else:
            print("Nothing done")
            singularity = False

if singularity:
    singularity_stuff = [
        'module --force purge\nmodule load Stages/2022 GCCcore/.11.2.0 Apptainer-Tools/2022 GCC/11.2.0 ParaStationMPI/5.5.0-1',
        'singularity exec lfpykernels.sif']
else:
    singularity_stuff = [None, None]


########
# PS0
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
wt = 180 + 2 + 10 + 240 + 720 / 2000 * params.networkParameters['tstop'] + 180
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
{}
srun --mpi=pmi2 {} python -u hay2011_network.py {}
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
                    singularity_stuff[0],
                    singularity_stuff[1],
                    md5
                ).replace('None', ''))
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
        cmd = 'python hay2011_network.py {}'.format(md5)
        print(cmd)
        sp.run(cmd.split(' '))


#######
# PS1
######

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
{}
srun --mpi=pmi2 {} python -u hay2011_network_reconstruction.py {}
"""

# estimate run times
# (s) buffer + setup + create + connect + simulation + save
wt = 180 + 2 + 10 + 60 + 3000 / 6000 * params.networkParameters['tstop'] + 180
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
                    singularity_stuff[0],
                    singularity_stuff[1],
                    md5
                ).replace('None', ''))

            # figure out job dependency:
            pset_0 = pset.copy()
            for key in ['biophys', 'i_syn', 'g_eff', 'perseg_Vrest']:
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

        cmd = 'python -u hay2011_network_reconstruction.py {}'.format(md5)
        print(cmd)
        sp.run(cmd.split(' '))


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
{}
srun --mpi=pmi2 {} python -u hay2011_network_kernel.py {}
"""

# estimate run times
# (s) buffer + setup + create + connect + simulation + save
wt = 180 + 5 + 10 + 60 + 720 / 600 * 600 + 60
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
                    singularity_stuff[0],
                    singularity_stuff[1],
                    md5
                ).replace('None', ''))

            # figure out job dependency:
            pset_0 = pset.copy()
            for key in ['biophys', 't_E', 't_I', 'g_eff', 'perseg_Vrest']:
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

        cmd = 'python -u hay2011_network_kernel.py {}'.format(md5)
        print(cmd)
        sp.run(cmd.split(' '))
