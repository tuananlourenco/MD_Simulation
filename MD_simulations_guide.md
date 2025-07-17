# How to Perform MD Simulations - update July 2025

In this file we show the procedures and organization that we use to perform our MD simulations. Keep in mind that this file was created at June, 2025. Then, it can be **outdated**.
First, we recommend the reading of the follow articles, which summarizes best practices and a short introduction about the methodology:

- (i) [Foundations of MD simulations](https://livecomsjournal.org/index.php/livecoms/article/view/v1i1e5957)
- (ii) [Errors and limitations]()
- (iii) [Best Practices in Transport Properties](https://livecomsjournal.org/index.php/livecoms/article/view/v1i1e6324)
- (iv) [Best Practices in Error Quantification](https://livecomsjournal.org/index.php/livecoms/article/view/v1i1e5067)
- (iv) [Size effects on transport properties](https://pubs.aip.org/aip/jcp/article/156/13/134705/2841128/Investigating-finite-size-effects-in-molecular)
- (iv) [Force Field development (The case of CLP&P)](https://link.springer.com/article/10.1007/s00214-012-1129-7)


## Directory Organization

The idea of this file is to provide overview, recommendation and scripts to create the main directories and files.
First, you need to understand the directory organization, which is organized as follows:

- Main directory: **Project_{overleaf_name}_{start_year}_MD_simulations**
  - Sub-directories-first-class: 
    - **Protocol_input_file** (which have the *.mdp files)
    - **Force_Field_input_files** (all *.itp files)
    - **User_Notes_References** (you can use this repository to organize your schedule, add references and project notes of the meetings)
    - **BOX** (here you have the packmol.inp based on your choices in the user_variables.txt input. But keep in mind that you need to add the molecules.pdb files by yourself)
    - **RUN_{system}_{salt_concentration}_{system_size}_{T_production}** (directory for simulations)
      - Sub-directories-second-class:
        - **1-min_integrator-{integrator}_emtol-{emtol}_nsteps-{nsteps}_emstep-{emstep}**
        - **2-NVT_thermalization_integrator-{integrator}_dt_{dt}_time_{time_ps}ps_T-{ref-t}_{tcoupl}**
        - **3-NPT_ergodicity_integrator-{integrator}_dt_{dt}_time_{time_ps}ps_T-{t-initial}-{t-final}-K_P-{ref-p}-bar_{tcoupl}_{pcoupl}**
        - **4-NPT_equilibration_integrator-{integrator}_dt_{dt}_time_{time_ps}ps_T-{ref-t}-K_P-{ref-p}-bar_{tcoupl}_{pcoupl}**
        - **6-NVT_re-equilibrium_integrator-{integrator}_dt_{dt}_time_{time_ps}ps_T-{ref-t}-K_{tcoupl}**
        - **7-NVT_production_integrator-{integrator}_dt_{dt}_time_{time_ps}ps_T-{ref-t}-K_{tcoupl}_dump-{nstxout-compressed}**
            - Each of these directories has a resuls and repository directory to organize the outputs.

Note that the name of the directories are based on the options and methodes used in the MD simulations, this is done to ensure the understanding of the protocol used in each MD simulation Step.

## Description of MD Simulation Workflow

Now let's talk about how we should proceed to prepare and run the simulations. We are not talking about the scripts here, only discussing the procedures.

### First step, force field and molecules/ions chemical structures:

In general we use OPLS-based force fields, such as the CL&P and the OPLS-2009IL. To create the force files (*itp) we recommend the use of the following tools and databases. Try to use this order in your change, for example, if you did not find your molecule/ion in the first tool go to the second and so on.

- [DLPGEN/OPLS/CL&P](https://webpages.ciencias.ulisboa.pt/~cebernardes/dlpgen_prog/Software_dlpgen.html) (Carlos Bernardes Group) - Creates *itp and *pdb inputs based on OPLS and CL&P force field, the software is available in its website, has several inputs for usage and there are some youtube videos from Dr. Carlos Bernardes explain how to use the code.

- [fftool/clandp (Padua Group)](https://github.com/paduagroup) - Creates *itp and *pdb inputs based on OPLS and CL&P force field, the github has a good tutorial and it is well organized.
    
- [OPLS-2009IL](https://github.com/orlandoacevedo?tab=repositories) - Provides several *itp and *gro/pdb files for ionic liquids systems. The files were created by Orlando Acevedo group using the default OPLS parameters.

- [LIGPARGEN](https://zarbi.chem.yale.edu/ligpargen/) - Online website that uses your *pdb file to create gromacs *itp files based on default OPLS force field. Here you need to ensure that your *pdb file is relaxed using ab initio methods and also that the atom numbers in the output and input are the same. We recommed to use *gro or *pdb outputs from the LIGPARGEN when running your simulations. This choice should be used only as final option because the *itp create are more "complex" to understand.

If you need to create *pdb files from scratch you need to remember to relax the structure using ab initio methods as recommeded in the force field papers. Also in almost all tools above you need to convert *xyz files to *pdb files, which in general we do using openbabel and create some changes in the atom names (openbabal only consider periodic table). Then, it helps to prevent problems if you ensure manually that both *itp/*top files and your final *pdb have the same atom-name or atomtype, this prevent the tons of warnings in the beggining of the simulation. Also, always draw the connections by yourself to check if everything is OK.

### Second step, initial MD boxes structures:

In this step we use [PACKMOL from Prof. Dr. Leandro Martinez](https://m3g.github.io/packmol/) to create the cubic boxes using the desired system composition. PACKMOL has a simple usage and ensures a random starting point. To run PACKMOL we need the *pdb files for each molecule/ion and the packmol.inp file, which has the following structure:

    tolerance 2.5 -> minimal distance between two non-bonded atoms
    filetype pdb -> output and inputs type, which must be the same from the inputs
    randominitialpoint -> Ensures a random configuration
    add_box_sides 0.2 -> Extra distances added to the box sides. This helps to prevent atom overposition
    output box.pdb -> output type, which must be the same from the inputs

    structure molecule.pdb # PDB File name for you molecule/ion
        number XXX # Number of the above species
        inside box x_initial y_initial z_initial x_final y_final z_final # Box size in angstrom
    end structure -> Just finishes the loop

    seed -1 -> PACKMOL seed (-1 ensures random numbers every time)

If you have more than one molecule/ion in your simulation you need to repeat the structure loop for every molecule. For example, if you have 2 molecules, you need to have two structure loops, describing the file, number and the size of the box. If you are working with bulk systems the "inside box" option must be the same for every molecule to ensure that you have a real random mixture.

Keep in mind that general PACKMOL creates boxes that are not perfectly cubic, then we recommend to check all the final box.pdb outputs to ensure that all the box sides have the same lenght. Also, it is good to add an extra 0.5 angstrom in the final sides to ensure that atom overlaps are not going to happen.

Second, if PACKMOL find troubles to create your box, maybe you should use larger box sides. A good recommendation is to use the [**PACKMOL VOLUME GUESSER**](https://m3g.github.io/packmol/nmols.shtml) from PACKMOL site to find a good initial box side. Also, it is always recommended to start the MD simulation with boxes that are larger than the expected equilibrated system. Then, I recommend to add 5 angstrom to the values suggested from **PACKMOL VOLUME GUESSER** tool. If you do not have any information about your equilibrium density or experimental, try to use values from similar systems.

Third, remember that in MD we relly on statystics, then, very small systems will have problems or not good results.

### Third step, MD simulation protocol:

Our MD protocol is based on a seven steps approach:

**1-min** - Initial hotspot eliminiation using Steep Descent algorithm. - Here we use an energy minimization algorithm to eliminate any hotspot in the box, i.e., atoms with high residual forces.
**2-NVT_thermanlization** - Initial NVT Thermalization. - In this step we aim to thermalize the system at a higher temperature to remove any remanied hotspot and ensures that we are not using the structure from PACKMOL anymore.
**3-NPT_annealing.** - NPT annealing equilibration using a stronger barostat and thermostat. Here we are going to increase and decrase the temperature around our target value. This help the system to traverse the potential energy surface. Also, some systems like ionic liquids can present metastable phases depending of the temperature rate equilibration.
**4-NPT_Equilibration** - NPT equilibration at fixed temperature using a soft barostat/thermostat, such as Parrinello-Rahman and Nose-Hoover. Here we are going to ensure the equilibration process and also the real ensemble (read the articles recommended in the beggining of this file). This step provide us the equilibrated volume for production.
**5-Box_equilibrium_resize** - Before proceed to the next step you need to check the equilibration of density, temperature, potential energy and pressure. For that you can use gmx energy, the equilibration_detection.py, xmgrace and look to the final box structure. Remeber that you must look for half of the data produced in the previous step. If the system is not equilibrated, you need to extend the previous step for longer times using the final configuration. Once the system is equilibrated, you need to obtain the equilibrated box density and use the gmx editconf tool resize the final NPT Equilibrated box to the average density.
**6-NVT_requilibration** - Here we perform a NVT requilibration and thermalization to the resized box obtained from gmx editconf. The idea here is to prepare the box to the production stage.
**7-NVT_production** - Now you start your NVT production simulation in which you are going to calculate the properties and use it to write your article. Keep in mind that the simulation lenght must be enough to converge all the desired properties. 

### How to use the codes and softwares in the MD procedure

#### PACKMOL

To run packmol you only need the *pdb files and the *packmol.inp* input. To execute the software you use the command:

    packmol < packmol.inp > packmol.out

Once the loop is over and you see the follwoing message, everything is ok:

################################################################################

                                 Success! 
              Final objective function value: .15800E-01
              Maximum violation of target distance:   0.008085
              Maximum violation of the constraints: .77707E-02

--------------------------------------------------------------------------------

              Please cite this work if Packmol was useful: 

           L. Martinez, R. Andrade, E. G. Birgin, J. M. Martinez, 
         PACKMOL: A package for building initial configurations for
                   molecular dynamics simulations. 
        Journal of Computational Chemistry, 30(13) pp. 2157-2164, 2009.
                  https://doi.org/10.1002/jcc.21224

################################################################################


Now, lets ensure that the box is really cubic. Open the output and check CRYST1 line, for example:

    CRYST1    60.21    60.27    60.21  90.00  90.00  90.00 P 1           1

You can see that the box sides have different lenght, then you should make all of them (based on the largest one) and add a 0.5 angstrom:

    CRYST1    60.77    60.77    60.77  90.00  90.00  90.00 P 1           1

#### GROMACS

Once you have all the inputs, you need to start and run the MD. To start, we need to compile everything using the *gmx grompp* as follows:

    gmx grompp -f MDP_FILE.mdp -c BOX_CONFIG.pdb/gro -p topology.top -o OUTPUT_NAME

The above command create the *OUTPUT_NAME.tpr* that will be used to run the simulation using the *gmx mdrun*:

    gmx mdrun -deffnm OUTPUT_NAME 

The above command runs the simulation and all the outputs will have the name *OUTPUT_NAME****, for example *OUTPUT_NAME.gro*, *OUTPUT_NAME.edr*, *OUTPUT_NAME.xtc*, *OUTPUT_NAME.cpt* and etc. Also, you should add optimization and paralelization flags in this command depending of your HPC resources.

To check the energies and convergence we use the gmx energy, which has a big menu in which you can chose different energies. To run we do:

    gmx energy -f OUTPUT_NAME.edr -o ENERGY_OUTPUT -xvg none

From that you will see something as:

Let's say that you are going to check the density and potential energy from your NPT, then, you should type the number of the options and once the calculation finishes (in general some seconds) you will have a *ENERGY_OUTPUT.xvg* file that you can open in xmgrace. 

To check the energy convergence you can use the **equilibrium_detection.py** script based on [*detect_equilibration.py* script](https://github.com/choderalab/automatic-equilibration-detection/tree/master), which uses the [methodology developed by John D. Chodera](https://pubs.acs.org/doi/10.1021/acs.jctc.5b00784). This script is provided in my github and to run you just need to type:

    python **equilibrium_detection.py** file.xvg or use --help for more information

Once you know the average equilibrated density for your system, then, you need to rescale the final *4-NPT_equilibrium.gro* configuration. To do that we use the gmx editconf as follows:

    gmx editconf -f *4-NPT_equilibrium.gro* -density XXXX -o 5-Box_rescaled-average-density.gro

I recommend to check the *gmx editconf -h* to understand what is being done. The 5-Box_rescaled-average-density.gro file will be your starting point in the NVT production.

Ok, now that we have discussed the main points when running and preparing your MD simulations, we can proceed to understand how the MD_setup.py script works.

## MD Seting-up Script

The script **MD_setting_up.py** prepare the directories for your simulations, write the *mdp files using the procedure described above, create the PBS job.md files (which are separated by stages), creates the packmol input and a latex table containing all the variables used in each mdp file. However, keep in mind that you **MUST CHECK ALL THE FILES AND DIRECTORIES** by yourself. Also, we are not going to create the *itp and individual *pdb files here, this is your job and can be done using the tools described above.

So, first, lets understand how the script works:

- First, you define all the importante variables in your work, like system name, size, temperature, simulation lenght, pressure, and etc. That means, we need all the variables described within the  **user_variables.txt** file. Also, to chose the variables you need to understand your goals and know a prior information of your systems. For that, check the previous literature about the material or similar materials and also similar works to your topic. I mean, general literature, not only your scientific group.

- Once you ran the script, checked all the directories, *mdp files, job.md.N files and packmol.inp create in BOX folder. You need to add all your *itp and *pdb files to the ***Force_Field_input_files*** and ***BOX*** and run packmol as suggested above.
- Check all the job.md.N files regarding the queue, ncpus, -ntomp and -np 96 gmx_mpi. Remember that **-ntomp X -np MUST BE EQUAL ncpus** this is a critical part and if wrong can create problem to other jobs from other users in the cluster.
- Finally, you can run your simulations. Remember, in each step you need to check how the simulation performed, look to your box structure, to the energies. Some files are created on the fly and are storaged in the *-results folder within the MD stage directory.
