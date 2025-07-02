# How to Perform MD Simulations - update June 2025

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

- Main directory: **system_md_ff_software_size_special-case**
  - Sub-directories-first-class: 
    - **Protocol** (which have the *.mdp files)
    - **Force_Field** (all *.itp files)
    - **RUN_system_salt-concentration_temperature-special-case** (directory for simulations)
      - Sub-directories-second-class:
        - **min_etol_nstep**
        - **NVT_thermalization_dt_tim**
        - **NPT_annealing_dt_time_thermostat_barostat_Tinitial_Tfinal_Pressure**
        - **NPT_equilibrium_dt_time_thermostat_barostat_Tfinal_Pressure**
        - **NVT_production_dt_time_thermostat_T_dump**
            - Sub-directories-third-class
              - **Analysis**

Note that the name of the directories are based on the options and methodes used in the MD simulations, this is done to ensure the understanding of the protocol used in each MD simulation Step.

## Description of MD Simulation Workflow

Now let's talk about how we should proceed to prepare and run the simulations. We are talking about the scripts here, only discussing the procedures.

### First step, force field and molecules/ions chemical structures:

In general we use OPLS-based force fields, such as the CL&P and the OPLS-2009IL. To create the force files (*itp) we recommend the use of the following tools and databases:

- [fftool/clandp (Padua Group)](https://github.com/paduagroup) - Creates *itp and *pdb inputs based on OPLS and CL&P force field, the github has a good tutorial and it is well organized.
    
- [DLPGEN/OPLS/CL&P](https://webpages.ciencias.ulisboa.pt/~cebernardes/dlpgen_prog/Software_dlpgen.html) (Carlos Bernardes Group) - Creates *itp and *pdb inputs based on OPLS and CL&P force field, the software is available in its website, has several inputs for usage and there are some youtube videos from Dr. Carlos Bernardes explain how to use the code.

- [OPLS-2009IL](https://github.com/orlandoacevedo?tab=repositories) - Provides several *itp and *gro/pdb files for ionic liquids systems. The files were created by Orlando Acevedo group using the default OPLS parameters.

- [LIGPARGEN](https://zarbi.chem.yale.edu/ligpargen/) - Online website that uses your *pdb file to create gromacs *itp files based on default OPLS force field. Here you need to ensure that your *pdb file is relaxed using ab initio methods and also that the atom numbers in the output and input are the same. We recommed to use *gro or *pdb outputs from the LIGPARGEN when running your simulations. This choice should be used only as final option because the *itp create are more "complex" to understand.

If you need to create *pdb files from scratch you need to remember to relax the structure using ab initio methods as recommeded in the force field papers.

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

Keep in mind that general PACKMOL creates boxes that are not perfectly cubic, then we recommend to check all the final box.pdb outputs to ensure that all the box sides have the same lenght. Also, it is good to add an extra 0.2 angstrom in the final sides to ensure that atom overlaps are not going to happen.

Second, if PACKMOL find troubles to create your box, maybe you should use larger box sides. A good recommendation is to use the [volume guesser tool](https://m3g.github.io/packmol/nmols.shtml) from PACKMOL site to find a good initial box side. Also, it is always recommended to start the MD simulation with boxes that are larger than the expected equilibrated system. If you do not have any information about your equilibrium density or experimental, try to use values from similar systems.

Third, remember that in MD we relly on statystics, then, very small systems will have problems or not so good results.

### Third step, MD simulation protocol:

Our MD protocol is based on a five steps approach:

1. Potential Energy Minimization using Steep Descent algorithm. - The idea here is just avoid a bad initial starting point.
2. Initial NVT Thermalization. - In this step we aim to thermalize the system at a higher temperature to remove any hot spot and ensures that we are not using the structure from PACKMOL anymore.
3. NPT Annealing Equilibration. - NPT annealing equilibration using a stronger barostat and thermostat. Here we are going to increase and decrase the temperature around our target value. This help the system to traverse the potential energy surface. Also, some systems like ionic liquids can present metastable phases depending of the temperature rate equilibration.
4. NPT Equilibration at target Temperature. - NPT equilibration at fixed temperature using a soft barostat/thermostat, such as Parrinello-Rahman and Nose-Hoover. Here we are going to ensure the equilibration process and also the real ensemble (read the articles recommended in the beggining of this file).
        Before proceed to the next step you need to check the equilibration of density, temperature, potential energy and pressure. For that you can use gmx energy, detect_equilibration.py, xmgrace and look to the final box structure. If the system is not equilibrated, you need to extend the XXX for longer times using the final configuration. Once the system is equilibrated, you need to obtain the equilibrated box volume and use the gmx editconf tool resize the final NPT Equilibrated box to the average volume.
5. NVT Production at target Temperature. - Now you start your NVT production simulation in which you are going to calculate the properties and use it to write your article. Keep in mind that the simulation lenght must be enough to converge all the transport properties. In addition, you need to add extra 5 ns in this time, which will be used as "thermalization NVT equilibration" and discarded in the end. For example, if you need 20 ns of production, your simulation must have 25 ns and you are going to use only from 5 - 25 ns to calculate the properties.

### How to use the codes and softwares in the MD procedure

#### PACKMOL

To run packmol you only need the *pdb files and the *packmol.inp* input. To execute the software you use the command:

    packmol < packmol.inp > packmol.out

Once the loop is over and you see the follwoing message, everything is ok:

Now, lets ensure that the box is really cubic. Open the output and check CRYST1 line, for example:

    CRYST1    60.21    60.27    60.21  90.00  90.00  90.00 P 1           1

You can see that the box sides have different lenght, then you should correct and add a 0.2 angstrom:

    CRYST1    60.47    60.47    60.47  90.00  90.00  90.00 P 1           1

#### GROMACS

Once you have all the inputs, you need to start and run the MD. To start, we need to compile everything using the *gmx grompp* as follows:

    gmx grompp -f MDP_FILE.mdp -c BOX_CONFIG.pdb/gro -p topology.top -o OUTPUT_NAME

The above command create the *OUTPUT_NAME.tpr* that will be used to run the simulation using the *gmx mdrun*:

    gmx mdrun -deffnm OUTPUT_NAME 

The above command runs the simulation and all the outputs will have the name *OUTPUT_NAME****, for example *OUTPUT_NAME.gro*, *OUTPUT_NAME.edr*, *OUTPUT_NAME.xtc*, *OUTPUT_NAME.cpt* and etc. Also, you should add optimization and paralelization flags in this command depending of your HPC resources.

To check the energies and convergence we use the gmx energy, which has a big menu in which you can chose different energies. To run we do:

    gmx energy -f OUTPUT_NAME.edr -o ENERGY_OUTPUT

From that you will see something as:


Let's say that you are going to check the density and potential energy from your NPT , then, you should type the number of the options and once the calculation finishes (in general some seconds) you will have a *ENERGY_OUTPUT.xvg* file that you can open in xmgrace. If you want, you can add the flag *-xvg none* to have only the values of the properties without the xmgrace plot configuration.

To check the energy convergence you can use first the xmgrace and look to "linear" regions, that can be confirmed by a linear regression and also the [*detect_equilibration.py* script](https://github.com/choderalab/automatic-equilibration-detection/tree/master), which uses the [methodology developed by John D. Chodera](https://pubs.acs.org/doi/10.1021/acs.jctc.5b00784). This script is provided in my github and to run you need first have the *equilibration.py*, *detect_equilibration.py* in the same directory. To run you only need to call the detect_equilibration and provide your *ENERGY_OUTPUT.xvg* file containing only the property that you want to check the equilibration.

    python detect_equilibration.py energy.xvg

Once your system is equilibrated, you must go back to the gmx energy and calculate the density average value from your 4-NPT_equilibrium_fixed_temperature step within the equilibrated time range obtained from the detect_equilibration.py script. To do that you need to use the flags *-b* and *-e* in *gmx energy* to provide your time range. For example, let's say that you are going to calculate only the values between the 2 and 4 ns of the simulation. Then you do:

    gmx energy -f OUTPUT_NAME.edr -o ENERGY_OUTPUT -b 2000 -e 4000 

Ok, now you know the average equilibrated density for your system and you need to rescale the final *4-NPT_equilibrium_fixed_temperature.gro* configuration. To do that we use the gmx editconf as follows:

    gmx editconf -f ENERGY_OUTPUT.gro -density XXXX -o NPT_equilibrium_fixed_temperature-average-density.gro

I recommend to check the *gmx editconf -h* to understand what is being done. The NPT_equilibrium_fixed_temperature-average-density.gro file will be your starting point in the NVT production.


Ok, now that we have discussed the main points when running and preparing your MD simulations, we can proceed to the scripts and work.

## MD Seting-up Script

The script **MD_setting_up.py** prepare the directories for your simulations, write the *mdp files using the procedure described above and also organize some of the inputs. However, keep in mind that you **MUST CHECK ALL THE FILES AND DIRECTORIES** by yourself. Also, we are not going to create the *itp files here, this is your job and can be done using the tools described above.

So, first, lets understand how the script works:

- First, you define all the importante variables in your work, like system name, size, temperature, simulation lenght, pressure, and etc. That means, we need all the variables described within the  **MD_setting_up.py** file. Also, to chose the variables you need to understand your goals and know a prior information of your systes. For that, check the previous literature about the material or similar materials and also similar works to your topic. I mean, general literature, not only your scientific group.

- Second, you need to add all your *itp and *pdb files to the ***Original_FF_Folder*** and ***Original_PDB_Folder*** and provide to the script the location of these directories. The script will copy these files from the folders and create everything based that. Then, ensures that everything is correct and that you are using the same name for all the connected files. For example, if you water pdb file is caleed H2O.pdb, the itp must be H2O.itp and you need to have H2O in your system name.

- Finally, you can run the **MD_setting_up.py** and check all the information. First, check if the *mdp files in the Protocol folder represents what you want. Then, do the same for the Force_Field files, topology.top files and all the directories.

- In the folder HPC_Script, you can find some examples of PBS job submission files, but, you need to change it based on your pourposes.  















