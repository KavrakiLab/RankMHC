# RankMHC
A class-I pMHC binding mode identification tool

## Prelims

### Cloning

To set up RankMHC locally, download the repo using `git clone`.

### Frequency matrices download

Before you move on to anything else, am important step must be performed: Download the frequency matrices generated using MHCFlurry2.0 through the link below:

https://rice.box.com/s/duxshqxtkykg7u9y3b0tl79j8byc2dp2

These will be crucial for the software to run. After these are downloaded, go ahead and put those on the `helper_files` folder.

### NACCESS installation (optional but IMPORTANT!)

Part of our featurization module involves SASA/RSA values. These are extracted through `NACCESS`. Download NACCESS from here:

http://www.bioinf.manchester.ac.uk/naccess/

Follow the instructions from the NACCESS README and contact the authors of NACCESS in order to access a working key for the tool. Extract NACCESS inside the `feature_extraction` folder. As a result, you will have an `naccess2.1.1` folder inside the `feature_extraction` folder. Make sure to follow the instructions so that NACCESS is set up correctly. 

#### FreeSASA alternative:

Since NACCESS is locked behind key access, to facilitate use of RankMHC, we also provide a FreeSASA alternative. However, some RSA values will differ between NACCESS and FreeSASA, which might slightly alter results. See the following link for more information:

https://freesasa.github.io/doxygen/CLI.html

We recommend using NACCESS if possible, as RankMHC was trained on NACCESS values specifically. If not possible however, FreeSASA should provide an adequate alternative for SASA/RSA calculation. You can alternate usage of NACCESS or FreeSASA using the `--sasa_method` flag (it is set to `naccess` by default). 

## Installation

### Installing the RankMHC environment

#### Installing mamba

As a first step, you will need `mamba`, a conda-based package manager. The reason we propose mamba instead of conda is because mamba is faster and resolves package dependencies much quicker and efficient than baseline conda (conda's solvers are very slow). To install mamba in you system, follow the instructions found here:

https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html

#### Build the environment

After mamba is set, go to the RankMHC directory and type:

```
mamba env create -f installation/environment.yml
```

This will create an environment named `RankMHC` that contains most of the packages needed for running RankMHC. To activate the environment, simply execute:

```
mamba activate RankMHC
```

#### Installing PyRosetta

A final step before running RankMHC would be to install `PyRosetta` (provided that the above steps + the `NACCESS` installation (optional) were successfully completed). Obtain the free academic license for PyRosetta by applying over this link: https://els2.comotion.uw.edu/product/pyrosetta (it's free!). You will receive a username and password for downloading PyRosetta. Now you can install PyRosetta using mamba (replace the USERNAME and PASSWORD with the ones you received):

```
mamba install pyrosetta=2023.38 --channel https://USERNAME:PASSWORD@conda.rosettacommons.org
```

After this, you should be all set to use RankMHC!

## Run example 

In the `example` folder, there is an ensemble of peptide conformations modeled by `APE-Gen2.0`, our offering on class-I peptide-MHC prediction. The peptide-MHC complex in question is `8TBW` as found in the PDB database. 

To run RankMHC on this dataset, just do:

```
python RankMHC.py example/
```

After RankMHC finishes running, ranking results will appear on screen, and also as a .csv file on the `--filestore` folder (called `intermediate_files` by default but location name can change through the arguments passed to RankMHC). **IMPORTANT**: You can ignore the multiple PyRosetta outputs that you see, this is just PyRosetta being initialized multiple times (still trying to figure out a way to mute PyRosetta instantiations...).

If you're confident enough on the peptide-MHC input, you can provide those too and the peptide-MHC identification step will be bypassed (however, be certain on the inputs you provide!):

```
python RankMHC.py example/ --peptide KLSHQPVLL --MHC HLA-A*01:01
```

You can also use the FreeSASA implementation instead of NACCESS as previously discussed:

```
python RankMHC.py example/ --sasa_method freesasa
```

Finally, the alternative models as discussed in the papers can be also used, there is a set of different arguments that allows access to those. For instance, if you'd prefer to use the padding model instead of the pooling model, you can use it as so:

```
python RankMHC.py example/ --feat_transform padding
```

Please check the `helper_scripts/argparser.py` file for the available arguments and the main script for the available options. All the available models can be found in the `models` folder. 