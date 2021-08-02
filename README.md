# Balrog-GalSim

Balrog is an image simulation tool for making randoms and characterizing systematics by embedding fake objects into real imaging to accurately characterize measurement biases. This is an extension of the work done by [Suchyta et al.](https://arxiv.org/abs/1507.08336) for science verification data of the Dark Energy Survey (DES). Check out the original repo [here](https://github.com/emhuff/Balrog).

While Balrog is not specific to DES, the current version of this code only *supports* DES chip images and [ngmix](https://github.com/esheldon/ngmix) galaxy catalogs for injection inputs. Hopefully this generalization will be added soon; in the meantime, you can start a branch and one of the developers will help guide you on what small changes need to be made. Most changes can be made in the input GalSim config file.

The archived repository for the DES Year 3 (Y3) cosmology analysis calibration can be found below:
https://zenodo.org/badge/DOI/10.5281/zenodo.5154103.svg

## Installation

The only non-standard Balrog dependency is [GalSim](https://github.com/GalSim-developers/GalSim). While it's a non-trivial install, the [installation wiki found here](https://github.com/GalSim-developers/GalSim/blob/master/INSTALL.md) outlines the steps in great detail. (edit - Galsim 2.0+ can now be installed with pip/conda!)

Once GalSim is installed, simply clone this repo:

```
git clone git@github.com:sweverett/Balrog-GalSim.git
```

## Running Balrog

(Updated as of 4/5/18) Balrog is run with a python call to the script `balrog_injection.py`, with a few extra (required and optional) input arguments that can either be passed by command line:

```
python balrog_injection.py [config_file] [-l tile_list] [-g geom_file] [-t tile_dir] [-c config_dir] \
                           [-p psf_dir] [-o output_dir] [-v --verbose]

```
or explicitly in the `config_file` described below:

* `config_file`: A GalSim config file (for now .yaml, but more types in future) that defines universal parameters for all injections. Chip-specific config parameters are appended to this file during processing. Each DES tile will produce a separate yaml config that handles all injections for all chips in all bands that overlap with the tile ([see GalSim's Demo 6](https://github.com/GalSim-developers/GalSim/blob/master/examples/demo6.yaml) to see how appended configs work / look). An example balrog config file is given in `configs/bal_config.yaml`. It also includes a description of how to pass these command-line input arguments in the config itself.
* `tile_list`: A txt or csv file that contains a list of all tiles that are to be processed. Delimiter can be commas or newlines.
* `geom_file`: The fits file containing tile geometry (e.g. `Y3A2_COADDTILE_GEOM.fits`).
* `tile_dir`: Directory location that contains all desired DES tile files/folders (i.e. for N tiles, the location would contain at least the N directories `/DES0001-0001`, `/DES0001-0002`, etc.). Set to `.` by default.
* `config_dir`: Directory location of the global GalSim config file. tile list file, and geometry file if not given in inputted filenames (the files only must be in the same directory if this is passed). Set to `.` by default. Cannot be set in `config_file` (for obvious reasons).
* `psf_dir`: Relative directory path of a tile's psf *from* a given tile location (e.g. `{tile_dir}/{TILE}/{psf_dir}`). Set to `psfs` by default.
* `output_dir`: Location of parent output directory for Balrog images and config files. Images are saved as `{output_dir}/balrog_images/{realization}/{tilename}/{band}/{chipname}_balrog_inj.fits`, and tile configs are saved to `{output_dir}/configs/bal_config_{realization}_{tilename}.yaml`.
* `verbose`: Use -v for more verbose messages.

### Required Directory Structures

Note that the current version of Balrog only supports DES injections. Part of this requirement is due to an assumption of how tiles and single exposure chips are are named and structured with respect to one another. Inside the `tile_dir` there should be a collection of directories with DES tilenames with each housing the nullweight chip images in their respective `{tilename}/{nullwt-{band}}` directories. Visually:

```
Balrog-Galsim
│   README.md    
│
└───balrog
|      ...
|
└───config
|      ...
└───inputs
│   │
│   └───tiles
|       |   
│       └───DES2329-5622
|       |    |
|       |    └───nullwt-g
|       |    |      gchip1.fits
|       |    |      gchip2.fits
|       |    |      ...
|       |    └───nullwt-r
|       |    |      rchip1.fits
|       |    |      rchip2.fits
|       |    |      ...
|       |    └───nullwt-i
|       |    |      ...
|       |    └───nullwt-z
|       |    |      ... 
|       |    └───psfs
|       |           psf1.fits 
|       |           psf2.fits 
|       |           ...
|       └───DES2349+1334
|       |       ...
|       └───DES0744+1126
|       |       ...
|       └───...
|   
└───output_dir
    |
    └───conifgs
    |       ...
    └───balrog_images
        |
        └───DES2329-5622
        └───DES2349+1334
        ...
```

### Example Usage

Let's say you are running the standard DES Y3 setup for Balrog. Then running from the repo home you would use the following values for the above inputs:
* `config_file = configs/bal_config.yaml`
* `tile_list = inputs/tilelist.csv`
* `geom_file = inputs/Y3A2_COADDTILE_GEOM.fits`
* `tile_dir = inputs/tiles`
* `psf_dir = psfs` (the default option, so not needed)
* `output_dir = outputs/` (or whatever you would like!)

and so the terminal command would be

```
python balrog/balrog_injection.py config/bal_config.yaml inputs/tilelist.csv inputs/Y3A2_COADDTILE_GEOM.fits \
-t inputs/tiles -o outputs/
```

Alternatively, these values could have been set in the `image` field of the `config_file` (except for `output_dir`; see example file) in which case you would simply type

```
python balrog/balrog_injection.py config/bal_config.yaml
```

## Balrog Config

The global balrog config file is very similar to yaml config files used for GalSim. Check out the [GalSim demo page](https://github.com/GalSim-developers/GalSim/wiki/Tutorials) for lots of examples of how to use the yaml files and how they are translated to python scripts. The important parts are summarized in the comments of the example config [bal_config.yaml](https://github.com/sweverett/Balrog-GalSim/blob/master/config/bal_config.yaml). (**Note**: You will need to change a few path names in the config file to work for your local machine!) 

However, note that a Balrog config will **not** run successfully if called by the `galsim` executable; it only houses the global simulation variables while the individual chip injections parameters are set during the `balrog_injection.py` script. Each tile produces its own complete multi-output yaml config file that is sent to `galsim` at the end of processing.

There are a few config inputs specific to the `image` field of Balrog configs that are worth highlighting here:
* `n_objects`: The total number of objects to be injected per DES tile **per realization**.
* `object_density`: The injected object density in the tile field **per realization**.
* `n_realizations`: The number of injection realizations used to reach the desired galaxy count or density.

Two things to note: (1) **Only one** of `n_objects` or `object_density` is allowed as an input; not both! (2) Either input should give the desired count or density **per realization**!. An older version of the code had this set to the desired final count/desnity, but this was counter-intuitive for users.

## Input Catalogs

(more later - for now, `ngmix_catalog`, `meds_catalog`, `des_star_catalog`. See `balinput.py` and `balobject.py`)

(Fits file containing input objects to be injected into chip images. For now only ngmix catalogs are supported (gauss, cm, or mof photometry), but the code is designed to allow other input types in future including galaxy postage stamps. Some of the standard GalSim inputs may also work, but arent' currently supported.)

## More to come...

## Contributors

* Spencer Everett (contact at sweveret@ucsc.edu)
* Yuanyuan Zhang
* Brian Yanny
* Nikolay Kuropatkin
* Erin Sheldon
* Eric Huff
* Ami Choi
* Vinicious Busti
* Eli Rykoff
* Megan Splettstoesser
* More to add!
