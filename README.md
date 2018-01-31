# Balrog-GalSim

Balrog is an image simulation tool for making randoms and characterizing sytematics by embedding fake objects into real imaging to accurately characterize measurement biases. This is an extension of the work done by [Suchyta et al.](https://arxiv.org/abs/1507.08336) for science verification data of the Dark Energy Survey (DES). Check out the original repo [here](https://github.com/emhuff/Balrog).

While Balrog is not specific to DES, the current version of this code only *supports* DES chip images and [ngmix](https://github.com/esheldon/ngmix) galaxy catalogs for injection inputs. Hopefully this generalization will be added soon; in the meantime, you can start a branch and one of the developers will help guide you on what small changes need to be made. Most changes can be made in the input GalSim config file.

## Features

Coming soon!

## Installation

The only non-standard Balrog dependancy is [GalSim](https://github.com/GalSim-developers/GalSim). While it's a non-trivial install, the [installation wiki found here](https://github.com/GalSim-developers/GalSim/blob/master/INSTALL.md) outlines the steps in great detail.

Once GalSim is installed, simply clone this repo:

```
git clone git@github.com:sweverett/Balrog-GalSim.git
```

## Running Balrog

Balrog is run with a python call to the script `balrog_injection.py`, with a few extra (required and optional) input arguments:

```
python balrog_injection.py [config_file] [input_catalog] [tile_list] [geom_file] [-t tile_dir] / 
                           [-c config_dir] [-p psf_dir] [-o output_dir] [-v --verbose]

```
(Still in progress)

* `config_file`: A GalSim config file (for now .yaml, but more types in future) that defines universal parameters for all injections. Chip-specific config parameters are appended to this file during processing. Each DES tile will have  ([See GalSim's Demo 6](https://github.com/GalSim-developers/GalSim/blob/master/examples/demo6.yaml) to see how appended configs work / look)
* `input_catalog`: Fits file containing input objects to be injected into chip images. For now only ngmix catalogs are supported (gauss, cm, or mof photometry), but the code is designed to allow other input types in future including galaxy postage stamps. Some of the standard GalSim inputs may also work, but arent' currently supported.
* `tile_list`: A txt or csv file that contains a list of all tiles that are to be processed. Delimiter can be commas or newlines.
* `geom_file`: The fits file containing tile geometry (e.g. `Y3A2_COADDTILE_GEOM.fits`).
* `tile_dir`: Directory location that contains all desired DES tile files/folders (i.e. for N tiles, the location would contain at least the N directories `/DES0001-0001`, `/DES0001-0002`, etc.). Set to `.` by default.
* `config_dir`: Directory location of the global GalSim config file. tile list file, and geometry file if not given in inputted filenames (the files only must be in the same directory if this is passed). Set to `.` by default.
* `psf_dir`: Relative directory path of a tile's psf *from* a given tile location (e.g. `{tile_dir}/{TILE}/{psf_dir}`). Set to `psfs` by default.
* `output_dir`: Location of parent output directory for Balrog images and config files. Images are saved as `{output_dir}/balrog_images/{realization}/{tilename}/{band}/{filenameTBD.fits}`, and tile configs are saved to `{output_dir}/configs/bal_config_{realization}_{tilename}.yaml`.
* `verbose`: Use -v for more verbose messages.

(more incoming...)

## Example Usage

## Contributors

* Spencer Everett
* Yuanyuan Zang
* Brian Yanny
* Nikolay Kuropatkin
* Erin Sheldon
* Eric Huff
* Ami S
* Vinicious B
* Eli Rykoff
* Megan Splettstoesser
* More to add!
