#! /bin/bash
source /cvmfs/des.opensciencegrid.org/eeups/startupcachejob21i.sh
setup python 2.7.9+1
setup easyaccess 1.2.0+2
setup fitsio 0.9.8rc1+3
setup swarp 2.36.2+3
setup mofpsfex 0.3.2
setup psycopg2 2.4.6+7
setup despyastro 0.3.9+2
setup pil 1.1.7+13
setup pyfits 3.3+3
setup cfitsio 3.370+0
setup fitsverify
setup numpy 1.9.1+8
setup pathos
setup healpy
setup esutil 0.6.2rc1+1
setup meds 0.9.3rc2
setup pyyaml 3.11+2
setup astropy 1.1.2+4
setup pixcorrect 0.5.3+12
setup sextractor 2.23.2+4
setup despydb
setup IntegrationUtils 2.0.9+1
setup ngmix
setup covmatrix 0.9.0+1 
setup galsim
#
export BALROG_BASE=/data/des71.a/data/kuropat/balrog-base
#
export PYTHONPATH=$PYTHONPATH:/data/des71.a/data/kuropat/balrog-base/lib/python2.7/site-packages/
export DESMEDS_CONFIG_DIR=${BALROG_BASE}/desmeds-config/
#
export MEDS_DATA=/data/des71.a/data/kuropat/meds_test
export DESDATA=${BALROG_BASE}/DESDATA
export PYTHONPATH=${BALROG_BASE}/mof/ngmixer/y3v0.9.4a+1/python:$PYTHONPATH
export PATH=${BALROG_BASE}/mof/ngmixer/y3v0.9.4a+1/bin:$PATH
export PATH=${BALROG_BASE}/bin:$PATH
export medsconf="y3v02"
