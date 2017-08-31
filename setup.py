import os
import glob
from distutils.core import setup

scripts=[
    'balrog',
]
scripts=[os.path.join('bin',s) for s in scripts]

conf_files=glob.glob('config/*.yaml')

data_files=[]
for f in conf_files:
    data_files.append( ('share/balrog_config',[f]) )


setup(name="balrog", 
      version="0.1.0",
      description="The balrog code",
      license = "GPL",
      author="many",
      author_email="erin.sheldon@gmail.com",
      scripts=scripts,
      data_files=data_files,
      packages=['balrog'])
