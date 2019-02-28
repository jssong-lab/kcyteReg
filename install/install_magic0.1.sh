#!/bin/bash

cwd=$(pwd)
installDir="./kcyteReg/src/"
mkdir -p $installDir
cd $installDir
git clone https://github.com/KrishnaswamyLab/MAGIC.git
cd MAGIC
git checkout b53d723c05a92405b257c07d4bb85b5ff122df55 # tag=0.1
cat src/magic/MAGIC.py | sed 's/L_t == None/L_t is None/g' > src/magic/MAGIC.py.tmp  ## fix a bug
mv src/magic/MAGIC.py.tmp src/magic/MAGIC.py
pip install . 
cd $cwd
