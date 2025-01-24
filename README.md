**1. Install the software**

Download the software from https://github.com/wlqg1983/MiRI_MiRIV_1.0/archive/refs/heads/main.zip

unzip MiRI_MiRIV_1.0-main.zip

cd  MiRI_MiRIV_1.0-main

conda env create -f  MiRI_MiRIV_1.0.yml

conda activate MiRI_MiRIV_1.0

sh Install.sh

rm Install.sh


**2. Script of searching subconfigurations**

python bin/MiRI.py -c MiRI.config.ini


**3. Draw the map of recombination organelle genome**

python bin/MiRIV.py -c MiRIV.config.ini




