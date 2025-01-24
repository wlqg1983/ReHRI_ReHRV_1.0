#!/bin/bash

cp -r bin/plasmidrender $(conda info --base)/envs/MiRI_MiRIV_1.0/lib/python3.12/site-packages/
rm -r bin/plasmidrender
chmod +x bin/*
cp bin/* $(conda info --base)/envs/MiRI_MiRIV_1.0/bin
rm MiRI_MiRIV_1.0.yml
