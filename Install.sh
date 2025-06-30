#!/bin/bash

cp -r bin/plasmidrender $(conda info --base)/envs/ReHRI_ReHRV_1.0/lib/python3.12/site-packages/
rm -r bin/plasmidrender
chmod +x bin/*
rm ReHRI_ReHRV_1.0.yml

