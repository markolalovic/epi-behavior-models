#!/bin/bash
cd src
python3 run_example.py --model baseline
python3 run_example.py --model behavior
cd ..
echo "Done"