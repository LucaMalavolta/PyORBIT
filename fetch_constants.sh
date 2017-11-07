#!/bin/sh
wget https://raw.githubusercontent.com/lucaborsato/trades/master/pytrades/constants.py
mv constants.py pyorbit/classes/constants.pyx
