#!/bin/sh
#
# Bash script for running the ADIF summarizing program adifstats.
#
PYTHON="python3"
adif_filename = "/home/frejon/.local/share/WSJT-X/wsjtx_log.adi"

$PYTHON adifstats.py $adif_filename
