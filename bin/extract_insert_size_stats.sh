#! /usr/bin/env bash
LOG=$1

grep 'insert size' $1 | tail -n 1 | grep -oP '(?<=insert size: )\d+\.\d+|(?<=sd: )\d+\.\d+'
