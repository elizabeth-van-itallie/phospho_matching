#!/usr/bin/env python3

from pathlib import Path

def set_up_folders(folders):

    for f in folders:
        Path(f).mkdir(exist_ok = True)

if __name__ == '__main__':

    folders = snakemake.params.folders
    set_up_folders(folders)
