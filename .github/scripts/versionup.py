#!/usr/bin/env python3

import argparse
import subprocess

# Runtime arguments
parser = argparse.ArgumentParser()
parser.add_argument("--major", action="store_true", help="Major version up", dest="major")
parser.add_argument("--minor", action="store_true", help="Minor version up", dest="minor")
parser.add_argument("--patch", action="store_true", help="Patch version up", dest="patch")
args = parser.parse_args()

# Get current version
process = subprocess.run("hatch version", shell=True, capture_output=True)  # (e.g.) '0.0.1'
cur_ver_str = process.stdout.decode("utf-8").strip()
ver_list = [int(i) for i in cur_ver_str.split('.')]  # (e.g.) [0, 0, 1]

# Create new version
if args.major:
    # Major version up
    ver_list[0] += 1
    # Minor version reset
    for i in range(1, len(ver_list)):
        ver_list[i] = 0
elif args.minor:
    # Minor version up
    ver_list[1] += 1
    # Patch version reset
    for i in range(2, len(ver_list)):
        ver_list[i] = 0
elif args.patch:
    # patch version up
    ver_list[2] += 1
else:  # Default: patch version up
    ver_list[2] += 1

new_ver_str = '.'.join([str(i) for i in ver_list])
# Update version with hatch command
subprocess.run(f"hatch version {new_ver_str}", shell=True)
