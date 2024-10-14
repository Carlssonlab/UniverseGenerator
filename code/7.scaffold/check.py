import sys

for line in open(sys.argv[1], 'r').read().splitlines():

    if "CC(=O)Nc1noc2c1C([1Re])C([Re])([1Re])C([Re])([2Re])C2[2Re]" in line:

        print(line.strip())
