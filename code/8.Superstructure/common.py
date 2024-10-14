import sys

seta = set(open(sys.argv[1], 'r').read().splitlines())
setb = set(open(sys.argv[2], 'r').read().splitlines())

if len(seta) > len(setb):
    big = seta
    small = setb

else:
    big = setb
    small = seta

for s in small:
    if s in big:
        big.remove(s)

for b in big:
    print(b)
