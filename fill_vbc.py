import numpy as np
from drft import load_snapshot

levelmin = 8
levelmax = 16
path = './'
base = load_snapshot(path, levelmin, 'vbc')
vbc = base.load_box()

for level in range(levelmin+1, levelmax+1):
    cur = load_snapshot(path, level, 'deltab')
    off = (np.array(cur.xoff) / base.dx).astype(int)
    ex = (np.array(cur.n) / 2 ** (level - levelmin)).astype(int)

    av_vbc = np.mean(vbc[off[0]:off[0]+ex[0],
                         off[1]:off[1]+ex[1],
                         off[2]:off[2]+ex[2]])

    print(cur.xoff, cur.n, base.dx, off, ex)
    cur_vbc = np.full(cur.n, av_vbc)
    cur.write_field(cur_vbc, 'vbc')

for level in range(levelmin+1, levelmax+1):
    cur = load_snapshot(path, level, 'vbc')
    cur.plot_slice(0)
