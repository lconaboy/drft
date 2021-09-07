"""
A tool to downsample an IC field and average over it.
"""

import numpy as np
import matplotlib.pyplot as plt

import grafic_tools as grafic

def gen_fake():
    """Generates a cube of fake test data, with 8 different values. To use
    as test, set level_s to 1 -- should reproduce exactly the same cube,
    just downsampled."""
    
    split = 2
    ss = n // split
    x = np.ones((ss, ss, ss))
    xtop = np.hstack((1*x, 2*x))
    xbot = np.hstack((3*x, 4*x))
    xx = np.vstack((xtop, xbot))
    box = np.dstack((xx, xx+(split**2)))

    return box
    

def test_plot(box, box_s):
    """Compares the first and last slices of the original cube (box) and
    the downsampled cube (box_s)."""
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(12, 12))
    axes = axes.ravel().tolist()
    # [ax.axis('off') for ax in axes]
    ims = [axes[0].imshow(box[:, :, 0]), axes[1].imshow(box_s[:, :, 0]),
           axes[2].imshow(box[:, :, -1]), axes[3].imshow(box_s[:, :, -1])]
    axes[0].set_ylabel(r'[:, :, 0]')
    axes[2].set_ylabel(r'[:, :, -1]')
    cbs = [fig.colorbar(im, ax=ax) for im, ax in zip(ims, axes)]
    # fig.savefig('test_grid.pdf')
    # plt.show()

    return fig, ax


def get_regions(d_tol=0.001, n_res=20):
    return
    

def downsample(path, level, field, level_s):
    assert level_s <= level, 'Sampling level must be smaller than or equal to IC level'

    ics = grafic.load_snapshot(path, level, field)
    n = 2**level      # number of cells in original cube
    n_s = 2**level_s  # number of cells in downsampled cube
    s = 2**(level - level_s)  # ratio of original to downsampled cells

    # Storage for the downsampled cube
    box_s = np.zeros((n_s, n_s, n_s), dtype=np.float32)

    box = ics.load_box()

    if (level == level_s):
        print('-- level == level_s, doing no downsampling')
        return box
    
    for ii in range(n_s):
        for jj in range(n_s):
            for kk in range(n_s):
                tmp = np.zeros(s**3)
                
                # Check whether we've looped around
                iis = np.arange(ii*s, ii*s + s)
                jjs = np.arange(jj*s, jj*s + s)
                kks = np.arange(kk*s, kk*s + s)

                l = 0
                for iii in iis:
                    for jjj in jjs:
                        for kkk in kks:
                            tmp[l] = box[iii, jjj, kkk]
                            l += 1
                            
                box_s[ii, jj, kk] = np.mean(tmp)

    return box_s


def get_cen(pp, n):
    """Gets the centre of the chosen cell in code units ([0., 1.]) for
    coordinates pp = [ii, jj, kk] and original grid size n."""

    s = 2**(level - level_s)
    
    return [((p*s + s//2) % n) / n for p in pp]


# path = '/home/lc589/projects/bd/test_ics/test_ics_p18'
# level = 7
# level_s = 4  # level to downsample to, i.e. n_samp = 2**level_s

path = './'
level = 7
level_s = 5

ics = grafic.load_snapshot(path, level, 'deltab')

fac = 1001. * ics.cosmo['aexp']  # scale v_bc to z=1000

d_tol = 0.001  # how close to mean density?
n_res = 20     # how many regions?

# Load up downsampled boxes
box = ics.load_box()
d_box = downsample(path, level, 'deltab', level_s)
v_box = downsample(path, level, 'vbc', level_s)

# Density contrast smaller than some specified value
d_idx = np.argwhere(np.abs(d_box) < d_tol)

print('---- found {0} regions with |\delta| < {1}'.format(d_idx.shape[0], d_tol))

# Find largest v_bc
v_idx = d_idx[np.argsort(v_box[d_idx[:, 0], d_idx[:, 1], d_idx[:, 2]]), :]

h = '-------- v_bc (km/s) | \delta'
s = '         {0:>7.5}     | {1:>7.5}'
print(h)

out = np.zeros((n_res, 5))
j = 0

# Sorts are done in ascending order, so read backwards
for i in range(-1, -1 - n_res, -1):
    v_i = v_box[v_idx[i, 0], v_idx[i, 1], v_idx[i, 2]]
    d_i = d_box[v_idx[i, 0], v_idx[i, 1], v_idx[i, 2]]
    cc = get_cen([v_idx[i, 0], v_idx[i, 1], v_idx[i, 2]], 2**level)

    print(s.format(v_i*fac, d_i))
    
    out[j, 0] = v_i*fac
    out[j, 1] = d_i
    out[j, 2:] = cc
    j += 1
    
np.savetxt('out.txt', out, fmt='%10.5f', header='v_bc,rec (km/s) \delta x y z')
