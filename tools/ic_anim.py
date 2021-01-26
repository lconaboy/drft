def main(path, level, field, cmap='viridis', vlim=None):
    import numpy as np
    import grafic_tools as grafic
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation

    def update(i): 
        im.set_data(x[:, :, i])

    # Load data
    x = grafic.load_snapshot(path, level, field).load_box()
    
    if vlim is None:
        vmin = x.min()
        vmax = x.max()
    else:
        vmin = vlim[0]
        vmax = vlim[1]

    # Set up figure
    fig, ax = plt.subplots(figsize=(6, 6))
    # plt.axis('off')
    
    # Remove whitespace
    # (https://stackoverflow.com/questions/15882395/matplotlib-animation-how-to-remove-white-margin)
    fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=None, hspace=None)

    # Call imshow now, then we can just update the data later
    im = ax.imshow(np.zeros_like(x[:, :, 0]), vmin=vmin, vmax=vmax, cmap=cmap)

    cax = fig.add_axes([0.2, 0.9, 0.6, 0.03])
    cb = fig.colorbar(im, cax=cax,orientation='horizontal')
    cb.ax.title.set_color('k')
    # cbar.ax.tick_params(axis='x',colors='Black')
    cb.ax.set_title(field)

    
    # This is the main work
    anim = FuncAnimation(fig, update, frames=np.arange(0, x.shape[2]), interval=200)
    anim.save(field + '.gif', dpi=100, writer='imagemagick')

    
if __name__ == '__main__':
    import sys
    
    if len(sys.argv) < 4:
        print('Usage: python ic_anim.py <path> <level> <field> [<cmap> <vmin> <vmax>]')
        sys.exit()
        
    path = sys.argv[1]
    level = int(sys.argv[2])
    field = sys.argv[3]
    cmap = 'viridis'
    vlim = None

    if len(sys.argv) > 4:
        cmap = sys.argv[4]
    
    if len(sys.argv) > 5:
        vlim = [float(sys.argv[5]), float(sys.argv[6])]

    main(path, level, field, cmap, vlim)
