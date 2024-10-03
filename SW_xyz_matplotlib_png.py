import glob
import matplotlib.pyplot as plt
import torch
import sys
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import imageio.v2 as imageio
import bond_analyze


dataset_info={
    'atom_encoder': {'H': 0, 'C': 1, 'N': 2, 'O': 3, 'F': 4},
    'colors_dic': ['#FFFFFF99', 'C7', 'C0', 'C3', 'C1'],
    'radius_dic': [0.46, 0.77, 0.77, 0.77, 0.77]
}


dataset_info={
     'atom_encoder': {'H': 0, 'B': 1, 'C': 2, 'N': 3, 'O': 4, 'F': 5, 'Al': 6, 'Si': 7,
    'P': 8, 'S': 9, 'Cl': 10, 'As': 11, 'Br': 12, 'I': 13, 'Hg': 14, 'Bi': 15},
    'atom_decoder': ['H', 'B', 'C', 'N', 'O', 'F', 'Al', 'Si', 'P', 'S', 'Cl', 'As', 'Br', 'I', 'Hg', 'Bi'],
    'colors_dic': ['#FFFFFF99','C2', 'C7', 'C0', 'C3', 'C1', 'C5','C6', 'C4', 'C8', 'C9', 'C10','C11', 'C12','C13', 'C14'],
    'radius_dic': [0.3, 0.6, 0.6, 0.6, 0.6,0.6, 0.6, 0.6, 0.6, 0.6,0.6, 0.6, 0.6, 0.6, 0.6,0.6]
}





def draw_sphere(ax, x, y, z, size, color, alpha):
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    xs = size * np.outer(np.cos(u), np.sin(v))
    ys = size * np.outer(np.sin(u), np.sin(v)) * 0.8  # Correct for matplotlib.
    zs = size * np.outer(np.ones(np.size(u)), np.cos(v))
    ax.plot_surface(x + xs, y + ys, z + zs, rstride=2, cstride=2, color=color, linewidth=0,
                    alpha=alpha)

def plot_molecule(ax, positions, atom_type,dataset_info, alpha, spheres_3d, hex_bg_color):
    x = positions[:, 0]
    y = positions[:, 1]
    z = positions[:, 2]
    colors_dic = np.array(dataset_info['colors_dic'])
    radius_dic = np.array(dataset_info['radius_dic'])
    area_dic = 1500 * radius_dic ** 2
    areas = [area_dic[int(i)] for i in atom_type]
    radii = [radius_dic[int(i)] for i in atom_type]
    colors = [colors_dic[int(i)] for i in atom_type]
    if spheres_3d:
        for i, j, k, s, c in zip(x, y, z, radii, colors):
            draw_sphere(ax, i.item(), j.item(), k.item(), 0.7 * s, c, alpha)
    else:
        ax.scatter(x, y, z, s=areas, alpha=0.9 * alpha,c=colors) 
    for i in range(len(x)):
        for j in range(i + 1, len(x)):
            p1 = np.array([x[i], y[i], z[i]])
            p2 = np.array([x[j], y[j], z[j]])
            dist = np.sqrt(np.sum((p1 - p2) ** 2))
            s = sorted((atom_type[i], atom_type[j]))
            pair = (dataset_info['atom_decoder'][int(s[0])],
                    dataset_info['atom_decoder'][int(s[1])])
            draw_edge_int = bond_analyze.geom_predictor(pair, dist)
            line_width = 2
            draw_edge = draw_edge_int > 0
            if draw_edge:
                if draw_edge_int == 4:
                    linewidth_factor = 1.5
                else:
                    linewidth_factor = 1
                ax.plot([x[i], x[j]], [y[i], y[j]], [z[i], z[j]],linewidth=line_width*linewidth_factor,c=hex_bg_color, alpha=alpha)

def plot_data3d(positions, atom_type,dataset_info,camera_elev=133, camera_azim=121,save_path=None,spheres_3d=False,bg='black',alpha=1.):
    black = (0, 0, 0)
    white = (1, 1, 1)
    hex_bg_color = '#FFFFFF' if bg == 'black' else '#666666'
    fig=plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.set_aspect('auto')
    ax.view_init(elev=camera_elev, azim=camera_azim)
    if bg == 'black':
        ax.set_facecolor(black)
    else:
        ax.set_facecolor(white)
        ax.xaxis.pane.set_alpha(0)
    ax.yaxis.pane.set_alpha(0)
    ax.zaxis.pane.set_alpha(0)
    ax._axis3don = False
    if bg == 'black':
        ax.xaxis.line.set_color("black")
    else:
        ax.xaxis.line.set_color("white")
    plot_molecule(ax, positions, atom_type,dataset_info, alpha, spheres_3d,hex_bg_color)
    # set x,y,z axis limit
    max_value = positions.abs().max().item()
    axis_lim = min(40, max(max_value / 1.5 + 0.3, 3.2))
    ax.set_xlim(-axis_lim, axis_lim)
    ax.set_ylim(-axis_lim, axis_lim)
    ax.set_zlim(-axis_lim, axis_lim)
    dpi = 120 if spheres_3d else 50
    if save_path is not None:
        plt.savefig(save_path, bbox_inches='tight', pad_inches=0.0, dpi=dpi)
        if spheres_3d:
            img = imageio.imread(save_path)
            img_brighter = np.clip(img * 1.4, 0, 255).astype('uint8')
            imageio.imsave(save_path, img_brighter)
    else:
        plt.show()
    plt.close()


folder = './xyz_example'
files = glob.glob(folder+'/*txt')

with open(files[0],'r') as f:
    data=f.readlines()
    n_atoms = int(data[0])
    extract_atom=data[2:]
    atom_type=torch.zeros(n_atoms,1)
    positions=torch.zeros(n_atoms,3)
    for item in range(n_atoms):
        atom=extract_atom[item].strip().split('    ')
        atom_type[item]=dataset_info['atom_encoder'][atom[0]]
        position=torch.Tensor([float(i) for i in atom[1:]])
        positions[item]=position
    atom_type=atom_type.squeeze().numpy()

plot_data3d(positions, atom_type,dataset_info,camera_elev=102, camera_azim=130,save_path="test.png",spheres_3d=True,bg='black',alpha=1.)