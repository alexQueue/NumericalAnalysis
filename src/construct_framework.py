import numpy as np
import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(0,0)
ax.set(xlim=(0,10),ylim=(0,10))
ax.grid()

coords = []

def onclick(event):
    ix, iy = event.xdata, event.ydata

    if ix is None or iy is None:
        return
    coords.append((ix, iy))

    ax.plot(ix,iy,marker="o",color="black")
    ax.set(xlim=(0,10),ylim=(0,10))
    fig.canvas.draw()

cid = fig.canvas.mpl_connect('button_press_event', onclick)

plt.show()


def find_closest(coords, point):
    min_length = None
    index = None
    for i,coord in enumerate(coords):
        length = np.linalg.norm(np.array(coord)-np.array(point))
        if min_length is None or length < min_length:
            min_length = length
            index = i
    return index

fig = plt.figure()
ax = fig.add_subplot(111)
x = [x[0] for x in coords]
y = [y[1] for y in coords]
ax.plot(x,y,"o",color="black")
ax.set(xlim=(0,10),ylim=(0,10))
for i,(x,y) in enumerate(coords):
    ax.annotate(i+1, (x+0.1,y+0.1))

edges = []
line_coords = []
line_indices = []
node_indices = []

def onclick(event):
    ix, iy = event.xdata, event.ydata
    try:
        node_index = find_closest(coords, (ix,iy))
    except TypeError:
        return
    node_indices.append(node_index + 1)
    line_coords.append(coords[node_index])
    if len(line_coords) == 2:
        if sorted([x for x in node_indices]) not in [sorted(x) for x in edges]:
            ax.plot([line_coords[0][0],line_coords[1][0]],[line_coords[0][1],line_coords[1][1]],color="black")
            ax.set(xlim=(0,10),ylim=(0,10))
            edges.append((node_indices[0],node_indices[1]))
        
        line_coords.pop()
        line_coords.pop()
        node_indices.pop()
        node_indices.pop()

    fig.canvas.draw()

cid = fig.canvas.mpl_connect('button_press_event', onclick)
plt.show(block=False)

print('Determine type of each node:\n 1 -> FIXED\n 2 -> FREE\n 3 -> FORCE\n 4 -> MOVABLE')
allowed_vals = ['1','2','3','4']
rest = {'1':'FIXED','2':'FREE'}
types = []
for i in range(1,len(coords)+1):
    nr = ''
    while nr not in allowed_vals:
        nr = input(f'Node {i}: ')
    if nr == '3':
        Fx = input('F_x = ')
        Fy = input('F_y = ')
        M = input('M = ')
        types.append(f'{i} FORCE [{Fx} {Fy}] [{M}]\n')
    elif nr == '4':
        x = input('x-val = ')
        y = input('y-val = ')
        types.append(f'{i} MOVABLE [{x} {y}]\n')
    else:
        types.append(f'{i} {rest[nr]}\n')


parameters = ''
print('Determine parameter functions for E,I,A and mu')
pars = ['E','I','A','mu']
for par in pars:
    fnc = input(f'{par}(x) = ')
    parameters += f'{par} {fnc}\n'
parameters = parameters[0:-1]

save_file = input('Save to file (.ne will be added): ')

with open(f'../data/{save_file}.ne', 'w') as f:
    f.write('NODES\n')
    for coord in coords:
        f.write(f'{coord[0]} {coord[1]}\n')
    f.write('EDGES\n')
    for edge in edges:
        f.write(f'{edge[0]} {edge[1]}\n')
    f.write('TYPE\n')
    for type_ in types:
        f.write(type_)
    f.write('PARAMETERS\n')
    f.write(parameters)