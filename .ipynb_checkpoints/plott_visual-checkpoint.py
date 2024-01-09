import os
import numpy as np
import plotly.graph_objs as go
import numpy as np
import pylab as pl
import random
import re


import os
import shutil
print(os.getcwd())

from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation


from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1 import ImageGrid

from IPython.display import display, HTML
location_home = os.getcwd()
location_Save = "Save"

location_Save = os.path.join(location_home,location_Save)

backup_list=os.listdir(location_Save)
print(backup_list)


intrestFolder = 'BackupRun_3'
location_backup=os.path.join(location_Save,intrestFolder)
print(location_backup)
files = os.listdir(location_backup)
print(files)
search = "HI_map_"
intrestFiles = []
for i in files:
    if search in i:
        intrestFiles.append(i)
        print(i)


intrestFiles = sorted(intrestFiles, key=lambda s: [int(text) if text.isdigit() else text.lower() for text in re.split('([0-9]+)', s)])[::-1]
print(intrestFiles)

# Load data from the file
filename = intrestFiles[0]  # Assuming intrestFiles is a list of filenames
f = open(os.path.join(location_backup, filename))
temp_mesh = np.fromfile(f, count=3, dtype='int32')
N1, N2, N3 = temp_mesh
data_num = 1
print('No of grids in the density file')
print(N1)
print(N2)
print(N3)
datatype = np.float32
data = np.fromfile(f, dtype=datatype, count=N1 * N2 * N3 * data_num)
f.close()

# Create 3D grid
x, y, z = np.meshgrid(np.arange(N1), np.arange(N2), np.arange(N3))

# Normalize the data
normalized_data = (data - np.min(data)) / (np.max(data) - np.min(data))

# Create a color scale for transparency
opacity_scale = [[value, value] for value in normalized_data.flatten()]

# Create a 3D scatter plot
scatter = go.Scatter3d(
    x=x.flatten(),
    y=y.flatten(),
    z=z.flatten(),
    mode='markers',
    marker=dict(
        size=4,
        color=normalized_data.flatten(),
        colorscale='Viridis',
        opacity=1,  # This opacity is for the entire plot
        colorbar=dict(title='Density')
    ),
    hoverinfo='none',  # Disable hover info for better performance
)

# Layout settings
layout = go.Layout(
    scene=dict(
        xaxis=dict(title='X-axis'),
        yaxis=dict(title='Y-axis'),
        zaxis=dict(title='Z-axis'),
    ),
    title='3D Density Plot',
)

# Create the figure
fig = go.Figure(data=[scatter], layout=layout)

# Show the interactive plot
fig.show()
