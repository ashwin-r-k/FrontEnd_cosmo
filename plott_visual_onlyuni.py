import re
import os
import pylab as pl
import numpy as np

from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation

#plt.rcParams["figure.figsize"] = [7.50, 3.50]
#plt.rcParams["figure.autolayout"] = True
#plt.rcParams['animation.ffmpeg_path']='/gpfs-home/m220590ph/softwares/ffmpeg/bin/ffmpeg'

print(os.getcwd())
os.chdir("../")
location_home = os.getcwd()
location_Save = "Save"
location_Save = os.path.join(location_home,location_Save)

print(location_Save)

#latest
#location_backup = max([os.path.join(location_Save,d) for d in os.listdir(location_Save)], key=os.path.getmtime)
#intrestFolder = os.path.basename(location_backup)
backup_list=os.listdir(location_Save)
search = "BackupRun"
intrestFolders = []
for i in backup_list:
    if search in i:
        intrestFolders.append(i)

intrestFolder_latest = sorted(intrestFolders, key=lambda s: [int(text) if text.isdigit() else text.lower() for text in re.split('([0-9]+)', s)])[::-1]
location_backup = intrestFolder_latest[0]
intrestFolder = os.path.basename(location_backup)
location_backup = os.path.join(location_Save,intrestFolder)

files = os.listdir(location_backup)

search = "HI_map_"
intrestFiles = []
for i in files:
    if search in i:
        intrestFiles.append(i)
        #print(i)
intrestFiles = sorted(intrestFiles, key=lambda s: [int(text) if text.isdigit() else text.lower() for text in re.split('([0-9]+)', s)])[::-1]
print(intrestFiles)


filename=intrestFiles[0]
f = open( os.path.join(location_backup,filename))
temp_mesh = np.fromfile(f,count=3,dtype='int32')
N1,N2,N3 = temp_mesh
data_num =1
print('No of grids in the desnity file')
print(N1)
print(N2)
print(N3)
datatype = dtype=(np.float32)
data = np.fromfile(f, dtype=datatype,count=N1*N2*N3*data_num)
f.close()

#data = np.loadtxt(filename)
#print(data.shape)
dT = data.reshape((N1,N2,N3,data_num), order='C') # Reshaping the data into N1xN2xN3 grid format
 
#print(dT[0,0,9,0])
Tpl=dT[:,:,:,0]  # extracting out the first 5 slices
#print(np.shape(Tpl))
av_Tpl=np.mean(Tpl, axis=0) # stacking (averaging) the 5 slices along x-axis(axis=0) 
#print(np.shape(av_Tpl))

#fig = pl.figure(1, (5., 5.))
#pl.imshow(av_Tpl,origin='lower')
#pl.colorbar(pad=0.01,fraction=0.047)
#pl.clim(0,5)

#pl.save()

fig, ax = plt.subplots(figsize=(10.24,10.24))
plt.axis('off')
fig.tight_layout(pad=0,h_pad=None, w_pad=None)
im = ax.imshow(av_Tpl, aspect="equal", cmap='viridis')
#cbar = fig.colorbar(im, pad=0.01, fraction=0.047)



def update(frame):
    filename = intrestFiles[frame]
    f = open(os.path.join(location_backup, filename), 'rb')
    temp_mesh = np.fromfile(f, count=3, dtype='int32')
    N1, N2, N3 = temp_mesh
    data = np.fromfile(f, dtype=datatype, count=N1*N2*N3*data_num)
    f.close()

    dT = data.reshape((N1, N2, N3, data_num), order='C')
    Tpl = dT[:, :, :, 0]
    av_Tpl = np.mean(Tpl, axis=0)


  #  pl.imshow(av_Tpl,origin='lower')
    im.set_array(av_Tpl)

#    z=re.search(r'([0-9]+)', filename).group()
#    ax.set_title(f'Z= {z}')
    return im,

anim = FuncAnimation(fig, update, frames=len(intrestFiles), interval=1000, blit=True)

#animation.save(f'Save/animation/mov_{intrestFolder}_{search}.mp4', writer='ffmpeg', fps=1)

f = f'Save/animation/mov_{intrestFolder}_{search}_clean.mp4'
writervideo = animation.FFMpegWriter(fps=5) 
anim.save(f, writer=writervideo)

f = f'Save/animation/mov_{intrestFolder}_{search}_clean.gif'
writergif = animation.PillowWriter(fps=30) 
anim.save(f, writer=writergif)


print("Animination saved as "+f'Save/animation/mov_{intrestFolder}_{search}_clean.mp4')
