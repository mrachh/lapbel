from numpy import *
from pylab import *


norder = 3
npat = 192
npts = npat*10
skiprow1 = 11
maxrow1 = 1
fname = 'ressum-'+str(npat)+'-'+"{:02d}".format(norder)+'.dat'
niters = loadtxt(fname,skiprows=skiprow1,max_rows=maxrow1)
niters = niters.astype(int)
niterall = sum(niters)

maxrow2 = niterall
skiprow2 = skiprow1+3
errs192 = loadtxt(fname,skiprows=skiprow2,max_rows=maxrow2)

skiprow3 = skiprow2+maxrow2+2
svals192 = loadtxt(fname,skiprows=skiprow3)

xx192 = linspace(0,1,npts)

npat = 768
npts = npat*10
skiprow1 = 11
maxrow1 = 1
fname = 'ressum-'+str(npat)+'-'+"{:02d}".format(norder)+'.dat'
niters = loadtxt(fname,skiprows=skiprow1,max_rows=maxrow1)
niters = niters.astype(int)
niterall = sum(niters)

maxrow2 = niterall
skiprow2 = skiprow1+3
errs768 = loadtxt(fname,skiprows=skiprow2,max_rows=maxrow2)

skiprow3 = skiprow2+maxrow2+2
svals768 = loadtxt(fname,skiprows=skiprow3)

xx768 = linspace(0,1,npts)


markers = ['*','^','s','x']
cc = ['darkgreen','navy','orangered','magenta']



figure(1)
clf()
semilogy(errs192[0:niters[0]],'-*',color=cc[0],markersize=5)
semilogy(errs192[niters[0]:-1],'-^',color=cc[1],markersize=5)
semilogy(errs768[0:niters[0]],'-s',color=cc[2],markersize=3)
semilogy(errs768[niters[0]:-1],'-x',color=cc[3],markersize=5)
tuse1 = [70, 1e-10]
tuse2 = [70, 1e-11]
tuse3 = [70, 1e-12]
tuse4 = [70, 1e-13]

tt = [68,72]
semilogy(tt,ones(2)*tuse1[1],'-',color=cc[0],linewidth=1,alpha=0.8)
semilogy(tuse1[0],tuse1[1],'*',color=cc[0],markersize=5)
semilogy(tt,ones(2)*tuse2[1],'-',color=cc[1],linewidth=1,alpha=0.8)
semilogy(tuse2[0],tuse2[1],'^',color=cc[1],markersize=5)
semilogy(tt,ones(2)*tuse3[1],'-',color=cc[2],linewidth=1,alpha=0.8)
semilogy(tuse3[0],tuse3[1],'s',color=cc[2],markersize=3)
semilogy(tt,ones(2)*tuse4[1],'-',color=cc[3],linewidth=1,alpha=0.8)
semilogy(tuse4[0],tuse4[1],'x',color=cc[3],markersize=5)
savefig('gmres_iteration.pdf',bbox_inches='tight')



figure(2)
clf()
semilogy(xx192,svals192[:,0],'*',color=cc[0],markersize=5,alpha=0.7)
semilogy(xx192,svals192[:,1],'^',color=cc[1],markersize=5,alpha=0.5)
semilogy(xx768,svals768[:,0],'s',color=cc[2],markersize=3,alpha=0.5)
semilogy(xx768,svals768[:,1],'x',color=cc[3],markersize=5,alpha=0.5)

tuse1 = [0.05, 1e-4]
tuse2 = [0.05, 3e-5]
tuse3 = [0.05, 1e-5]
tuse4 = [0.05, 3e-6]
semilogy(tuse1[0],tuse1[1],'*',color=cc[0],markersize=5)
semilogy(tuse2[0],tuse2[1],'^',color=cc[1],markersize=5)
semilogy(tuse3[0],tuse3[1],'s',color=cc[2],markersize=3)
semilogy(tuse4[0],tuse4[1],'x',color=cc[3],markersize=5)
tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False) # labels along the bottom edge are off
savefig('svd.pdf',bbox_inches='tight')
show()
