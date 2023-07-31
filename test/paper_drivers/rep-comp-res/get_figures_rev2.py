from numpy import *
from pylab import *


norder = 5 
npat = 192
npol = (norder+1)*(norder+2)/2
npts = int(npat*npol)
skiprow1 = 11
maxrow1 = 1
fname = 'ressum-'+str(npat)+'-norder-'+"{:02d}".format(norder)+'.dat1'
fname2 = 'ssing_npat'+str(npat)+'_norder'+"{0:02d}".format(norder)+".txt"
niters = loadtxt(fname,skiprows=skiprow1,max_rows=maxrow1)
niters = niters.astype(int)
niterall = sum(niters)

maxrow2 = niterall
skiprow2 = skiprow1+3
errs192 = loadtxt(fname,skiprows=skiprow2,max_rows=maxrow2)

svals192 = loadtxt(fname2,delimiter=',')

xx192 = linspace(0,1,npts)

npat = 768
npts = int(npat*npol)
skiprow1 = 11
maxrow1 = 1
fname = 'ressum-'+str(npat)+'-norder-'+"{:02d}".format(norder)+'.dat1'
fname2 = 'ssing_npat'+str(npat)+'_norder'+"{0:02d}".format(norder)+".txt"
niters768 = loadtxt(fname,skiprows=skiprow1,max_rows=maxrow1)
niters768 = niters768.astype(int)
niterall768 = sum(niters768)

maxrow2 = niterall768
skiprow2 = skiprow1+3
errs768 = loadtxt(fname,skiprows=skiprow2,max_rows=maxrow2)

svals768 = loadtxt(fname2,delimiter=',')

xx768 = linspace(0,1,npts)


markers = ['*','^','s','x']
cc = ['darkgreen','navy','orangered','magenta']



figure(1)
clf()
semilogy(errs192[0:niters[0]],'-*',color=cc[0],markersize=5)
semilogy(errs192[niters[0]:(niters[0]+niters[1])],'-*',color=cc[1],markersize=5)
semilogy(errs192[(niters[0]+niters[1])::],'-*',color=cc[2],markersize=5)

semilogy(errs768[0:niters768[0]],'-^',color=cc[0],markersize=5)
semilogy(errs768[niters768[0]:(niters768[0]+niters768[1])],'-^',color=cc[1],markersize=5)
semilogy(errs768[(niters768[0]+niters768[1])::],'-^',color=cc[2],markersize=5)
tuse1 = [70, 1e-11]
tuse2 = [70, 1e-12]

tuse3 = [120, 5e-11]
tuse4 = [120, 5e-12]
tuse5 = [120, 5e-13]

tt = [68,72]
tt2 = [119,121]
semilogy(tt,ones(2)*tuse1[1],'-',color='black',linewidth=1,alpha=0.8)
semilogy(tuse1[0],tuse1[1],'*',color='black',markersize=5)
semilogy(tt,ones(2)*tuse2[1],'-',color='black',linewidth=1,alpha=0.8)
semilogy(tuse2[0],tuse2[1],'^',color='black',markersize=5)

semilogy(tt2,ones(2)*tuse3[1],'-',color=cc[0],linewidth=2,alpha=0.8)
semilogy(tt2,ones(2)*tuse4[1],'-',color=cc[1],linewidth=2,alpha=0.8)
semilogy(tt2,ones(2)*tuse5[1],'-',color=cc[2],linewidth=2,alpha=0.8)
savefig('gmres_iteration_rev'+str(norder)+'.pdf',bbox_inches='tight')


figure(2)
clf()
semilogy(xx192,svals192[:,0],'*',color=cc[0],markersize=5,alpha=0.7)
semilogy(xx192,svals192[:,1],'*',color=cc[1],markersize=5,alpha=0.5)
semilogy(xx192,svals192[:,2],'*',color=cc[2],markersize=5,alpha=0.5)

semilogy(xx768,svals768[:,0],'^',color=cc[0],markersize=5,alpha=0.7)
semilogy(xx768,svals768[:,1],'^',color=cc[1],markersize=5,alpha=0.5)
semilogy(xx768,svals768[:,2],'^',color=cc[2],markersize=5,alpha=0.5)

tuse1 = [0.05, 3e-7]
tuse2 = [0.05, 1e-7]

tt2 = [0.31,0.32]
tuse3 = [0.315,5.5e-7]
tuse4 = [0.315,1.8e-7]
tuse5 = [0.315,5.5e-8]

semilogy(tuse1[0],tuse1[1],'*',color='black',markersize=5)
semilogy(tuse2[0],tuse2[1],'^',color='black',markersize=5)

semilogy(tt2,ones(2)*tuse3[1],'-',color=cc[0],linewidth=2,alpha=0.8)
semilogy(tt2,ones(2)*tuse4[1],'-',color=cc[1],linewidth=2,alpha=0.8)
semilogy(tt2,ones(2)*tuse5[1],'-',color=cc[2],linewidth=2,alpha=0.8)
ylim([1e-8,2])
tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False) # labels along the bottom edge are off
savefig('svd_rev'+str(norder)+'.pdf',bbox_inches='tight')
show()
