import utilities.io as util
import numpy
import matplotlib.pyplot as plt
x = numpy.linspace(0, 10, 10)
y = numpy.cos(x)

# f = open('asdf', 'w')
# for i in range(len(x)):
#     f.write('%f %f \n' % (x[i], y[i]))
# f.close()
stl1 = ['3000', '5000', '7000']
stl2 = ['tran_ImEx', 'tran_ReEx', 'refl_ImEx', 'refl_ReEx']

for st1 in stl1:
    for st2 in stl2:
        util.rebinFile('../!TMM_Reconstruction/02/'+st1+'/'+st2+'_'+st1+'.txt', 'data/'+st2+'_'+st1+'.txt', 400, 600, 2)
        print('rebinned ../!TMM_Reconstruction/02/'+st1+'/'+st2+'_'+st1+'.txt')
# nx, ny = util.rebin(x, y, 2, 8, 0.01)
# plt.plot(x,y,'--',nx,ny,'-')
# plt.show(block=True)