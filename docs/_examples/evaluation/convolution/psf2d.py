# for evaluation/convolution.rst
# save in docs/_static/models/combine

from matplotlib import pyplot as plt
import numpy as np

from sherpa.models.basic import Box2D
from sherpa.instrument import PSFModel
from sherpa.data import Data2D

k = np.asarray([[0, 1, 0], [1, 0, 1], [0, 1, 0]])
yg, xg = np.mgrid[:3, :3]
kernel = Data2D('kdata', xg.flatten(), yg.flatten(), k.flatten(),
                shape=k.shape)
psf = PSFModel(kernel=kernel)

print("# print(psf)")
print(psf)

pt = Box2D('pt')
box = Box2D('box')
unconvolved_mdl = pt + box

pt.xlow = 1.5
pt.xhi = 2.5
pt.ylow = 2.5
pt.yhi = 3.5
pt.ampl = 10

box.xlow = 4
box.xhi = 10
box.ylow = 6.5
box.yhi = 7.5
box.ampl = 10
print("# print(unconvolved_mdl)")
print(unconvolved_mdl)

convolved_mdl = psf(unconvolved_mdl)
print("# print(convolved_mdl)")
print(convolved_mdl)

from sherpa.models.basic import Const2D
bgnd = Const2D('bgnd')
bgnd.c0 = 0.25

final_mdl = convolved_mdl + bgnd
print("# print(final_mdl)")
print(final_mdl)

yg, xg = np.mgrid[:10, :10]
xg1d, yg1d = xg.flatten(), yg.flatten()
m1 = unconvolved_mdl(xg1d, yg1d).reshape(xg.shape)

blank = Data2D('blank', xg1d, yg1d, np.ones(xg1d.shape), xg.shape)
m1 = blank.eval_model(unconvolved_mdl).reshape(xg.shape)

psf.fold(blank)
m2 = blank.eval_model(convolved_mdl).reshape(xg.shape)

# TODO: set the axes? Actually, be design they are the same as
#       that shown by matplotlib
#
plt.imshow(m1, origin='lower', cmap='viridis')
plt.colorbar()
plt.title('Unconvolved')
plt.savefig('convolution_psf2d_evaluate_unconv.png')
print('Created: convolution_psf2d_evaluate_unconv.png')
plt.clf()

plt.imshow(m2, origin='lower', cmap='viridis', vmin=0, vmax=10)
plt.colorbar()
plt.title('Convolved')
plt.savefig('convolution_psf2d_evaluate_conv.png')
print('Created: convolution_psf2d_evaluate_conv.png')
plt.clf()

print("Sum: m1 = {}".format(m1.sum()))
print("Sum: m2 = {}".format(m2.sum()))

print("m2[7] =\n{}\n".format(repr(m2[7])))
