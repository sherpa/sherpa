# for evaluation/examples#model_evaluate_example_pha_directly

import os

import numpy as np
import matplotlib.pyplot as plt


# as we jump to the data directory need to know the working directory so
# can write the PNG files there
#
savedir = os.getcwd()


def report(name):
    print("# print({})".format(name))
    print(eval(name))
    print("----------------------------------------")


def dump(name):
    print("# dump")
    print("{}".format(name))
    print(repr(eval(name)))
    print("----------------------------------------")


def savefig(name):
    outfile = os.path.join(savedir, name)
    plt.savefig(outfile)
    print("# Created: {}".format(name))


os.chdir('../../../../sherpa-test-data/sherpatest/')

from sherpa.astro.io import read_arf, read_rmf

arf = read_arf('3c273.arf')
rmf = read_rmf('3c273.rmf')
dump("rmf.detchans")

from sherpa.models.basic import PowLaw1D
mdl = PowLaw1D()

from sherpa.astro.instrument import RSPModelNoPHA
inst = RSPModelNoPHA(arf, rmf, mdl)

dump("inst")
report("inst")

from sherpa.models.model import ArithmeticModel
dump("isinstance(inst, ArithmeticModel)")
dump("inst.pars")

dump("inst(np.arange(1, 1025))")
dump("inst([0.1, 0.2, 0.3])")
dump("inst([0.1, 0.2, 0.3]).size")
dump("inst([10, 20]) == inst([])")

dump("inst([]).sum()")

chans = np.arange(rmf.offset, rmf.offset + rmf.detchans)
ydet = inst(chans)
plt.plot(chans, ydet)
plt.xlabel('Channel')
plt.ylabel('Count / s')
savefig('rspmodelnopha_channel.png')
plt.clf()

report("rmf")
report("arf")

# emid = (rmf.e_min + rmf.e_max) / 2
# de = rmf.e_max - rmf.e_min
# plt.plot(emid, arf.exposure * ydet / de)
# plt.bar(rmf.e_min, arf.exposure * ydet / de, width=de, align='edge',
#         log=True)

x = np.vstack((rmf.e_min, rmf.e_max)).T.flatten()
y = arf.exposure * ydet / (rmf.e_max - rmf.e_min)
y = y.repeat(2)
plt.plot(x, y, '-')

plt.yscale('log')
plt.ylim(1e4, 3e6)
plt.xlim(0, 10)
plt.xlabel('Energy (keV)')
plt.ylabel('Count / keV')
savefig('rspmodelnopha_energy.png')


from sherpa.astro.io import read_pha
from sherpa.astro.instrument import RSPModelPHA

pha2 = read_pha('3c273.pi')
arf2 = pha2.get_arf()
rmf2 = pha2.get_rmf()

mdl2 = PowLaw1D('mdl2')
inst2 = RSPModelPHA(arf2, rmf2, pha2, mdl2)
report("inst2")

dump("inst2([]).size")

pha2.set_analysis('energy')
report('pha2.get_filter()')
report('pha2.get_filter_expr()')

pha2.notice(0.5, 7.0)
report('pha2.get_filter()')
report('pha2.get_filter_expr()')

dump("pha2.grouped")
# pha2.ungroup()
# report('pha2.get_filter_expr()')

pha2.ignore(2.0, 3.0)
report('pha2.get_filter_expr()')

y1 = inst2([])
inst2.startup()
y2 = inst2([])
inst2.teardown()
report("y1.size, y2.size")
report("np.all(y1 == y2)")

plt.clf()
plt.plot(pha2.channel, y1, label='all')
plt.plot(pha2.channel, y2, label='filtered')
plt.xscale('log')
plt.yscale('log')
plt.ylim(0.001, 1)
plt.xlim(5, 1000)
plt.legend(loc='center')
savefig('rspmodelpha_compare.png')
