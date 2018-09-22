from sherpa.ui import *

x = [100, 200, 300, 400]
y = [10, 12, 9, 13]

load_arrays(1, x, y)

print("# list_data_ids")
print(list_data_ids())

print("# get_data()")
print(repr(get_data()))

print("# get_data()")
print(get_data())

print("# get_stat_name/get_method_name")
print(get_stat_name())
print(get_method_name())

set_stat('cash')
set_method('simplex')

set_source('const1d.mdl')

print("# mdl")
print(mdl)

print("# get_source")
print(get_source())

print("# fit")
fit()

print("# get_fit_results")
r = get_fit_results()
print(r)

get_data_plot_prefs()['yerrorbars'] = False
print("--> call")
print("plot_fit()")
