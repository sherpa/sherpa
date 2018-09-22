from sherpa.ui.utils import Session
import sherpa.models.basic

x = [100, 200, 300, 400]
y = [10, 12, 9, 13]

s = Session()
s._add_model_types(sherpa.models.basic)

s.load_arrays(1, x, y)

print("# list_data_ids")
print(s.list_data_ids())

print("# get_data()")
print(repr(s.get_data()))

print("# get_data()")
print(s.get_data())

print("# get_stat_name/get_method_name")
print(s.get_stat_name())
print(s.get_method_name())

s.set_stat('cash')
s.set_method('simplex')

s.set_source('const1d.mdl')

print("# mdl")
print(mdl)

print("# get_source")
print(s.get_source())

print("# fit")
s.fit()

print("# get_fit_results")
r = s.get_fit_results()
print(r)

s.get_data_plot_prefs()['yerrorbars'] = False
print("--> call")
print("s.plot_fit()")
