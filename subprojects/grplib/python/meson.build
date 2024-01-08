
sources = [
  'pygrplib.c',
  'pygrplib.h',
  ]

includes = [
  incdir_numpy,
  include_directories('.'),
  ]

grplib = py.extension_module(
  'group',
  sources,
  c_args: ['-DSTUB_ERROR_LIB=1'],
  dependencies: grp_dep,
  include_directories: includes,
  install: true,
)

grplib_dep = declare_dependency(link_with: grplib)
