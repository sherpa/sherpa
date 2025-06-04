
sources = [
  'grplib.c',
  'grp_priv.c',
  'stub_dserr.c',
  'stub_dserr.h',
  'grp_priv.h',
  ]

grp = static_library(
  'grp',
  sources,
  c_args: '-DSTUB_ERROR_LIB=1',
  install: false,
)

grp_dep = declare_dependency(
  link_with: grp,
  include_directories: include_directories('.')
)  
