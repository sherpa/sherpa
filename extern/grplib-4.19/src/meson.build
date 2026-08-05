
sources = [
  'grplib.c',
  'grp_priv.c',
  'stub_dserr.c',
  'stub_dserr.h',
  'grp_priv.h',
  ]

extra_args = [
  '-DSTUB_ERROR_LIB=1',
  '-Wno-unused-but-set-variable'
]

grp = static_library(
  'grp',
  sources,
  c_args: extra_args,
  install: false,
)

grp_dep = declare_dependency(
  link_with: grp,
  include_directories: include_directories('.')
)
