
sources = [
  'pystk.c',
  'pystk.h',
  'logical.h',
  'stack.h',
  'stack.c'
  ]

stk = py.extension_module(
  'stk',
  sources,
  install: true
)

stklib_dep = declare_dependency(link_with: stk)
