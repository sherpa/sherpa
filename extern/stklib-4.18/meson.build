# stklib

project(
  'stklib',
  ['c'],
  version: '4.11',

  meson_version: '>=1.3.0',
  default_options: [
    'c_std=c11',
    'python.allow_limited_api=false'
  ]
)

# If a limited API is used, what is the version to use (this needs to
# be updated to match the minimum-supported Python version in pyproject.toml).
#
limited_api_version = '3.11'

py = import('python').find_installation(pure: false)

subdir('src')
