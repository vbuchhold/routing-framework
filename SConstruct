import os

# Directories containing .cpp source files.
srcDirs = ['Launchers', 'RawData']

# Common compiler options.
cflags = '-std=c99'
cxxflags = '-std=c++14'
ccflags = (
  '-Werror -Wfatal-errors -Wpedantic -pedantic-errors -Wall -Wextra '
  '-Wno-missing-braces -Wno-unknown-pragmas -Wno-strict-overflow -Wno-sign-compare ')

# Additional compiler options per build variant.
debug = '-O0 -g '
devel = '-O3 '
release = '-O3 -DNDEBUG '

# Additional compiler options per instruction set extension.
sse = '-msse4'
avx = '-mavx2 -mfma'

# Get options from the command line.
variant = ARGUMENTS.get('variant', 'Devel')
simd = ARGUMENTS.get('simd', 'sse')

# Assemble ccflags according to the build variant.
if variant == 'Debug':
  ccflags += debug
elif variant == 'Devel':
  ccflags += devel
elif variant == 'Release':
  ccflags += release
else:
  print('Invalid value for option variant: ' + variant + '.')
  print("Valid values are: ('Debug', 'Devel', 'Release')")
  Exit(1)

# Assemble ccflags according to the SIMD extension.
if simd == 'sse':
  ccflags += sse
elif simd == 'avx':
  ccflags += avx
else:
  print('Invalid value for option simd: ' + simd + '.')
  print("Valid values are: ('sse', 'avx')")
  Exit(1)

# Forward macro definitions to the compiler.
cppdefines = []
for keyword, value in ARGLIST:
  if keyword == 'def':
    cppdefines.append(value.split('=', 1))

# Set construction variables in the default environment.
DefaultEnvironment(
  CFLAGS = cflags,
  CXXFLAGS = cxxflags,
  CCFLAGS = ccflags,
  CPPDEFINES = cppdefines,
  CPPPATH = '#',
  ENV = {
    'CPLUS_INCLUDE_PATH': os.environ.get('CPLUS_INCLUDE_PATH'),
    'LIBRARY_PATH': os.environ.get('LIBRARY_PATH'),
    'PATH': os.environ.get('PATH')})

# Execute subsidiary SConscript files.
for dir in srcDirs:
  targets = SConscript(dirs=dir, variant_dir='Build/' + variant + '/' + dir, duplicate=0)
  Alias(dir.split('/')[-1], targets)
