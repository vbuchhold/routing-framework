# Directories containing .cpp source files.
src_dirs = ['RawData']

# Common compiler options.
cflags   = '-std=c11'
cxxflags = '-std=c++11'
ccflags  = '-Werror -Wfatal-errors -Wall -Wextra '

# Additional compiler options per build variant.
debug   = '-O0 -g'
devel   = '-O3'
release = '-O3 -DNDEBUG'

# Get build variant from the command line (defaults to 'Devel').
variant = ARGUMENTS.get('variant', 'Devel')

# Assemble ccflags according to the build variant.
if variant == 'Debug':
  ccflags += debug
elif variant == 'Devel':
  ccflags += devel
elif variant == 'Release':
  ccflags += release
else:
  print 'Invalid value for option variant: ' + variant + '.'
  print "Valid values are: ('Debug', 'Devel', 'Release')"
  Exit(1)

# Forward macro definitions to the compiler.
cppdefines = []
for keyword, value in ARGLIST:
  if keyword == 'def':
    cppdefines.append(value.split('=', 1))

# Set construction variables in the default environment.
DefaultEnvironment(
    CFLAGS     = cflags,
    CXXFLAGS   = cxxflags,
    CCFLAGS    = ccflags,
    CPPDEFINES = cppdefines,
    CPPPATH    = '#')

# Execute subsidiary SConscript files.
for dir in src_dirs:
  targets = SConscript(dirs = dir, variant_dir = 'Build/' + variant + '/' + dir, duplicate = 0)
  Alias(dir.split('/')[-1], targets)

