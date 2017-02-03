# Routing Framework

## Prerequisites

You need to have some tools and libraries installed. On Debian and its derivatives (such as Ubuntu)
the `apt-get` tool can be used:

```
$ sudo apt-get install build-essential
$ sudo apt-get install scons
$ sudo apt-get install python3
$ sudo apt-get install libboost-all-dev
$ sudo apt-get install libnuma-dev
$ sudo apt-get install zlib1g-dev
```

Next, you need to clone, build and install the libraries in the `External` subdirectory. To do so,
type the following commands at the top-level directory of the framework:

```
$ git submodule update --init
$ cd External/RoutingKit
$ make
$ sudo cp -r include lib /usr/local
```

## Building

Once you installed the packages, simply type `scons` at the top-level directory of the framework:

```
$ scons
```

## SCons Integration for Eclipse

The plugin SConsolidator provides tool integration for SCons in Eclipse.
Install it via its update site `http://www.sconsolidator.com/update`.

At the time of writing, there is a bug in the SConsolidator sources that you have to fix by hand.
If you are lucky, the bug will have been fixed in the official sources at the time you read this.
If not, try the following:

1. Find out where Eclipse installed the SConsolidator sources.
2. Open file `ch.hsr.ifs.sconsolidator.core/scons_files/BuildInfoCollector.py` with a text editor.
3. Change the definition of the function `get_compiler_flags` from

```
def get_compiler_flags(lang, environ):
    return (environ['CXXFLAGS'] if lang == 'c++' else environ['CCFLAGS'])
```

to

```
def get_compiler_flags(lang, environ):
    return [environ['CXXFLAGS'] if lang == 'c++' else environ['CCFLAGS']]
```

To import the framework as a SCons project into Eclipse, follow these steps:

1. Choose the menu `File`, `Import...`.
2. Select `C/C++`, `New SCons project from existing source`.
3. Choose the top-level directory of the framework.
4. Choose the menu `Project`, `Properties`, `C/C++ Build`, `Settings`, `Binary Parsers`.
5. Enable `Elf Parser`.

