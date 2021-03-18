# Routing Framework

This repository contains the C++17 source code used in the following publications/submissions:

* Valentin Buchhold, Peter Sanders, and Dorothea Wagner. Real-Time Traffic Assignment Using Fast
  Queries in Customizable Contraction Hierarchies. In Gianlorenzo D'Angelo, editor, *Proceedings of
  the 17th International Symposium on Experimental Algorithms (SEA'18)*, volume 103 of *Leibniz
  International Proceedings in Informatics (LIPIcs)*, pages 27:1–27:15. Schloss Dagstuhl, 2018.
  [doi:10.4230/LIPIcs.SEA.2018.27](http://dx.doi.org/10.4230/LIPIcs.SEA.2018.27).
* Valentin Buchhold, Peter Sanders, and Dorothea Wagner. Efficient Calculation of Microscopic
  Travel Demand Data with Low Calibration Effort. In Farnoush Banaei-Kashani, Goce Trajcevski, Ralf
  Hartmut Güting, Lars Kulik, and Shawn D. Newsam, editors, *Proceedings of the 27th ACM SIGSPATIAL
  International Conference on Advances in Geographic Information Systems (SIGSPATIAL'19)*, pages
  379-388. ACM Press, 2019.
  [doi:10.1145/3347146.3359361](http://dx.doi.org/10.1145/3347146.3359361).
* Valentin Buchhold, Peter Sanders, and Dorothea Wagner. Real-time Traffic Assignment Using
  Engineered Customizable Contraction Hierarchies. ACM Journal of Experimental Algorithmics,
  24(2):2.4:1-2.4:28, 2019. [doi:10.1145/3362693](http://dx.doi.org/10.1145/3362693).
* Valentin Buchhold and Dorothea Wagner. Nearest-Neighbor Queries in Customizable Contraction
  Hierarchies and Applications. In David Coudert and Emanuele Natale, editors, *Proceedings of the
  19th International Symposium on Experimental Algorithms (SEA'21)*, pages 18:1-18:18. Schloss
  Dagstuhl, 2021. [doi:10.4230/LIPIcs.SEA.2021.18](http://dx.doi.org/10.4230/LIPIcs.SEA.2021.18).

## License

All files in this repository except the files in the directory `External` are licensed under the MIT
license. External libraries are licensed under their respective licenses. Note that the compiled
programs `CalculateDemand` and `ComputeUnionBoundary` use libraries that are released under the GNU
GPLv3, and thus the compiled programs `CalculateDemand` and `ComputeUnionBoundary` have to be under
the GNU GPLv3.

## Prerequisites

You need to have some tools and libraries installed. On Debian and its derivatives (such as Ubuntu)
the `apt-get` tool can be used:

```
$ sudo apt-get install build-essential
$ sudo apt-get install cmake
$ sudo apt-get install python3
$ sudo apt-get install libboost-all-dev
$ sudo apt-get install libcairo2-dev
$ sudo apt-get install libcgal-dev
$ sudo apt-get install libproj-dev
$ sudo apt-get install zlib1g-dev
```

Next, you need to clone, build and install the libraries in the `External` subdirectory. To do so,
type the following commands at the top-level directory of the framework:

```
$ git submodule update --init
$ cd External
$ cd fast-cpp-csv-parser && sudo cp *.h /usr/local/include && cd ..
$ cd randomc && sudo mkdir /usr/local/include/randomc && sudo cp *.h $_ && cd ..
$ cd rapidxml && sudo cp *.hpp /usr/local/include && cd ..
$ cd RoutingKit && make && sudo cp -r include lib /usr/local && cd ..
$ cd stocc && sudo mkdir /usr/local/include/stocc && sudo cp *.h $_ && cd ..
$ cd vectorclass && sudo mkdir /usr/local/include/vectorclass && sudo cp *.h $_ && cd ..
```

## Building

Once you installed the packages, simply type the following commands at the top-level directory of
the framework:

```
$ cmake -S . -B Build/Devel
$ cmake --build Build/Devel
```

## Experiments in *Real-Time Traffic Assignment Using Fast Queries in Customizable Contraction Hierarchies*

To get the version of the source code used in the publication, check out the `SEA18` tag:

```
$ git checkout SEA18
```

To run the experiments presented in the publication, enter the following commands at the top-level
directory of the framework:

```
$ cd Publications/FastTA/SEA18
$ ./PrepareP2P <path-to-xatf-data> ~
$ ./PrepareTA <path-to-visum-data> <path-to-mobitopp-data> ~
$ ./ConductP2P ~
$ ./ConductTA ~
```

*Note: This repository contains the full source code used in the publication. However,
unfortunately the data used in the experiments is proprietary and not publicly available.*

## Experiments in *Real-Time Traffic Assignment Using Engineered Customizable Contraction Hierarchies*

To get the version of the source code used in the publication, check out the `JEA18` tag:

```
$ git checkout JEA18
```

To run the experiments presented in the publication, enter the following commands at the top-level
directory of the framework:

```
$ cd Publications/FastTA/JEA
$ ./PrepareP2P <path-to-stuttgart-data> <path-to-london-data> <path-to-xatf-data>     ~
$ ./PrepareTA  <path-to-stuttgart-data> <path-to-london-data> <path-to-mobitopp-data> ~
$ ./ConductP2P ~
$ ./ConductTA  ~
```

*Note: This repository contains the full source code used in the publication. However,
unfortunately the data used in the experiments is proprietary and not publicly available.*

## Experiments in *Efficient Calculation of Microscopic Travel Demand Data with Low Calibration Effort*

To get the version of the source code used in the publication, check out the `SIGSPATIAL19` tag:

```
$ git checkout SIGSPATIAL19
```

To run the experiments presented in the publication, create the following directory structure in a
freely chosen base directory `<base-dir>` and fill the directories with the respective raw data:

```
<base-dir>/Inputs/RawData/GermanCensus2011/
<base-dir>/Inputs/RawData/XATF/
<base-dir>/Inputs/RawData/GEOSTAT2011/
<base-dir>/Inputs/RawData/mobiTopp/
<base-dir>/Inputs/RawData/Visum/London/
<base-dir>/Inputs/RawData/Visum/Stuttgart/
```

*Note: This repository contains the full source code used in the publication. However,
unfortunately most data used in the experiments is proprietary and not publicly available.
Exceptions are the [German population grid](https://www.zensus2011.de/DE/Home/Aktuelles/DemografischeGrunddaten.html)
and the [European population grid](https://ec.europa.eu/eurostat/web/gisco/geodata/reference-data/population-distribution-demography/geostat).*

Now, enter the following command at the top-level directory of the framework:

```
$ Publications/DemandCalculation/SIGSPATIAL19/RunExperiments <base-dir>
```

## Experiments in *Nearest-Neighbor Queries in Customizable Contraction Hierarchies and Applications*

To get the version of the source code used in the publication, check out the `SEA21` tag:

```
$ git checkout SEA21
```

To run the experiments presented in the submission, enter the following command at the top-level
directory of the framework, with a freely chosen output directory `<base-directory>`:

```
$ Publications/NearestNeighbors/SEA21/RunExperiments \
>     <path-to-xatf-data> <path-to-euro-grid> <base-directory>
```
