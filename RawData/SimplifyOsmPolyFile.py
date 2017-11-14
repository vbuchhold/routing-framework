#!/usr/bin/env python3

import argparse

from shapely.geometry import MultiPolygon

def parseCommandLine():
  'Parses the command line.'
  ap = argparse.ArgumentParser(
      description='This program reads an OSM POLY file from disk, simplifies the polygons and '
      'outputs the result as an OSM POLY file.')
  ap.add_argument('-t', '--tolerance', default=0.0, type=float,
                  help='the maximum tolerance for the simplification algorithm', metavar='FLOAT')
  ap.add_argument('infile', type=argparse.FileType('rt'), help='the input POLY file')
  ap.add_argument('outfile', type=argparse.FileType('wt'), help='the output POLY file')
  args = ap.parse_args()
  return args

def parseOsmPolyFile(infile):
  'Reads the specified OSM POLY file from disk.'
  inPolygon = False
  coordinates = []
  fileAndSectionNames = []
  fileAndSectionNames.append(next(infile).strip())
  for line in infile:
    line = line.strip()
    if not inPolygon and line == 'END':
      # End-of-file occurred.
      area = MultiPolygon(coordinates)
      assert area.is_valid
      return area, iter(fileAndSectionNames)
    elif not inPolygon:
      # Section header occurred.
      inPolygon = True
      fileAndSectionNames.append(line)
      if line[0] == '!':
        coordinates[-1][1].append([])
        currentPolygon = coordinates[-1][1][-1]
      else:
        coordinates.append(([], []))
        currentPolygon = coordinates[-1][0]
    elif inPolygon and line == 'END':
      # End-of-section occurred.
      inPolygon = False
    elif inPolygon:
      # Vertex record occurred.
      tokens = line.split()
      assert len(tokens) == 2
      currentPolygon.append((float(tokens[0]), float(tokens[1])))
  assert False

def writeClosedChainToOsmPolyFile(outfile, closedChain, name):
  'Writes the specified closed polygonal chain to the given OSM POLY file.'
  print(name, file=outfile)
  if not closedChain is None:
    for vertex in closedChain.coords:
      print('  {lng:E}  {lat:E}'.format(lng=vertex[0], lat=vertex[1]), file=outfile)
  print('END', file=outfile)

def main():
  args = parseCommandLine()
  area, fileAndSectionNames = parseOsmPolyFile(args.infile)
  simplifiedPolygons = []
  for polygon in area:
    simplifiedPolygons.append(polygon.simplify(args.tolerance))
  area = MultiPolygon(simplifiedPolygons)
  assert area.is_valid
  numVertices = 0
  print(next(fileAndSectionNames), file=args.outfile)
  for polygon in area:
    writeClosedChainToOsmPolyFile(args.outfile, polygon.exterior, next(fileAndSectionNames))
    numVertices += len(polygon.exterior.coords)
    for hole in polygon.interiors:
      writeClosedChainToOsmPolyFile(args.outfile, hole, next(fileAndSectionNames))
      numVertices += len(hole.coords)
  print('END', file=args.outfile)
  print('Number of Vertices:', numVertices)

if __name__ == '__main__':
  main()
