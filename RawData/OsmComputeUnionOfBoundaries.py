#!/usr/bin/python3

import argparse
import osmium
from shapely import geometry
from shapely import wkb

def parseCommandLine():
  ap = argparse.ArgumentParser(
    description='Extract the specified boundaries from an OSM file, compute their union, and '
    'output the largest component as a POLY file.')
  ap.add_argument('-o', '--output', default='-', type=argparse.FileType('w'),
                  help="The name of the output file. Default is '-' (STDOUT).", metavar='FILE')
  ap.add_argument('osm_file', help='The name of the OSM file.', metavar='osm-file')
  ap.add_argument('boundary_id', nargs='+', type=int,
                  help='The IDs of the boundaries that are to be joined.', metavar='boundary-id')
  return ap.parse_args()

class BoundaryHandler(osmium.SimpleHandler):
  
  def __init__(self, boundaryIds):
    osmium.SimpleHandler.__init__(self)
    self.boundaryIds = boundaryIds
    self.union = None
    self.wkbFactory = osmium.geom.WKBFactory()
  
  def area(self, osmMultipolygon):
    if osmMultipolygon.id // 2 in self.boundaryIds:
      self.boundaryIds.remove(osmMultipolygon.id // 2)
      bin = self.wkbFactory.create_multipolygon(osmMultipolygon)
      multipolygon = wkb.loads(bin, True)
      if not multipolygon.is_valid:
        print('Invalid boundary: ' + osmMultipolygon.id // 2)
      elif self.union is None:
        self.union = multipolygon
      else:
        self.union = self.union.union(multipolygon)
  
  def getLargestComponent(self):
    if self.union is None:
      return None
    elif self.union.geom_type == 'Polygon':
      return geometry.Polygon(self.union.exterior)
    elif self.union.geom_type == 'MultiPolygon':
      largestComp = None
      largestCompArea = 0
      for comp in self.union:
        compWithoutHoles = geometry.Polygon(comp.exterior)
        if abs(compWithoutHoles.area) > largestCompArea:
          largestComp = compWithoutHoles
          largestCompArea = abs(compWithoutHoles.area)
      return largestComp
    else:
      print('Unknown type: ' + self.union.geom_type)

def writePolyFile(polygon, file):
  print('(c) OpenStreetMap contributors', file=file)
  print('1', file=file)
  for vertex in polygon.exterior.coords:
    print('\t{lng:.7f}\t{lat:.7f}'.format(lng=vertex[0], lat=vertex[1]), file=file)
  print('END', file=file)
  print('END', file=file)

def main():
  args = parseCommandLine()
  handler = BoundaryHandler(set(args.boundary_id))
  handler.apply_file(args.osm_file, True)
  if handler.boundaryIds:
    print('Missing boundaries:', end='')
    for id in handler.boundaryIds:
      print(' {}'.format(id), end='')
    print()
  comp = handler.getLargestComponent()
  if comp is not None:
    writePolyFile(comp, args.output)

if __name__ == '__main__':
  main()
