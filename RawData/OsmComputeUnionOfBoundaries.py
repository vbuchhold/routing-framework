#!/usr/bin/python3

import argparse
import json
import osmium
from shapely import geometry
from shapely import wkb
from shapely.geometry import multipolygon

def parseCommandLine():
  ap = argparse.ArgumentParser(
    description='Extract the specified boundaries from an OSM file, compute their union, and '
    'output the largest component as a POLY file.')
  ap.add_argument('-c', '--coastline', type=argparse.FileType(), help="Compute the intersection of "
                  "the union boundary and the specified coastline data.", metavar='GEOJSON-FILE')
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

def readCoastlineDataFrom(file):
  coastlineData = []
  for feature in json.load(file)['features']:
    coastlineData.append(geometry.shape(feature['geometry']))
  return multipolygon.MultiPolygon(coastlineData)

def getLargestComponent(boundary):
  if boundary.geom_type == 'Polygon':
    return geometry.Polygon(boundary.exterior)
  elif boundary.geom_type == 'MultiPolygon':
    largestComp = None
    largestCompArea = 0
    for comp in boundary:
      compWithoutHoles = geometry.Polygon(comp.exterior)
      if abs(compWithoutHoles.area) > largestCompArea:
        largestComp = compWithoutHoles
        largestCompArea = abs(compWithoutHoles.area)
    return largestComp
  else:
    print('Unknown type: ' + boundary.geom_type)

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
  boundary = handler.union
  
  if boundary is not None and args.coastline is not None:
    coastline = readCoastlineDataFrom(args.coastline)
    if coastline is not None:
      boundary = boundary.intersection(coastline)
  
  if boundary is not None:
    writePolyFile(getLargestComponent(boundary), args.output)

if __name__ == '__main__':
  main()
