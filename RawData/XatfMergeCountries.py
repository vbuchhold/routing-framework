#!/usr/bin/env python3

import argparse
import csv
import os

def parseCommandLine():
  'Parses the command line.'
  # The countries contained in the standard DIMACS Europe instance.
  DEFAULT_COUNTRIES = [
    'aut', 'bel', 'che', 'deu', 'dnk', 'esp', 'fra',
    'gbr', 'ita', 'lux', 'nld', 'nor', 'prt', 'swe', 'ferries']

  # Parse the command line.
  ap = argparse.ArgumentParser(
    description='Merges the XATF files of the specified countries into a single XATF instance. '
    'When no countries are specified, outputs the standard DIMACS Europe instance.')
  ap.add_argument('-m', '--merge', action='store_true', help='merge coinciding vertices')
  ap.add_argument('-r', '--remove', action='store_true', help='remove self-loops')
  ap.add_argument('-I', '--input-directory', default='.',
                  help='XATF directory containing one subdirectory per country', metavar='DIR')
  ap.add_argument('-O', '--output-directory',
                  help='place the resulting instance into DIR', metavar='DIR')
  ap.add_argument('-o', '--output', help='the name of the resulting instance', metavar='NAME')
  ap.add_argument('country', nargs='*', default=DEFAULT_COUNTRIES,
                  help='three-letter ISO codes of the countries to be merged')
  args = ap.parse_args()

  # The --output-directory option defaults to the value of the --input-directory option.
  if args.output_directory is None:
    args.output_directory = args.input_directory
  return args

def main():
  args = parseCommandLine()

  # Associate field names with indices in the vertex/edge record.
  VERTEX_ID = 0
  LONGITUDE = 3
  LATITUDE = 4
  EDGE_ID = 0
  EDGE_TAIL = 1
  EDGE_HEAD = 2

  idsSeen = set()  # The vertex IDs seen so far.
  latLngsSeen = {} # The LatLngs seen so far, with a representative vertex ID for each.
  changedIds = {}  # The vertex IDs that need to be changed, with their respective new values.

  # Read the vertex records and discard unwanted ones.
  vertexRecords = []
  for country in sorted(args.country):
    filename = os.path.join(args.input_directory, country, country + '.nfx')
    if os.path.exists(filename):
      print('Reading', country + '.nfx...', end='', flush=True)
      with open(filename, 'r') as infile:
        for record in csv.reader(infile):
          # Eliminate duplicate vertex IDs.
          vertexId = record[VERTEX_ID]
          if vertexId in idsSeen:
            continue
          idsSeen.add(vertexId)

          # Merge coinciding vertices.
          latLng = int(record[LATITUDE]), int(record[LONGITUDE])
          if args.merge and latLng in latLngsSeen:
            changedIds[vertexId] = latLngsSeen[latLng]
            continue
          else:
            latLngsSeen[latLng] = vertexId
          vertexRecords.append(record)
      print(' done.')

  # Free up some variables that possibly occupy much storage.
  del latLngsSeen
  vertexRecords.sort(key=lambda record: int(record[VERTEX_ID]))

  if not args.output is None:
    # Create the directory containing the output files.
    out = os.path.join(args.output_directory, args.output)
    if not os.path.exists(out):
      os.mkdir(out)

    # Write the single sorted vertex file.
    print('Writing', args.output + '.nfx...', end='', flush=True)
    with open(os.path.join(out, args.output + '.nfx'), 'w') as outfile:
      csv.writer(outfile, lineterminator=os.linesep).writerows(vertexRecords)
    print(' done.')
  del vertexRecords

  # Read the edge records and change tail and head IDs that need to be changed.
  edgeRecords = []
  for country in sorted(args.country):
    print('Reading', country + '.sfx...', end='', flush=True)
    with open(os.path.join(args.input_directory, country, country + '.sfx'), 'r') as infile:
      for record in csv.reader(infile):
        if record[EDGE_TAIL] in idsSeen and record[EDGE_HEAD] in idsSeen:
          record[EDGE_TAIL] = changedIds.get(record[EDGE_TAIL], record[EDGE_TAIL])
          record[EDGE_HEAD] = changedIds.get(record[EDGE_HEAD], record[EDGE_HEAD])
          if not args.remove or record[EDGE_TAIL] != record[EDGE_HEAD]:
            edgeRecords.append(record)
    print(' done.')
  edgeRecords.sort(key=lambda record: int(record[EDGE_ID]))

  if not args.output is None:
    # Write the single sorted edge file.
    print('Writing', args.output + '.sfx...', end='', flush=True)
    with open(os.path.join(out, args.output + '.sfx'), 'w') as outfile:
      csv.writer(outfile, lineterminator=os.linesep).writerows(edgeRecords)
    print(' done.')

if __name__ == '__main__':
  main()
