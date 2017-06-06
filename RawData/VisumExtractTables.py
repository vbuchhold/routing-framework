#!/usr/bin/env python3

import argparse
import os

def parseCommandLine():
  'Parses the command line.'
  ap = argparse.ArgumentParser(
    description='This program reads a Visum NET file and stores its tables as single CSV files.')
  ap.add_argument('-o', '--output-directory', default='.',
                  help='place the output CSV files in DIR', metavar='DIR')
  ap.add_argument('infile', type=argparse.FileType(encoding='latin-1'),
                  help='the Visum network file (file extension: .net)')
  args = ap.parse_args()
  return args

def main():
  args = parseCommandLine()
  args.infile.readline()
  outfile = None
  for line in args.infile:
    if line[0] == '*':
      # Comment. Skip it.
      pass
    elif line[0] == '$':
      # Begin of a new table. Open a new output file.
      name, columns = tuple(line.split(':'))
      print('Writing', name[1:] + '.csv...', end='', flush=True)
      outfile = open(os.path.join(args.output_directory, name[1:] + '.csv'), 'w')
      outfile.write(columns)
    elif not line.strip():
      # End of the current table. Close the current output file.
      outfile.close()
      print(' done.')
    else:
      # Normal record line. Write it to the current output file.
      outfile.write(line)

if __name__ == '__main__':
  main()
