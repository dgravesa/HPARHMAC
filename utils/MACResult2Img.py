# Daniel Graves
# MACResult2Img.py
#
# Script takes vectorized output of HMACTool and creates segmented image based on pixel modes

import sys
import getopt
import numpy
import cv2

def main(argv):
  usage = "usage: python MACResult2Img.py -i inputfile -d width,height -o outputimage\n\n"
  inputfile = ''
  dimension = [0, 0]
  dim_given = 0
  outputimage = ''

  # parse command line arguments
  try:
    opts, args = getopt.getopt(argv, "i:d:o:", ["input=", "dimension=", "output="])
  except getopt.GetoptError:
    sys.stderr.write(usage)
    sys.exit(1)
  for opt, arg in opts:
    if opt in ("-i", "--input"):
      inputfile = arg
    elif opt in ("-d", "--dimension"):
      dim_given = 1
      dimension[0] = int(arg.split(',')[0])
      dimension[1] = int(arg.split(',')[1])
    elif opt in ("-o", "--output"):
      outputimage = arg
    else:
      sys.stderr.write(usage)
      sys.exit(2)

  # check for arguments
  if (not inputfile) or (not outputimage) or (not dim_given):
    sys.stderr.write(usage)
    sys.exit(3)

  # open input file and get header
  fin = open(inputfile, "r")
  header = fin.readline().split()
  nmodes = int(header[0])
  ncoords = int(header[1])
  nfields = int(header[2])

  # check data dimension
  if (ncoords != dimension[0] * dimension[1]):
    sys.stderr.write("error: specified dimension does not match data size\n")
    sys.exit(4)

  # get modes from input file
  modes = []
  for i in range(0, nmodes):
    mode = [int(float(x)) for x in fin.readline().split()]
    modes.append(mode)

  # create image
  sys.stdout.write("creating image...\n")
  img = numpy.zeros((dimension[1], dimension[0], nfields), numpy.uint8)
  for i in range(img.shape[0]):
    for j in range(img.shape[1]):
      index = int(fin.readline().split()[0])
      img[i, j] = modes[index]

  # write image
  sys.stdout.write("writing '" + outputimage + "'...\n")
  cv2.imwrite(outputimage, img)

if __name__ == "__main__":
  main(sys.argv[1:])

