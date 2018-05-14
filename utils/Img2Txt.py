# Daniel Graves
# Img2Txt.py
#
# Script takes an image input and prints the values of each pixel

import sys
import getopt
import cv2

def main(argv):
  usage = "usage: python Img2Txt.py -i inputfile\n\n"
  inputfile = ''

  # parse command line arguments
  try:
    opts, args = getopt.getopt(argv, "i:", ["input="])
  except getopt.GetoptError:
    sys.stderr.write(usage)
    sys.exit(1)
  for opt, arg in opts:
    if opt in ("-i", "--input"):
      inputfile = arg
    else:
      sys.stderr.write(usage)
      sys.exit(2)

  # check for input file
  if not inputfile:
    sys.stderr.write(usage)
    sys.exit(3)

  # open image
  img = cv2.imread(inputfile)
  if img is None:
    error_string = "error: could not open image: " + inputfile + "\n"
    sys.stderr.write(error_string)
    sys.exit(4)

  # write pixel values
  for i in range(img.shape[0]):
    for j in range(img.shape[1]):
      for k in range(img.shape[2]):
        sys.stdout.write(str(img.item(i, j, k)) + " ")
      sys.stdout.write("\n")

if __name__ == "__main__":
  main(sys.argv[1:])

