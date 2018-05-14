# Daniel Graves
# ImgDiff.py
#
# Script takes 2 images as input and calculates the pixel-wise percentage difference

import sys
import getopt
import cv2

def main(argv):
  #usage = "usage: python ImgDiff.py -i <input1,input2>\n\n"
  usage = "usage: python ImgDiff.py -1 input1 -2 input2\n\n"
  inputfile1 = ''
  inputfile2 = ''

  # parse command line arguments
  try:
    #opts, args = getopt.getopt(argv, "i:", ["input="])
    opts, args = getopt.getopt(argv, "1:2:", ["input1=", "input2="])
  except getopt.GetoptError:
    sys.stderr.write(usage)
    sys.exit(1)
  for opt, arg in opts:
    #if opt in ("-i", "--input"):
    #  inputfile1 = arg.split(',')[0]
    #  inputfile2 = arg.split(',')[1]
    if opt in ("-1", "--input1"):
      inputfile1 = arg
    elif opt in ("-2", "--input2"):
      inputfile2 = arg
    else:
      sys.stderr.write(usage)
      sys.exit(2)

  # check for input file
  if (not inputfile1) or (not inputfile2):
    sys.stderr.write(usage)
    sys.exit(3)

  # open images
  img1 = cv2.imread(inputfile1)
  img2 = cv2.imread(inputfile2)
  if img1 is None or img2 is None:
    error_input = inputfile1
    if img1 is not None:
      error_input = inputfile2
    error_string = "error: could not open image: " + error_input + "\n"
    sys.stderr.write(error_string)
    sys.exit(4)

  # check image sizes
  if (img1.shape != img2.shape):
    sys.stderr.write("error: image shapes do not match\n")
    sys.exit(5)

  # count differing pixel values
  count = 0
  for i in range(img1.shape[0]):
    for j in range(img1.shape[1]):
      for k in range(img1.shape[2]):
        if img1[i, j, k] != img2[i, j, k]:
          count = count + 1
          break

  pdiff = 100 * count / (img1.shape[0] * img1.shape[1])
  print("image difference =", pdiff, "%")

if __name__ == "__main__":
  main(sys.argv[1:])

