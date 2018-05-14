# Daniel Graves
# MACResultDiff.py
#
# Script takes two vectorized outputs of HMACTool and prints classification-wise difference

import sys
import getopt

def main(argv):
  #usage = "usage: MACResultDiff.py -i <input1,input2>\n\n"
  usage = "usage: MACResultDiff.py -1 input1 -2 input2\n\n"
  input1 = ''
  input2 = ''

  # parse command line arguments
  try:
    opts, args = getopt.getopt(argv, "1:2:", ["input1=", "input2="])
  except getopt.GetoptError:
    sys.stderr.write(usage)
    sys.exit(1)
  for opt, arg in opts:
    if opt in ("-1", "--input1"):
      input1 = arg
    elif opt in ("-2", "--input2"):
      input2 = arg
    else:
      sys.stderr.write(usage)
      sys.exit(2)

  # check for arguments
  if (not input1) or (not input2):
    sys.stderr.write(usage)
    sys.exit(3)

  # open input files
  fin1 = open(input1, "r")
  fin2 = open(input2, "r")

  # check headers
  header1 = [int(i) for i in fin1.readline().split()]
  header2 = [int(i) for i in fin2.readline().split()]
  if header1[1] != header2[1]:
    sys.stderr.write("error: data sizes do not match\n")
    sys.exit(4)
  elif header1[2] != header2[2]:
    sys.stderr.write("error: data dimensions do not match\n")
    sys.exit(5)
  elif header1[0] != header2[0]:
    sys.stderr.write("error: cluster counts must match in order to run difference\n")
    sys.exit(6)

  # ignore mode lines (assume roughly equal although slight differences are okay)
  for i in range(header1[0]):
    fin1.readline()
    fin2.readline()

  # count number of different classifications
  count = 0
  for i in range(header1[1]):
    value1 = int(fin1.readline().split()[0])
    value2 = int(fin2.readline().split()[0])
    if value1 != value2:
      count = count + 1

  # print difference
  pdiff = 100 * count / header1[1]
  print("difference =", pdiff, "%")

if __name__ == "__main__":
  main(sys.argv[1:])

