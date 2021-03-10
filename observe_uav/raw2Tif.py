import numpy as np
import matplotlib.pyplot as plt
import cv2
import os
import struct
import PIL.Image as Image

def getFileName(file_dir,ref):
        L=[]
        for root, dirs, files in os.walk(file_dir):
            for file in files:
                if os.path.splitext(file)[1] == ref:
                    L.append(os.path.join(root, file))
        return L

def readimage(infile,ns,nl,nb):
    fb = open(infile,'rb')
    mydata = np.zeros((ns,nl,nb))
    for kb in range(nb):
        for ks in range(ns):
            for kl in range(nl):
                arr = fb.read(2)
                elem = struct.unpack('h',arr)[0]
                mydata[ks][kl][kb] = elem
    return mydata

fileDir = r"D:\Data\UAVLST\huailai083010\338\\"

# read and save as jpg
fileNames=getFileName(fileDir,'.raw')
fileNum = len(fileNames)
#fileNum = 3
im=[]
off = 0

# for i in range(fileNum):
for i in range(400,1000):

    infile = fileNames[i+off]
    print(infile)
    fsize = os.path.getsize(infile)
    fsize = fsize / float(1024)
    if (fsize < 25) :
        continue
    img0 = readimage(infile, 512, 640, 1)[:, :, 0]
    outfile = os.path.splitext(infile)[0] + '.tif'
    img0 = np.uint16(img0)
    cv2.imwrite(outfile,img0)

print(1)
#read as image
# fileNames=getFileName(fileDir,'.jpg')
# num = fileNum
# fileNum = len(fileNames)
# if num < fileNum:
#     fileNum = num
# for i in range(fileNum):
#     infile = fileNames[i]
#     print(infile)
#     img0 = cv2.imread(infile)
#     im.append(img0)
# print(2)
#
# stitcher = cv2.createStitcher(False)
# fg,result=stitcher.stitch(im)
# cv2.imwrite(fileDir+'\\result.jpg',result)
# print(np.shape(result))
# print(3)

# cv2.imshow("Image", result)
# cv2.waitKey (0)
# cv2.destroyAllWindows()


'''
infile1 = r'D:\BaiduNetdiskDownload\dali\DIR004\MID000\00112.raw'
img1 = readimage(infile1,512,640,1)[:,:,0]
plt.imshow(img1)
plt.axis('off')
plt.show()
'''