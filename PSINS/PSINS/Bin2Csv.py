import csv
import struct
import math
DIR = r'D:\Data'
InputFile = open(DIR + r'\kf.bin', 'rb')
OutputFile = csv.writer(open(DIR + r'\kf.csv', 'w', newline=''))

while True:
    datas = []
    data = InputFile.read(8)
    if not data:
        break
    data = struct.unpack('d', data)[0]
    datas.append(data)
    for i in range(32):
        data = InputFile.read(8)
        data = struct.unpack('d', data)[0]
        # if i == 5 or i == 6:
        #     data = math.degrees(data)
        datas.append(data)
    OutputFile.writerow(datas)