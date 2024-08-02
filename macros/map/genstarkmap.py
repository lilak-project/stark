import math

#input order
#det_type,det_no,cobo_no,asad_no,zap_junc_no,zap_ohmic_no,mark_no,mark_id,FB (F=1,B=2),ring_type,dE-E(E=1,dE=2),ring_id,ring_diameter (cm),z_distance (cm)

#output order
#cobo,asad,aget,chanfpn,chan,det_type,det_no,junc_ohmic,strip,updown,ring_radius,z_distance

#detector parameters for X6

outputlist = []
outputline = []
outputlist2 = []
outputlist3 = []
dimensionlist = []

my_file = open("starkmap_input.txt", "r") 
  
# reading the file 
data = my_file.readlines()
  
for line in data:
  outputlist = []
  outputlist2 = []
  dimensionlist = []
  columns = line.replace('\n','').split("\t")

  # printing the data 
  # print(columns) 

  channelbegin=0
  det_diameter = float(columns[12])
  det_distance = float(columns[13])
  det_radius = det_diameter/2.0

  #if int(columns[8]) == 1:
  #  thetabegin = 0
  #elif int(columns[8]) == 2:
  #  thetabegin = 90
  #angle_theta = thetabegin + float("{:.2f}".format(math.degrees(math.atan(det_radius/det_distance))))
  #phibegin = (int(columns[11])-1) * (360/int(columns[9]))
  #angle_phi = phibegin + float("{:.2f}".format(int(i/2)*(360/int(columns[9])/(total_channels/2))))
  #angle_phi = phibegin + i*(360/int(columns[9])/total_channels)

  if columns[0] == "x6":
    total_channels = 16


    outputlist.append(int(columns[2]))
    outputlist.append(int(columns[3]))
    if columns[4] == "a":
      outputlist.append(3)
    elif columns[4] == "b":
      outputlist.append(2)
    elif columns[4] == "c":
      outputlist.append(1)
    if int(columns[7]) == 1:
      channelbegin = 0
    elif int(columns[7]) == 2:
      channelbegin = 16 
    elif int(columns[7]) == 3:
      channelbegin = 32 
    elif int(columns[7]) == 4:
      channelbegin = 48 

    #print junction channels
    outputlist.append(channelbegin) #channel_wt_fpn
    outputlist.append(channelbegin) # channel
    outputlist2.append(columns[0]) # det_type
    outputlist2.append(int(columns[1])) # det_no
    outputlist2.append(1) # junction=1, ohmic=2
    outputlist2.append(0) # stripno
    outputlist2.append(0) # updown, up=0, down=1
    dimensionlist.append(det_radius) # det_radius
    dimensionlist.append(det_distance) # det_distance
    for i in range(total_channels):
      outputlist3 = []
      if (channelbegin+i) < 11: 
        outputlist[-2] = channelbegin+i
      elif (channelbegin+i) < 21:
        outputlist[-2] = channelbegin+i+1
      elif (channelbegin+i) < 43:
        outputlist[-2] = channelbegin+i+2
      elif (channelbegin+i) < 53:
        outputlist[-2] = channelbegin+i+3
      else:
        outputlist[-2] = channelbegin+i+4
      outputlist[-1] = channelbegin+i

      outputlist2[-2] = int(i/2)
      outputlist2[-1] = i%2

      outputlist3 = outputlist + outputlist2 + dimensionlist
      outputline.append(outputlist3)

    #print ohmic channels
    if columns[5] == "a":
      outputlist[-3] = 0
      channelbegin = 0
    elif columns[5] == "b":
      outputlist[-3] = 0
      channelbegin = 16 
    elif columns[5] == "c":
      outputlist[-3] = 0
      channelbegin = 32
    if int(columns[7]) == 1:
      channelbegin += 0
    elif int(columns[7]) == 2:
      channelbegin += 4
    elif int(columns[7]) == 3:
      channelbegin += 8
    elif int(columns[7]) == 4:
      channelbegin += 12
    outputlist2 = []
    outputlist2.append(columns[0]) # det_type
    outputlist2.append(int(columns[1])) # det_no
    outputlist2.append(2) # junction=1, ohmic=2
    outputlist2.append(0) # stripno
    for i in range(4):
      if (channelbegin+i) < 11: 
        outputlist[-2] = channelbegin+i
      elif (channelbegin+i) < 21:
        outputlist[-2] = channelbegin+i+1
      elif (channelbegin+i) < 43:
        outputlist[-2] = channelbegin+i+2
      elif (channelbegin+i) < 53:
        outputlist[-2] = channelbegin+i+3
      else:
        outputlist[-2] = channelbegin+i+4
      outputlist[-1] = channelbegin+i
      outputlist3 = []
      if i == 0:
        outputlist2[-1] = int(3)
      else:
        outputlist2[-1] = int(i-1)
      outputlist3 = outputlist + outputlist2 + [0] + dimensionlist # add ohmic updown wt 0
      outputline.append(outputlist3)

  elif columns[0] == "csd":
    total_channels = 8

    outputlist.append(int(columns[2]))
    outputlist.append(int(columns[3]))
    if columns[4] == "a":
      outputlist.append(3)
    elif columns[4] == "b":
      outputlist.append(2)
    elif columns[4] == "c":
      outputlist.append(1)
    if columns[6] == "mix1":
      if int(columns[7]) == 1:
        channelbegin = 0
      elif int(columns[7]) == 2:
        channelbegin = 8 
      elif int(columns[7]) == 3:
        channelbegin = 16
      elif int(columns[7]) == 4:
        channelbegin = 24 
    elif columns[6] == "mix2":
      if int(columns[7]) == 1:
        channelbegin = 32
      elif int(columns[7]) == 2:
        channelbegin = 40 
      elif int(columns[7]) == 3:
        channelbegin = 48 
      elif int(columns[7]) == 4:
        channelbegin = 56 

    #print junction channels
    outputlist.append(channelbegin) #channel_wt_fpn
    outputlist.append(channelbegin) # channel
    outputlist2.append(columns[0]) # det_type
    outputlist2.append(int(columns[1])) # det_no
    outputlist2.append(1) # junction=1, ohmic=2
    outputlist2.append(0) # stripno
    dimensionlist.append(det_radius) # det_radius
    dimensionlist.append(det_distance) # det_distance
    for i in range(total_channels):
      outputlist3 = []
      if (channelbegin+i) < 11: 
        outputlist[-2] = channelbegin+i
      elif (channelbegin+i) < 21:
        outputlist[-2] = channelbegin+i+1
      elif (channelbegin+i) < 43:
        outputlist[-2] = channelbegin+i+2
      elif (channelbegin+i) < 53:
        outputlist[-2] = channelbegin+i+3
      else:
        outputlist[-2] = channelbegin+i+4
      outputlist[-1] = channelbegin+i

      outputlist2[-1] = i

      outputlist3 = outputlist + outputlist2 + [0] + dimensionlist # add juntion updown wt 0 and dimensionlist
      outputline.append(outputlist3)

    #print ohmic channels
    if columns[5] == "a":
      outputlist[-3] = 0
      channelbegin = 0
    elif columns[5] == "b":
      outputlist[-3] = 0
      channelbegin = 16 
    elif columns[5] == "c":
      outputlist[-3] = 0
      channelbegin = 32
    if columns[6] == "mix1":
      if int(columns[7]) == 1:
        channelbegin += 2
      elif int(columns[7]) == 2:
        channelbegin += 0
      elif int(columns[7]) == 3:
        channelbegin += 6
      elif int(columns[7]) == 4:
        channelbegin += 4
    elif columns[6] == "mix2":
      if int(columns[7]) == 1:
        channelbegin += 10 
      elif int(columns[7]) == 2:
        channelbegin += 8
      elif int(columns[7]) == 3:
        channelbegin += 14
      elif int(columns[7]) == 4:
        channelbegin += 12
    outputlist2 = []
    outputlist2.append(columns[0]) # det_type
    outputlist2.append(int(columns[1])) # det_no
    outputlist2.append(2) # junction=1, ohmic=2
    outputlist2.append(0) # stripno
    outputlist2.append(0) # ohmic updown wt 0
    if (channelbegin) < 11: 
      outputlist[-2] = channelbegin
    elif (channelbegin) < 21:
      outputlist[-2] = channelbegin+1
    elif (channelbegin) < 43:
      outputlist[-2] = channelbegin+2
    elif (channelbegin) < 53:
      outputlist[-2] = channelbegin+3
    else:
      outputlist[-2] = channelbegin+4
    outputlist[-1] = channelbegin
    outputlist3 = []
    outputlist3 = outputlist + outputlist2 + dimensionlist
    outputline.append(outputlist3)
my_file.close() 

with open('starkmap_output.txt', 'w') as f:
    for line in outputline:
        f.write("%s\n" % str(line)[1:-1].replace("'","").replace(",","\t"))

print("Done")
