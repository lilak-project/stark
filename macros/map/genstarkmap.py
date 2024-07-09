import math

#input order
#   0 det_type
#   1 det_no
#   2 cobo_no
#   3 asad_no
#   4 zap_junc_no
#   5 zap_ohmic_no
#   6 mark_no
#   7 mark_id
#   8 FB (F=1,B=2)
#   9 ring_type
#  10 dE-E(E=1,dE=2)
#  11 ring_id
#  12 z_distance (cm)

#output order
#  cobo
#  asad
#  aget
#  chanfpn
#  chan
#  det_type
#  det_no
#  junc_ohmic
#  strip
#  updown
#  angle_theta
#  angle_phi

#detector parameters for X6

det_width = 40.3 # mm
det_height = 75. # mm
det_thickness = 5 # mm

outputlist = []
outputline = []
outputlist2 = []
outputlist3 = []
anglelist = []

my_file = open("starkmap_input.txt", "r") 
det_file = open("starkmap_detector.txt", "w") 
  
# reading the file 
data = my_file.readlines()

count_id = 0
  
for line in data:
  outputlist = []
  outputlist2 = []
  anglelist = []
  columns = line.replace('\n','').split("\t")

  # printing the data 
  # print(columns) 

  channelbegin=0
  if int(columns[8]) == 1:
    thetabegin = 0
  elif int(columns[8]) == 2:
    thetabegin = 90

  det_distance = float(columns[12])
  if int(columns[9]) == 12:
    if int(columns[10]) == 1:
      #det_radius = 15.2
      det_radius = 101.5
      layer_no = 1
    if int(columns[10]) == 2:
      det_radius = 90.
      #det_radius = 14.8
      layer_no = 0
  elif int(columns[9]) == 16:
    if int(columns[10]) == 1:
      #det_radius = 202.
      det_radius = 150.
      layer_no = 2
      #layer_no = 3
    if int(columns[10]) == 2:
      det_radius = 198.
      layer_no = 2
  theta1 = math.degrees(math.atan(det_radius/det_distance))
  angle_theta = thetabegin + float("{:.2f}".format(theta1))

  det_type = columns[0]
  det_id = count_id
  count_id = count_id+1
  #det_id = columns[1]
  ring_id = int(columns[11]) 
  if int(columns[9])==12:
      if   ring_id==1:  ring_index = 0
      elif ring_id==2:  ring_index = 11
      elif ring_id==3:  ring_index = 10
      elif ring_id==4:  ring_index = 9
      elif ring_id==5:  ring_index = 8
      elif ring_id==6:  ring_index = 7
      elif ring_id==7:  ring_index = 6
      elif ring_id==8:  ring_index = 5
      elif ring_id==9:  ring_index = 4
      elif ring_id==10: ring_index = 3
      elif ring_id==11: ring_index = 2
      elif ring_id==12: ring_index = 1
  elif int(columns[9])==16:
      if   ring_id==1:  ring_index = 0
      elif ring_id==2:  ring_index = 15
      elif ring_id==3:  ring_index = 14
      elif ring_id==4:  ring_index = 13
      elif ring_id==5:  ring_index = 12
      elif ring_id==6:  ring_index = 11
      elif ring_id==7:  ring_index = 10
      elif ring_id==8:  ring_index = 9
      elif ring_id==9:  ring_index = 8
      elif ring_id==10: ring_index = 7
      elif ring_id==11: ring_index = 6
      elif ring_id==12: ring_index = 5
      elif ring_id==13: ring_index = 4
      elif ring_id==14: ring_index = 3
      elif ring_id==15: ring_index = 2
      elif ring_id==16: ring_index = 1
  dphi = math.degrees(math.atan((det_width/2.)/det_radius))
  phi0 = 90 + ring_index * (360/int(columns[9]))
  while phi0>=360: phi0 = phi0 - 360
  phi1 = 90 + ring_index * (360/int(columns[9])) - dphi
  phi2 = 90 + ring_index * (360/int(columns[9])) + dphi
  if columns[0] == "x6":
      num_sides = 2
      num_junction_strips = 8
      num_ohmic_strips = 4
      use_junction_lr = 1
      use_ohmic_lr = 0
  if columns[0] == "csd":
      num_sides = 2
      num_junction_strips = 8
      num_ohmic_strips = 1
      use_junction_lr = 0
      use_ohmic_lr = 0
  det_file.write(f"{det_type:<4} {det_id:>4} {num_sides:>4} {num_junction_strips:>4} {num_ohmic_strips:>4} {use_junction_lr:>4} {use_ohmic_lr:>4} {det_distance:>5} {layer_no:>4} {det_radius:>5} {ring_index:>4} {phi0:>8} {det_width:>5} {det_height:>5} {det_thickness:>4}\n")

  if columns[0] == "x6":
    total_channels = 16

    phibegin = (int(columns[11])-1) * (360/int(columns[9]))

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
    outputlist2.append(int(det_id)) # det_no
    outputlist2.append(1) # junction=1, ohmic=2
    outputlist2.append(0) # stripno
    #outputlist2.append(0) # updown, up=0, down=1
    outputlist2.append('u') # updown, up=0, down=1
    anglelist.append(angle_theta) # angle_theta
    anglelist.append(phibegin) # angle_phi
    for i in range(total_channels):
      outputlist3 = []
      if (channelbegin+i) < 11: 
        outputlist[-2] = channelbegin+i
      elif (channelbegin+i) < 21:
        outputlist[-2] = channelbegin+i+1
      elif (channelbegin+i) < 43:
        outputlist[-2] = channelbegin+i+2
      elif (channelbegin+i) < 54:
        outputlist[-2] = channelbegin+i+3
      else:
        outputlist[-2] = channelbegin+i+4
      outputlist[-1] = channelbegin+i

      outputlist2[-2] = int(i/2)
      outputlist2[-1] = i%2

      angle_phi = phibegin + float("{:.2f}".format(int(i/2)*(360/int(columns[9])/(total_channels/2))))
      anglelist[-1] = angle_phi
      outputlist3 = outputlist + outputlist2 + anglelist
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
    outputlist2.append(int(det_id)) # det_no
    outputlist2.append(2) # junction=1, ohmic=2
    outputlist2.append(0) # stripno
    for i in range(4):
      if (channelbegin+i) < 11: 
        outputlist[-2] = channelbegin+i
      elif (channelbegin+i) < 21:
        outputlist[-2] = channelbegin+i+1
      elif (channelbegin+i) < 43:
        outputlist[-2] = channelbegin+i+2
      elif (channelbegin+i) < 54:
        outputlist[-2] = channelbegin+i+3
      else:
        outputlist[-2] = channelbegin+i+4
      outputlist[-1] = channelbegin+i
      outputlist3 = []
      outputlist2[-1] = int(i)
      anglelist[-1] = phibegin
      outputlist3 = outputlist + outputlist2 + [0] + anglelist # add ohmic updown wt 0
      outputline.append(outputlist3)

  elif columns[0] == "csd":
    total_channels = 8
    phibegin = (int(columns[11])-1) * (360/int(columns[9]))

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
    outputlist2.append(int(det_id)) # det_no
    outputlist2.append(1) # junction=1, ohmic=2
    outputlist2.append(0) # stripno
    anglelist.append(angle_theta) # angle_theta
    anglelist.append(phibegin) # angle_phi
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

      angle_phi = phibegin + i*(360/int(columns[9])/total_channels)
      anglelist[-1] = angle_phi

      outputlist3 = outputlist + outputlist2 + [0] + anglelist # add juntion updown wt 0 and anglelist
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
    outputlist2.append(int(det_id)) # det_no
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
    anglelist[-1] = phibegin
    outputlist3 = []
    outputlist3 = outputlist + outputlist2 + anglelist
    outputline.append(outputlist3)
my_file.close() 

with open('starkmap_output.txt', 'w') as f:
    for line in outputline:
        f.write("%s\n" % str(line)[1:-1].replace("'","").replace(",","\t"))

print("Done")
