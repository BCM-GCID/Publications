"""
Created on Tue Aug 24 13:12:06 2021

@author: piedra
"""
import numpy
import random
from random import randint

events = 10000
Dscan = 1        # = random walk speed of non-transcribing pol
k_transc = 5 # = linear speed of transcribing pol
Dbias = 1

genome_size = 11152                                      # genome length in nucleotides
pol_footprint = 28                                    # polymerase footprint in nucleotides
genome_chunks = round(genome_size/pol_footprint)    # genome length in polymerase footprints

GS = numpy.array([51, 1386, 2209, 3049, 4723])
GS_list = (['N', 'P', 'M', 'G', 'L'])
transc_prob = numpy.array([0.1, 0.5, 0.5, 0.5, 0.5])

GE = numpy.array([1383, 2196, 3036, 4716, 11102])
GE_list = (['N', 'P', 'M', 'G', 'L'])
transc_term_prob = numpy.array([0.99, 0.99, 0.99, 0.99, 0.99])

GS_chunks = numpy.round(GS/pol_footprint,0)
GE_chunks = numpy.round(GE/pol_footprint,0)

pol_ej = 0
pol_pos = numpy.array([0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.00])
pol_state = numpy.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

checks = numpy.size(pol_pos)-1

transc_events = numpy.array([0, 0, 0, 0, 0])
transc_term_events = numpy.array([0, 0, 0, 0, 0])
transc_RT_events = numpy.array([0, 0, 0, 0, 0])

for i in range(0, events):

    print()
    print('event #', i)
    print()

    for n in range(0, checks):  #arranges values of pol_pos array from lowest to highest and keeps pol_state values with associated pol_pos values
        for x in range(1, checks+1):
            if pol_pos[x] == 0 and pol_pos[x-1] != 0:
                pol_pos[x] = pol_pos[x-1]
                pol_pos[x-1] = 0
                pol_state[x] = pol_state[x-1]
                pol_state[x-1] = 0

    for x in range(0, numpy.size(pol_pos)):                                     #pol rebinding: checks pol_pos array from right to left for a 0. If found, changes the first 0 to 1
        if (False == numpy.logical_and(pol_pos > 0, pol_pos < 1).any()) and 1 not in pol_pos and pol_pos[numpy.size(pol_pos)-1-x] == 0:
            pol_pos[numpy.size(pol_pos)-1-x] = 1

    for i in range(0, numpy.size(pol_pos)):                                     #pol movement (scanning or transcription)
       GS_tobe = numpy.logical_and(GS_chunks > pol_pos[i], GS_chunks < pol_pos[i] + 1*Dscan*Dbias)
       GE_tobe = numpy.logical_and(GE_chunks > pol_pos[i], GE_chunks < pol_pos[i] + k_transc)
       if pol_pos[i] >= genome_chunks: #big block #1  #######################################
           pol_state[i] = 0
           pol_pos[i] = 0
           for n in range(0, checks):  #arranges values of pol_pos array from lowest to highest and keeps pol_state values with associated pol_pos values
               for x in range(0, checks+1):
                   if pol_pos[x] == 0 and pol_pos[x-1] != 0:
                       pol_pos[x] = pol_pos[x-1]
                       pol_pos[x-1] = 0
                       pol_state[x] = pol_state[x-1]
                       pol_state[x-1] = 0
       elif (pol_pos[i] >= 1) and (pol_pos[i] < genome_chunks) and (pol_state[i] == 0) and (pol_pos[i] not in GS_chunks) and (0 == numpy.logical_and(GS_chunks > pol_pos[i], GS_chunks < pol_pos[i] + 1*Dscan*Dbias).any()): #big block #2  #######################################
           a = randint(0,1)
           if (a == 0) and (pol_pos[i] - 1*Dscan not in pol_pos) and (pol_pos[i] - 1*Dscan >= 2):
               pol_pos[i] = pol_pos[i] - 1*Dscan
           elif (a == 1) and (pol_pos[i] + 1*Dscan*Dbias not in pol_pos) and (0 == numpy.logical_and(pol_pos > pol_pos[i], pol_pos < pol_pos[i] + 1*Dscan*Dbias).any()):
               pol_pos[i] = pol_pos[i] + 1*Dscan*Dbias
       elif (pol_pos[i] >= 1) and (pol_pos[i] < genome_chunks) and (pol_state[i] == 0) and (pol_pos[i] not in GS_chunks) and (1 == numpy.logical_and(GS_chunks > pol_pos[i], GS_chunks < pol_pos[i] + 1*Dscan*Dbias).any()) and (0 == numpy.logical_and(pol_pos > pol_pos[i], pol_pos < GS_chunks[list(GS_tobe).index(1)]).any()): #big block #3  #######################################
           b = randint(0,1)
           if (b == 0) and (pol_pos[i] - 1*Dscan not in pol_pos) and (pol_pos[i] - 1*Dscan >= 2):
               pol_pos[i] = pol_pos[i] - 1*Dscan
           elif (b == 1):
               d = random.choices(['yes','no'], [transc_prob[list(GS_tobe).index(1)], 1-transc_prob[list(GS_tobe).index(1)]])
               if (d == ['yes']):
                   pol_state[i] = 1
                   print('transcription_initiated @', GS_list[list(GS_tobe).index(1)])
                   transc_events[list(GS_tobe).index(1)] = transc_events[list(GS_tobe).index(1)] + 1
                   pol_pos[i] = GS_chunks[list(GS_tobe).index(1)] + k_transc
                   for n in range(0, checks):
                       for x in range(0, checks):
                           if pol_pos[numpy.size(pol_pos)-1-x] <= pol_pos[numpy.size(pol_pos)-2-x] and pol_pos[numpy.size(pol_pos)-1-x] != 0 and pol_state[numpy.size(pol_pos)-2-x] == 1 and pol_state[numpy.size(pol_pos)-1-x] == 0:
                               print('pol ejected from pos', pol_pos[numpy.size(pol_pos)-1-x])
                               pol_ej = pol_ej +1
                               pol_pos[numpy.size(pol_pos)-1-x] = pol_pos[numpy.size(pol_pos)-2-x]
                               pol_pos[numpy.size(pol_pos)-2-x] = 0
                               pol_state[numpy.size(pol_pos)-2-x] = 0
                               pol_state[numpy.size(pol_pos)-1-x] = 1
                           elif pol_pos[numpy.size(pol_pos)-1-x] < pol_pos[numpy.size(pol_pos)-2-x] and pol_pos[numpy.size(pol_pos)-1-x] != 0 and pol_state[numpy.size(pol_pos)-1-x] == 0 and pol_state[numpy.size(pol_pos)-2-x] == 0:
                               pol_pos[numpy.size(pol_pos)-2-x:numpy.size(pol_pos)-x].sort()
                               pol_state[numpy.size(pol_pos)-2-x] = 0
                               pol_state[numpy.size(pol_pos)-1-x] = 0
                           elif pol_pos[numpy.size(pol_pos)-1-x] <= pol_pos[numpy.size(pol_pos)-2-x] and pol_pos[numpy.size(pol_pos)-1-x] != 0 and pol_state[numpy.size(pol_pos)-2-x] == 0:
                               pol_state[numpy.size(pol_pos)-2-x] = 0
                               pol_state[numpy.size(pol_pos)-1-x] = 0
               elif (d == ['no']) and (0 == numpy.logical_and(pol_pos > pol_pos[i], pol_pos < pol_pos[i] + 1*Dscan*Dbias).any()):
                   pol_pos[i] = pol_pos[i] + 1*Dscan*Dbias
       elif (pol_pos[i] >= 1) and (pol_pos[i] < genome_chunks) and (pol_state[i] == 0) and (pol_pos[i] in GS_chunks): #big block #4  #######################################
           b = random.choices(['yes','no'], [transc_prob[list(GS_chunks).index(numpy.round(pol_pos[i]))], 1-transc_prob[list(GS_chunks).index(numpy.round(pol_pos[i]))]])
           if (b == ['yes']):
                 pol_state[i] = 1
                 print('transcription_initiated @', GS_list[list(GS_chunks).index(numpy.round(pol_pos[i]))])
                 transc_events[list(GS_chunks).index(numpy.round(pol_pos[i]))] = transc_events[list(GS_chunks).index(numpy.round(pol_pos[i]))] + 1
                 pol_pos[i] = pol_pos[i] + k_transc
                 for n in range(0, checks):
                     for x in range(0, checks):
                         if pol_pos[numpy.size(pol_pos)-1-x] <= pol_pos[numpy.size(pol_pos)-2-x] and pol_pos[numpy.size(pol_pos)-1-x] != 0 and pol_state[numpy.size(pol_pos)-2-x] == 1 and pol_state[numpy.size(pol_pos)-1-x] == 0:
                             print('pol ejected from pos', pol_pos[numpy.size(pol_pos)-1-x])
                             pol_ej = pol_ej +1
                             pol_pos[numpy.size(pol_pos)-1-x] = pol_pos[numpy.size(pol_pos)-2-x]
                             pol_pos[numpy.size(pol_pos)-2-x] = 0
                             pol_state[numpy.size(pol_pos)-2-x] = 0
                             pol_state[numpy.size(pol_pos)-1-x] = 1
                         elif pol_pos[numpy.size(pol_pos)-1-x] < pol_pos[numpy.size(pol_pos)-2-x] and pol_pos[numpy.size(pol_pos)-1-x] != 0 and pol_state[numpy.size(pol_pos)-1-x] == 0 and pol_state[numpy.size(pol_pos)-2-x] == 0:
                             pol_pos[numpy.size(pol_pos)-2-x:numpy.size(pol_pos)-x].sort()
                             pol_state[numpy.size(pol_pos)-2-x] = 0
                             pol_state[numpy.size(pol_pos)-1-x] = 0
                         elif pol_pos[numpy.size(pol_pos)-1-x] <= pol_pos[numpy.size(pol_pos)-2-x] and pol_pos[numpy.size(pol_pos)-1-x] != 0 and pol_state[numpy.size(pol_pos)-2-x] == 0:
                             pol_state[numpy.size(pol_pos)-2-x] = 0
                             pol_state[numpy.size(pol_pos)-1-x] = 0
           elif (b == ['no']):
                 c = randint(0,1)
                 if (c == 0) and (pol_pos[i] - 1*Dscan not in pol_pos) and (pol_pos[i] - 1*Dscan >= 2):
                     pol_pos[i] = pol_pos[i] - 1*Dscan
                 elif (c ==1) and (pol_pos[i] + 1*Dscan*Dbias not in pol_pos) and (0 == numpy.logical_and(pol_pos > pol_pos[i], pol_pos < pol_pos[i] + 1*Dscan*Dbias).any()):
                     pol_pos[i] = pol_pos[i] + 1*Dscan*Dbias
       elif (pol_pos[i] >= 1) and (pol_pos[i] < genome_chunks) and (pol_state[i] == 1) and (pol_pos[i] in GE_chunks): #big block #5  #######################################
           b = random.choices(['yes','no'], [transc_term_prob[list(GE_chunks).index(numpy.round(pol_pos[i]))], 1-transc_term_prob[list(GE_chunks).index(numpy.round(pol_pos[i]))]])
           if b == ['yes']:
              pol_state[i] = 0
              print('transcription_TERMINATED @', GE_list[list(GE_chunks).index(pol_pos[i])])
              transc_term_events[list(GE_chunks).index(pol_pos[i])] = transc_term_events[list(GE_chunks).index(pol_pos[i])] + 1
           elif b == ['no']:
              pol_state[i] = 1
              print('transc RT @', GE_list[list(GE_chunks).index(pol_pos[i])])
              transc_RT_events[list(GE_chunks).index(pol_pos[i])] = transc_RT_events[list(GE_chunks).index(pol_pos[i])] + 1
              pol_pos[i] = pol_pos[i] + k_transc
              for n in range(0, checks):
                     for x in range(0, checks):
                         if pol_pos[numpy.size(pol_pos)-1-x] <= pol_pos[numpy.size(pol_pos)-2-x] and pol_pos[numpy.size(pol_pos)-1-x] != 0 and pol_state[numpy.size(pol_pos)-2-x] == 1 and pol_state[numpy.size(pol_pos)-1-x] == 0:
                             print('pol ejected from pos', pol_pos[numpy.size(pol_pos)-1-x])
                             pol_ej = pol_ej +1
                             pol_pos[numpy.size(pol_pos)-1-x] = pol_pos[numpy.size(pol_pos)-2-x]
                             pol_pos[numpy.size(pol_pos)-2-x] = 0
                             pol_state[numpy.size(pol_pos)-2-x] = 0
                             pol_state[numpy.size(pol_pos)-1-x] = 1
                         elif pol_pos[numpy.size(pol_pos)-1-x] < pol_pos[numpy.size(pol_pos)-2-x] and pol_pos[numpy.size(pol_pos)-1-x] != 0 and pol_state[numpy.size(pol_pos)-1-x] == 0 and pol_state[numpy.size(pol_pos)-2-x] == 0:
                             pol_pos[numpy.size(pol_pos)-2-x:numpy.size(pol_pos)-x].sort()
                             pol_state[numpy.size(pol_pos)-2-x] = 0
                             pol_state[numpy.size(pol_pos)-1-x] = 0
                         elif pol_pos[numpy.size(pol_pos)-1-x] <= pol_pos[numpy.size(pol_pos)-2-x] and pol_pos[numpy.size(pol_pos)-1-x] != 0 and pol_state[numpy.size(pol_pos)-2-x] == 0:
                             pol_state[numpy.size(pol_pos)-2-x] = 0
                             pol_state[numpy.size(pol_pos)-1-x] = 0
       elif (pol_pos[i] >= 1) and (pol_pos[i] < genome_chunks) and (pol_state[i] == 1) and (0 == numpy.logical_and(GE_chunks > pol_pos[i], GE_chunks < pol_pos[i] + k_transc).any()): #big block #6  #######################################
           pol_pos[i] = pol_pos[i] + k_transc
           for n in range(0, checks):
                     for x in range(0, checks):
                         if pol_pos[numpy.size(pol_pos)-1-x] <= pol_pos[numpy.size(pol_pos)-2-x] and pol_pos[numpy.size(pol_pos)-1-x] != 0 and pol_state[numpy.size(pol_pos)-2-x] == 1 and pol_state[numpy.size(pol_pos)-1-x] == 0:
                             print('pol ejected from pos', pol_pos[numpy.size(pol_pos)-1-x])
                             pol_ej = pol_ej +1
                             pol_pos[numpy.size(pol_pos)-1-x] = pol_pos[numpy.size(pol_pos)-2-x]
                             pol_pos[numpy.size(pol_pos)-2-x] = 0
                             pol_state[numpy.size(pol_pos)-2-x] = 0
                             pol_state[numpy.size(pol_pos)-1-x] = 1
                         elif pol_pos[numpy.size(pol_pos)-1-x] < pol_pos[numpy.size(pol_pos)-2-x] and pol_pos[numpy.size(pol_pos)-1-x] != 0 and pol_state[numpy.size(pol_pos)-1-x] == 0 and pol_state[numpy.size(pol_pos)-2-x] == 0:
                             pol_pos[numpy.size(pol_pos)-2-x:numpy.size(pol_pos)-x].sort()
                             pol_state[numpy.size(pol_pos)-2-x] = 0
                             pol_state[numpy.size(pol_pos)-1-x] = 0
                         elif pol_pos[numpy.size(pol_pos)-1-x] <= pol_pos[numpy.size(pol_pos)-2-x] and pol_pos[numpy.size(pol_pos)-1-x] != 0 and pol_state[numpy.size(pol_pos)-2-x] == 0:
                             pol_state[numpy.size(pol_pos)-2-x] = 0
                             pol_state[numpy.size(pol_pos)-1-x] = 0
       elif (pol_pos[i] >= 1) and (pol_pos[i] < genome_chunks) and (pol_state[i] == 1) and (1 == numpy.logical_and(GE_chunks > pol_pos[i], GE_chunks < pol_pos[i] + k_transc).any()): #big block #7  #######################################
           d = random.choices(['yes','no'], [transc_term_prob[list(GE_tobe).index(1)], 1-transc_term_prob[list(GE_tobe).index(1)]])
           if d == ['yes']:
                pol_state[i] = 0
                print('transcription_TERMINATED @', GE_list[list(GE_tobe).index(1)])
                transc_term_events[list(GE_tobe).index(1)] = transc_term_events[list(GE_tobe).index(1)] + 1
                pol_pos[i] = GE_chunks[list(GE_tobe).index(1)]
                for n in range(0, checks):
                     for x in range(0, checks):
                         if pol_pos[numpy.size(pol_pos)-1-x] <= pol_pos[numpy.size(pol_pos)-2-x] and pol_pos[numpy.size(pol_pos)-1-x] != 0 and pol_state[numpy.size(pol_pos)-2-x] == 1 and pol_state[numpy.size(pol_pos)-1-x] == 0:
                             print('pol ejected from pos', pol_pos[numpy.size(pol_pos)-1-x])
                             pol_ej = pol_ej +1
                             pol_pos[numpy.size(pol_pos)-1-x] = pol_pos[numpy.size(pol_pos)-2-x]
                             pol_pos[numpy.size(pol_pos)-2-x] = 0
                             pol_state[numpy.size(pol_pos)-2-x] = 0
                             pol_state[numpy.size(pol_pos)-1-x] = 1
                         elif pol_pos[numpy.size(pol_pos)-1-x] < pol_pos[numpy.size(pol_pos)-2-x] and pol_pos[numpy.size(pol_pos)-1-x] != 0 and pol_state[numpy.size(pol_pos)-1-x] == 0 and pol_state[numpy.size(pol_pos)-2-x] == 0:
                             pol_pos[numpy.size(pol_pos)-2-x:numpy.size(pol_pos)-x].sort()
                             pol_state[numpy.size(pol_pos)-2-x] = 0
                             pol_state[numpy.size(pol_pos)-1-x] = 0
                         elif pol_pos[numpy.size(pol_pos)-1-x] <= pol_pos[numpy.size(pol_pos)-2-x] and pol_pos[numpy.size(pol_pos)-1-x] != 0 and pol_state[numpy.size(pol_pos)-2-x] == 0:
                             pol_state[numpy.size(pol_pos)-2-x] = 0
                             pol_state[numpy.size(pol_pos)-1-x] = 0
           elif d == ['no']:
                pol_state[i] = 1
                print('transc RT @', GE_list[list(GE_tobe).index(1)])
                transc_RT_events[list(GE_tobe).index(1)] = transc_RT_events[list(GE_tobe).index(1)] + 1
                pol_pos[i] = pol_pos[i] + k_transc
                for n in range(0, checks):
                     for x in range(0, checks):
                         if pol_pos[numpy.size(pol_pos)-1-x] <= pol_pos[numpy.size(pol_pos)-2-x] and pol_pos[numpy.size(pol_pos)-1-x] != 0 and pol_state[numpy.size(pol_pos)-2-x] == 1 and pol_state[numpy.size(pol_pos)-1-x] == 0:
                             print('pol ejected from pos', pol_pos[numpy.size(pol_pos)-1-x])
                             pol_ej = pol_ej +1
                             pol_pos[numpy.size(pol_pos)-1-x] = pol_pos[numpy.size(pol_pos)-2-x]
                             pol_pos[numpy.size(pol_pos)-2-x] = 0
                             pol_state[numpy.size(pol_pos)-2-x] = 0
                             pol_state[numpy.size(pol_pos)-1-x] = 1
                         elif pol_pos[numpy.size(pol_pos)-1-x] < pol_pos[numpy.size(pol_pos)-2-x] and pol_pos[numpy.size(pol_pos)-1-x] != 0 and pol_state[numpy.size(pol_pos)-1-x] == 0 and pol_state[numpy.size(pol_pos)-2-x] == 0:
                             pol_pos[numpy.size(pol_pos)-2-x:numpy.size(pol_pos)-x].sort()
                             pol_state[numpy.size(pol_pos)-2-x] = 0
                             pol_state[numpy.size(pol_pos)-1-x] = 0
                         elif pol_pos[numpy.size(pol_pos)-1-x] <= pol_pos[numpy.size(pol_pos)-2-x] and pol_pos[numpy.size(pol_pos)-1-x] != 0 and pol_state[numpy.size(pol_pos)-2-x] == 0:
                             pol_state[numpy.size(pol_pos)-2-x] = 0
                             pol_state[numpy.size(pol_pos)-1-x] = 0

    for n in range(0, checks):  #arranges values of pol_pos array from lowest to highest and keeps pol_state values with associated pol_pos values
        for x in range(1, checks+1):
            if pol_pos[x] == 0 and pol_pos[x-1] != 0:
                pol_pos[x] = pol_pos[x-1]
                pol_pos[x-1] = 0
                pol_state[x] = pol_state[x-1]
                pol_state[x-1] = 0

    print()
    print('pol_pos is', numpy.round(pol_pos,2))
    print('pol_state is', pol_state)
    print('# of pols bound is', numpy.count_nonzero(pol_pos))
    print()

print()
print('transc_events =', transc_events)
print('total transc_events =', numpy.sum(transc_events))
print('transc_term_events =', transc_term_events)
print('Diff =', numpy.subtract(transc_events, transc_term_events))
print('transc_RT_events =', transc_RT_events)
print('# pols ejected =', pol_ej)
print()
N_transc = transc_events[0]/(numpy.add(numpy.sum(transc_events),numpy.sum(transc_RT_events)))
print(N_transc)
P_transc = (transc_events[1] + transc_RT_events[0])/(numpy.add(numpy.sum(transc_events),numpy.sum(transc_RT_events)))
print(P_transc)
M_transc = (transc_events[2] + transc_RT_events[1] + transc_RT_events[0]*(1-transc_term_prob[1]))/(numpy.add(numpy.sum(transc_events),numpy.sum(transc_RT_events)))
print(M_transc)
G_transc = (transc_events[3] + transc_RT_events[2] + transc_RT_events[1]*(1-transc_term_prob[2]) + transc_RT_events[0]*(1-transc_term_prob[1])*(1-transc_term_prob[2]))/(numpy.add(numpy.sum(transc_events),numpy.sum(transc_RT_events)))
print(G_transc)
L_transc = (transc_events[4] + transc_RT_events[3] + transc_RT_events[2]*(1-transc_term_prob[3]) + transc_RT_events[1]*(1-transc_term_prob[2])*(1-transc_term_prob[3]) + transc_RT_events[0]*(1-transc_term_prob[1])*(1-transc_term_prob[2])*(1-transc_term_prob[3]))/(numpy.add(numpy.sum(transc_events),numpy.sum(transc_RT_events)))
print(L_transc)
print()
print('error=', (((abs(N_transc-0.361))+(abs(P_transc-0.252))+(abs(M_transc-0.177))+(abs(G_transc-0.124))+(abs(L_transc-0.087)))/5))
