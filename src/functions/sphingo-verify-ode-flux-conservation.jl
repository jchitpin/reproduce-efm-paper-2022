# ODEs for wildtype condition
function model_wt!(du, u, p, t)

  #Km1 = nothing
  Km2 = 0.081
  Km3 = 0.171
  Km4 = 0.0034
  Km5 = 0.0385
  Km6 = 0.035
  #Km7 = nothing
  Km8 = 0.155
  Km9 = 0.126
  Km10 = 0.0601
  Km11 = 0.0034
  Km12 = 0.025
  #Km13 = nothing
  #Km14 = nothing
  #Km15 = nothing
  #Km16 = nothing
  Km17 = 0.0385
  Km18 = 0.0034
  Km19 = 0.0025
  Km20 = 0.149
  #Km21 = nothing
  #Km22 = nothing
  #Km23 = nothing
  #Km24 = nothing
  Km25 = 0.0455
  Km26 = 0.0601
  Km27 = 0.034
  Km28 = 0.148
  Km29 = 0.0601
  Km30 = 0.00506
  #Km31 = nothing
  #Km32 = nothing
  #Km33 = nothing
  #Km34 = nothing
  #Km35 = nothing
  Km36 = 0.0034
  #Km37 = nothing
  #Km38 = nothing
  #Km39 = nothing
  #Km40 = nothing
  Km41 = 0.04
  #Km42 = nothing
  Km43 = 0.003
  Km44 = 0.003
  #Km45 = nothing
  #Km46 = nothing
  Km47 = 0.019
  #Km48 = nothing
  Km49 = 0.0455
  Km50 = 0.149
  #Km51 = nothing
  #Km52 = nothing
  Km53 = 0.00506
  Km54 = 0.036
  Km55 = 0.02
  #Km56 = nothing
  Km57 = 0.148
  Km58 = 0.081
  Km59 = 0.107
  Km60 = 0.036
  #Km61 = nothing
  Km62 = 0.036
  Km63 = 0.107
  #Km64 = nothing
  #Km65 = nothing
  #Km66 = nothing
  #Km67 = nothing
  #Km68 = nothing
  Km69 = 0.036

  Vm1 = 0
  Vm2 = 0.48
  Vm3 = 2.4
  Vm4 = 0.175
  Vm5 = 3.1
  Vm6 = 0
  #Vm7 = nothing
  Vm8 = 0.109
  Vm9 = 0.00113
  Vm10 = 0.068
  Vm11 = 0.1
  Vm12 = 1.24
  #Vm13 = nothing
  #Vm14 = nothing
  #Vm15 = nothing
  #Vm16 = nothing
  Vm17 = 3
  Vm18 = 0.15
  Vm19 = 0.1
  Vm20 = 2.27
  #Vm21 = nothing
  #Vm22 = nothing
  #Vm23 = nothing
  #Vm24 = nothing
  Vm25 = 0.001
  Vm26 = 0.12
  Vm27 = 0.11
  Vm28 = 0.002
  Vm29 = 0.03
  Vm30 = 0.003
  #Vm31 = nothing
  #Vm32 = nothing
  #Vm33 = nothing
  #Vm34 = nothing
  #Vm35 = nothing
  Vm36 = 0.0036
  #Vm37 = nothing
  #Vm38 = nothing
  #Vm39 = nothing
  #Vm40 = nothing
  Vm41 = 0.01
  #Vm42 = nothing
  Vm43 = 0.001
  Vm44 = 0.001
  #Vm45 = nothing
  #Vm46 = nothing
  Vm47 = 0.0061
  #Vm48 = nothing
  Vm49 = 0.00183
  Vm50 = 0.1
  #Vm51 = nothing
  #Vm52 = nothing
  Vm53 = 0.13
  Vm54 = 2
  Vm55 = 0.3
  #Vm56 = nothing
  Vm57 = 0.05
  Vm58 = 0.8
  Vm59 = 2
  Vm60 = 0.75
  #Vm61 = nothing
  Vm62 = 0.78
  Vm63 = 2.1
  #Vm64 = nothing
  #Vm65 = nothing
  #Vm66 = nothing
  #Vm67 = nothing
  #Vm68 = nothing
  Vm69 = 0.43

  #k1 = nothing
  #k2 = nothing
  #k3 = nothing
  #k4 = nothing
  #k5 = nothing
  #k6 = nothing
  k7 = 0.8
  #k8 = nothing
  #k9 = nothing
  #k10 = nothing
  #k11 = nothing
  #k12 = nothing
  k13 = 0.12
  k14 = 0.5
  k15 = 0.44
  k16 = 0.25
  #k17 = nothing
  #k18 = nothing
  #k19 = nothing
  #k20 = nothing
  k21 = 0.43
  k22 = 23
  k23 = 0.01
  k24 = 0.006
  #k25 = nothing
  #k26 = nothing
  #k27 = nothing
  #k28 = nothing
  #k29 = nothing
  #k30 = nothing
  k31 = 1
  k32 = 1
  k33 = 1
  k34 = 0.4
  k35 = 3
  #k36 = nothing
  k37 = 4.5
  k38 = 1
  k39 = 5
  k40 = 0.08
  #k41 = nothing
  k42 = 1
  #k43 = nothing
  #k44 = nothing
  k45 = 0.2
  k46 = 0.03
  #k47 = nothing
  k48 = 0.002
  #k49 = nothing
  #k50 = nothing
  k51 = 1
  k52 = 0.2
  #k53 = nothing
  #k54 = nothing
  #k55 = nothing
  k56 = 0.015
  #k57 = nothing
  #k58 = nothing
  #k59 = nothing
  #k60 = nothing
  k61 = 2.5
  #k62 = nothing
  #k63 = nothing
  k64 = 0
  k65 = 0
  k66 = 0
  k67 = 0
  k68 = 0
  #k69 = nothing

  #r1 = nothing
  #r2 = nothing
  #r3 = nothing
  #r4 = nothing
  #r5 = nothing
  #r6 = nothing
  r7 = 0.1
  #r8 = nothing
  #r9 = nothing
  #r10 = nothing
  #r11 = nothing
  #r12 = nothing
  r13 = 0.001
  r14 = 0.23
  r15 = 0.25
  r16 = 0.044
  #r17 = nothing
  #r18 = nothing
  #r19 = nothing
  #r20 = nothing
  r21 = 0.034
  r22 = 0.03
  r23 = 0.005
  r24 = 0.001
  #r25 = nothing
  #r26 = nothing
  #r27 = nothing
  #r28 = nothing
  #r29 = nothing
  #r30 = nothing
  r31 = 1
  r32 = 3
  #r33 = nothing
  r34 = 0.25
  r35 = 0.0269
  #r36 = nothing
  r37 = 0.005
  r38 = 0.75
  #r39 = nothing
  #r40 = nothing
  #r41 = nothing
  #r42 = nothing
  #r43 = nothing
  #r44 = nothing
  #r45 = nothing
  #r46 = nothing
  #r47 = nothing
  #r48 = nothing
  #r49 = nothing
  #r50 = nothing
  #r51 = nothing
  r52 = 9
  #r53 = nothing
  #r54 = nothing
  #r55 = nothing
  #r56 = nothing
  #r57 = nothing
  #r58 = nothing
  #r59 = nothing
  #r60 = nothing
  #r61 = nothing
  #r62 = nothing
  #r63 = nothing
  #r64 = nothing
  r65 = 0
  #r66 = nothing
  r67 = 0
  r68 = 0
  #r69 = nothing

  # Inhibition (no values provided in paper. These values have been reverse-engineered)
  Ki_s1p = 0.0356898872566651
  Ki_c1p = 0.178449436283326

  du[1] = ( (Vm59 * u[6]) / (Km59 + u[6]) ) - ( (Vm60 * u[1]) / (Km60 + u[1]) ) - ( k61 * u[1] )
  du[2] = ( (k61 * u[1]) ) - ( (Vm62 * u[2]) / (Km62 + u[2]) ) + ( (Vm63 * u[10]) / (Km63 + u[10]) ) + ( (k65 - r65 * u[2]) )
  du[3] = ( (Vm1) / ((1 + u[17] / Ki_s1p) * (1 + u[1] / Ki_c1p)) ) - ( (Vm2 * u[3]) / (Km2 + u[3]) ) + ( (Vm3 * u[31]) / (Km3 + u[31]) ) - ( k7 * u[3] - r7 * u[9] ) - ( k38 * u[3] - r38 * u[8] ) - ( k39 * u[3] ) - ( k40 * u[3] )
  du[4] = ( (Vm28 * u[25]) / (Km28 + u[25]) ) - ( (Vm29 * u[4]) / (Km29 + u[4]) ) - ( k31 * u[4] - r31 * u[10] )
  du[5] = ( (Vm47 * u[13]) / (Km47 + u[13]) ) + ( (Vm49 * u[26]) / ((Km49 + u[26]) * (1 + u[1] / Ki_c1p) * (1 + u[17] / Ki_s1p)) ) - ( (Vm50 * u[5]) / (Km50 + u[5]) )
  du[6] = ( k39 * u[3] ) - ( (Vm55 * u[6]) / (Km55 + u[6]) ) + ( (Vm57 * u[27]) / (Km57 + u[27]) ) - ( (Vm58 * u[6]) / (Km58 + u[6]) ) - ( (Vm59 * u[6]) / (Km59 + u[6]) ) + ( (Vm60 * u[1]) / (Km60 + u[1]) )
  du[7] = (k40 * u[3]) - ( (Vm41 * u[7]) / (Km41 + u[7]) )
  du[8] = ( (Vm19 * u[35]) / (Km19 + u[35]) ) - ( (Vm20 * u[8]) / (Km20 + u[8]) ) + ( k38 * u[3] - r38 * u[8] )
  du[9] = ( k7 * u[3] - r7 * u[9] ) - ( (Vm8 * u[9]) / (Km8 + u[9]) ) + ( (Vm9 * u[28]) / (Km9 + u[28]) ) - ( (Vm10 * u[9]) / (Km10 + u[9]) )
  du[10] = ( (Vm25 * u[29]) / (Km25 + u[29]) ) - ( (Vm26 * u[10]) / (Km26 + u[10]) ) + ( k31 * u[4] - r31 * u[10] ) + ( (Vm62 * u[2]) / (Km62 + u[2]) ) - ( (Vm63 * u[10]) / (Km63 + u[10]) ) + (k66)
  du[11] =  ( k42 * u[12] ) - ( (Vm43 * u[11]) / (Km43 + u[11]) )
  du[12] = ( (Vm41 * u[7]) / (Km41 + u[7]) ) - ( k42 * u[12] )
  du[13] = ( k46 * u[15] ) - ( (Vm47 * u[13]) / (Km47 + u[13]) )
  du[14] = ( (Vm44 * u[16]) / (Km44 + u[16]) ) - ( k45 * u[14] )
  du[15] = ( k45 * u[14] ) - ( k46 * u[15])
  du[16] = ( (Vm43 * u[11]) / (Km43 + u[11]) ) - ( (Vm44 * u[16]) / (Km44 + u[16]) )
  du[17] = ( k15 * u[22] - r15 * u[17] ) - ( k16 * u[17] - r16 * u[21] ) + ( k34 * u[19] - r34 * u[17] ) + ( (Vm36 * u[30]) / (Km36 + u[30]) ) - ( k37 * u[17] - r37 * u[18] )
  du[18] = ( (Vm4 * u[31]) / (Km4 + u[31]) ) - ( (Vm5 * u[18]) / (Km5 + u[18]) ) - ( (Vm6 * u[18]) / (Km6 + u[18]) )  + ( k37 * u[17] - r37 * u[18])
  du[19] = ( (Vm30 * u[32]) / (Km30 + u[32]) ) - ( k33 * u[19] ) - ( k34 * u[19] - r34 * u[17] )
  du[20] = ( (Vm53 * u[34]) / (Km53 + u[34]) ) - ( (Vm54 * u[20]) / (Km54 + u[20]) )
  du[21] =  ( k16 * u[17] - r16 * u[21] ) - ( (Vm17 * u[21]) / (Km17 + u[21]) ) + ( (Vm18 * u[35]) / (Km18 + u[35]) )
  du[22] = ( (Vm11 * u[36]) / (Km11 + u[36]) ) - ( (Vm12 * u[22]) / (Km12 + u[22]) ) - ( k15 * u[22] - r15 * u[17] )
  du[23] = ( (Vm27 * u[37]) / (Km27 + u[37]) ) + ( k33 * u[19] ) + ( k68 - r68 * u[23] ) - ( (Vm69 * u[23]) / (Km69 + u[23]) )
  du[24] = ( k13 * u[28] - r13 * u[24] ) - ( k23 * u[24] - r23 * u[25] ) - ( k24 * u[24] - r24 * u[29] )
  du[25] = ( k23 * u[24] - r23 * u[25] ) - ( (Vm28 * u[25]) / (Km28 + u[25]) )
  du[26] = ( k48 * u[29] ) - ( (Vm49 * u[26]) / ((Km49 + u[26]) * (1 + u[1] / Ki_c1p) * (1 + u[17] / Ki_s1p)) )
  du[27] = ( (Vm55 * u[6]) / (Km55 + u[6]) ) - ( k56 * u[27] ) - ( (Vm57 * u[27]) / (Km57 + u[27]) )
  du[28] = ( (Vm8 * u[9]) / (Km8 + u[9]) ) - ( (Vm9 * u[28]) / (Km9 + u[28]) ) - ( k13 * u[28] - r13 * u[24] )
  du[29] = ( k24 * u[24] - r24 * u[29] ) - ( (Vm25 * u[29]) / (Km25 + u[29]) ) - ( k48 * u[29] ) + ( k56 * u[27] ) + ( k64 )
  du[30] = ( k14 * u[36] - r14 * u[30] ) + ( k21 * u[35] - r21 * u[30] ) - ( k22 * u[30] - r22 * u[31] ) + ( k35 * u[32] - r35 * u[30] ) - ( (Vm36 * u[30]) / (Km36 + u[30]) ) + ( k51 * u[33] ) - ( k52 * u[30] - r52 * u[34] )
  du[31] = ( (Vm2 * u[3]) / (Km2 + u[3]) ) - ( (Vm3 * u[31]) / (Km3 + u[31]) ) - ( (Vm4 * u[31]) / (Km4 + u[31]) ) + ( (Vm5 * u[18]) / (Km5 + u[18]) ) + ( k22 * u[30] - r22 * u[31])
  du[32] = ( (Vm29 * u[4]) / (Km29 + u[4]) ) - ( (Vm30 * u[32]) / (Km30 + u[32]) ) - ( k32 * u[32] - r32 * u[37]) - ( k35 * u[32] - r35 * u[30])
  du[33] = ( (Vm50 * u[5]) / (Km50 + u[5]) ) - ( k51 * u[33])
  du[34] = ( k52 * u[30] - r52 * u[34] ) - ( (Vm53 * u[34]) / (Km53 + u[34]) ) + ( (Vm54 * u[20]) / (Km54 + u[20]) ) + ( (Vm58 * u[6]) / (Km58 + u[6]) )
  du[35] = ( (Vm17 * u[21]) / (Km17 + u[21]) ) - ( (Vm18 * u[35]) / (Km18 + u[35]) ) - ( (Vm19 * u[35]) / (Km19 + u[35]) ) + ( (Vm20 * u[8]) / (Km20 + u[8]) ) - ( k21 * u[35] - r21 * u[30] )
  du[36] = ( (Vm10 * u[9]) / (Km10 + u[9]) ) - ( (Vm11 * u[36]) / (Km11 + u[36]) ) + ( (Vm12 * u[22]) / (Km12 + u[22]) ) - ( k14 * u[36] - r14 * u[30])
  du[37] = ( (Vm26 * u[10]) / (Km26 + u[10]) ) - ( (Vm27 * u[37]) / (Km27 + u[37]) ) + ( k32 * u[32] - r32 * u[37] ) + ( k67 - r67 * u[37]) + ( (Vm69 * u[23]) / (Km69 + u[23]) )
end

# ODEs for disease condition
function model_ad!(du, u, p, t)

  #Km1 = nothing
  Km2 = 0.081
  Km3 = 0.171
  Km4 = 0.0034
  Km5 = 0.0385
  Km6 = 0.035
  #Km7 = nothing
  Km8 = 0.155
  Km9 = 0.126
  Km10 = 0.0601
  Km11 = 0.0034
  Km12 = 0.025
  #Km13 = nothing
  #Km14 = nothing
  #Km15 = nothing
  #Km16 = nothing
  Km17 = 0.0385
  Km18 = 0.0034
  Km19 = 0.0025
  Km20 = 0.149
  #Km21 = nothing
  #Km22 = nothing
  #Km23 = nothing
  #Km24 = nothing
  Km25 = 0.0455
  Km26 = 0.0601
  Km27 = 0.034
  Km28 = 0.148
  Km29 = 0.0601
  Km30 = 0.00506
  #Km31 = nothing
  #Km32 = nothing
  #Km33 = nothing
  #Km34 = nothing
  #Km35 = nothing
  Km36 = 0.0034
  #Km37 = nothing
  #Km38 = nothing
  #Km39 = nothing
  #Km40 = nothing
  Km41 = 0.04
  #Km42 = nothing
  Km43 = 0.003
  Km44 = 0.003
  #Km45 = nothing
  #Km46 = nothing
  Km47 = 0.019
  #Km48 = nothing
  Km49 = 0.0455
  Km50 = 0.149
  #Km51 = nothing
  #Km52 = nothing
  Km53 = 0.00506
  Km54 = 0.036
  Km55 = 0.02
  #Km56 = nothing
  Km57 = 0.148
  Km58 = 0.081
  Km59 = 0.107
  Km60 = 0.036
  #Km61 = nothing
  Km62 = 0.036
  Km63 = 0.107
  #Km64 = nothing
  #Km65 = nothing
  #Km66 = nothing
  #Km67 = nothing
  #Km68 = nothing
  Km69 = 0.036

  #Vm1 = 0.004
  Vm1 = 0
  Vm2 = 0.32
  Vm3 = 2.4
  Vm4 = 0.875
  Vm5 = 3.1
  Vm6 = 0
  #Vm7 = nothing
  Vm8 = 0.109
  Vm9 = 0.1
  Vm10 = 0.00453
  Vm11 = 0.05
  Vm12 = 1.24
  #Vm13 = nothing
  #Vm14 = nothing
  #Vm15 = nothing
  #Vm16 = nothing
  Vm17 = 3
  Vm18 = 0.07
  Vm19 = 0.1
  Vm20 = 2
  #Vm21 = nothing
  #Vm22 = nothing
  #Vm23 = nothing
  #Vm24 = nothing
  Vm25 = 0.002
  Vm26 = 0.08
  Vm27 = 0.055
  Vm28 = 0.0025
  Vm29 = 0.02
  Vm30 = 0.0015
  #Vm31 = nothing
  #Vm32 = nothing
  #Vm33 = nothing
  #Vm34 = nothing
  #Vm35 = nothing
  Vm36 = 0.001
  #Vm37 = nothing
  #Vm38 = nothing
  #Vm39 = nothing
  #Vm40 = nothing
  Vm41 = 0.01
  #Vm42 = nothing
  Vm43 = 0.001
  Vm44 = 0.001
  #Vm45 = nothing
  #Vm46 = nothing
  Vm47 = 0.0061
  #Vm48 = nothing
  Vm49 = 0.01
  Vm50 = 0.00669
  #Vm51 = nothing
  #Vm52 = nothing
  Vm53 = 0.08
  Vm54 = 2
  Vm55 = 0.3
  #Vm56 = nothing
  Vm57 = 0.07
  Vm58 = 0.529
  Vm59 = 5
  Vm60 = 0.75
  #Vm61 = nothing
  Vm62 = 0.78
  Vm63 = 0.5
  #Vm64 = nothing
  #Vm65 = nothing
  #Vm66 = nothing
  #Vm67 = nothing
  #Vm68 = nothing
  Vm69 = 0.43

  #k1 = nothing
  #k2 = nothing
  #k3 = nothing
  #k4 = nothing
  #k5 = nothing
  #k6 = nothing
  k7 = 0.8
  #k8 = nothing
  #k9 = nothing
  #k10 = nothing
  #k11 = nothing
  #k12 = nothing
  k13 = 0.12
  k14 = 0.5
  k15 = 0.44
  k16 = 0.25
  #k17 = nothing
  #k18 = nothing
  #k19 = nothing
  #k20 = nothing
  k21 = 0.43
  k22 = 23
  k23 = 0.01
  k24 = 0.015
  #k25 = nothing
  #k26 = nothing
  #k27 = nothing
  #k28 = nothing
  #k29 = nothing
  #k30 = nothing
  k31 = 1
  k32 = 1
  k33 = 1
  k34 = 0.4
  k35 = 3
  #k36 = nothing
  k37 = 4.5
  k38 = 1
  k39 = 1
  k40 = 0.02
  #k41 = nothing
  k42 = 1
  #k43 = nothing
  #k44 = nothing
  k45 = 0.2
  k46 = 0.03
  #k47 = nothing
  k48 = 0.011
  #k49 = nothing
  #k50 = nothing
  k51 = 1
  k52 = 0.2
  #k53 = nothing
  #k54 = nothing
  #k55 = nothing
  k56 = 0.045
  #k57 = nothing
  #k58 = nothing
  #k59 = nothing
  #k60 = nothing
  k61 = 2.5
  #k62 = nothing
  #k63 = nothing
  k64 = 0
  k65 = 0
  k66 = 0
  k67 = 0
  k68 = 0
  #k69 = nothing

  #r1 = nothing
  #r2 = nothing
  #r3 = nothing
  #r4 = nothing
  #r5 = nothing
  #r6 = nothing
  r7 = 0.1
  #r8 = nothing
  #r9 = nothing
  #r10 = nothing
  #r11 = nothing
  #r12 = nothing
  r13 = 0.001
  r14 = 0.23
  r15 = 0.25
  r16 = 0.044
  #r17 = nothing
  #r18 = nothing
  #r19 = nothing
  #r20 = nothing
  r21 = 0.034
  r22 = 0.03
  r23 = 0.005
  r24 = 0.001
  #r25 = nothing
  #r26 = nothing
  #r27 = nothing
  #r28 = nothing
  #r29 = nothing
  #r30 = nothing
  r31 = 1
  r32 = 3
  #r33 = nothing
  r34 = 0.25
  r35 = 0.0269
  #r36 = nothing
  r37 = 0.005
  r38 = 0.75
  #r39 = nothing
  #r40 = nothing
  #r41 = nothing
  #r42 = nothing
  #r43 = nothing
  #r44 = nothing
  #r45 = nothing
  #r46 = nothing
  #r47 = nothing
  #r48 = nothing
  #r49 = nothing
  #r50 = nothing
  #r51 = nothing
  r52 = 9
  #r53 = nothing
  #r54 = nothing
  #r55 = nothing
  #r56 = nothing
  #r57 = nothing
  #r58 = nothing
  #r59 = nothing
  #r60 = nothing
  #r61 = nothing
  #r62 = nothing
  #r63 = nothing
  #r64 = nothing
  r65 = 0
  #r66 = nothing
  r67 = 0
  r68 = 0
  #r69 = nothing

  # Inhibition (no values provided in paper. These values have been reverse-engineered)
  Ki_s1p = 0.0356898872566651
  Ki_c1p = 0.178449436283326

  du[1] = ( (Vm59 * u[6]) / (Km59 + u[6]) ) - ( (Vm60 * u[1]) / (Km60 + u[1]) ) - ( k61 * u[1] )
  du[2] = ( (k61 * u[1]) ) - ( (Vm62 * u[2]) / (Km62 + u[2]) ) + ( (Vm63 * u[10]) / (Km63 + u[10]) ) + ( (k65 - r65 * u[2]) )
  du[3] = ( (Vm1) / ((1 + u[17] / Ki_s1p) * (1 + u[1] / Ki_c1p)) ) - ( (Vm2 * u[3]) / (Km2 + u[3]) ) + ( (Vm3 * u[31]) / (Km3 + u[31]) ) - ( k7 * u[3] - r7 * u[9] ) - ( k38 * u[3] - r38 * u[8] ) - ( k39 * u[3] ) - ( k40 * u[3] )
  du[4] = ( (Vm28 * u[25]) / (Km28 + u[25]) ) - ( (Vm29 * u[4]) / (Km29 + u[4]) ) - ( k31 * u[4] - r31 * u[10] )
  du[5] = ( (Vm47 * u[13]) / (Km47 + u[13]) ) + ( (Vm49 * u[26]) / ((Km49 + u[26]) * (1 + u[1] / Ki_c1p) * (1 + u[17] / Ki_s1p)) ) - ( (Vm50 * u[5]) / (Km50 + u[5]) )
  du[6] = ( k39 * u[3] ) - ( (Vm55 * u[6]) / (Km55 + u[6]) ) + ( (Vm57 * u[27]) / (Km57 + u[27]) ) - ( (Vm58 * u[6]) / (Km58 + u[6]) ) - ( (Vm59 * u[6]) / (Km59 + u[6]) ) + ( (Vm60 * u[1]) / (Km60 + u[1]) )
  du[7] = (k40 * u[3]) - ( (Vm41 * u[7]) / (Km41 + u[7]) )
  du[8] = ( (Vm19 * u[35]) / (Km19 + u[35]) ) - ( (Vm20 * u[8]) / (Km20 + u[8]) ) + ( k38 * u[3] - r38 * u[8] )
  du[9] = ( k7 * u[3] - r7 * u[9] ) - ( (Vm8 * u[9]) / (Km8 + u[9]) ) + ( (Vm9 * u[28]) / (Km9 + u[28]) ) - ( (Vm10 * u[9]) / (Km10 + u[9]) )
  du[10] = ( (Vm25 * u[29]) / (Km25 + u[29]) ) - ( (Vm26 * u[10]) / (Km26 + u[10]) ) + ( k31 * u[4] - r31 * u[10] ) + ( (Vm62 * u[2]) / (Km62 + u[2]) ) - ( (Vm63 * u[10]) / (Km63 + u[10]) ) + (k66)
  du[11] =  ( k42 * u[12] ) - ( (Vm43 * u[11]) / (Km43 + u[11]) )
  du[12] = ( (Vm41 * u[7]) / (Km41 + u[7]) ) - ( k42 * u[12] )
  du[13] = ( k46 * u[15] ) - ( (Vm47 * u[13]) / (Km47 + u[13]) )
  du[14] = ( (Vm44 * u[16]) / (Km44 + u[16]) ) - ( k45 * u[14] )
  du[15] = ( k45 * u[14] ) - ( k46 * u[15])
  du[16] = ( (Vm43 * u[11]) / (Km43 + u[11]) ) - ( (Vm44 * u[16]) / (Km44 + u[16]) )
  du[17] = ( k15 * u[22] - r15 * u[17] ) - ( k16 * u[17] - r16 * u[21] ) + ( k34 * u[19] - r34 * u[17] ) + ( (Vm36 * u[30]) / (Km36 + u[30]) ) - ( k37 * u[17] - r37 * u[18] )
  du[18] = ( (Vm4 * u[31]) / (Km4 + u[31]) ) - ( (Vm5 * u[18]) / (Km5 + u[18]) ) - ( (Vm6 * u[18]) / (Km6 + u[18]) )  + ( k37 * u[17] - r37 * u[18])
  du[19] = ( (Vm30 * u[32]) / (Km30 + u[32]) ) - ( k33 * u[19] ) - ( k34 * u[19] - r34 * u[17] )
  du[20] = ( (Vm53 * u[34]) / (Km53 + u[34]) ) - ( (Vm54 * u[20]) / (Km54 + u[20]) )
  du[21] =  ( k16 * u[17] - r16 * u[21] ) - ( (Vm17 * u[21]) / (Km17 + u[21]) ) + ( (Vm18 * u[35]) / (Km18 + u[35]) )
  du[22] = ( (Vm11 * u[36]) / (Km11 + u[36]) ) - ( (Vm12 * u[22]) / (Km12 + u[22]) ) - ( k15 * u[22] - r15 * u[17] )
  du[23] = ( (Vm27 * u[37]) / (Km27 + u[37]) ) + ( k33 * u[19] ) + ( k68 - r68 * u[23] ) - ( (Vm69 * u[23]) / (Km69 + u[23]) )
  du[24] = ( k13 * u[28] - r13 * u[24] ) - ( k23 * u[24] - r23 * u[25] ) - ( k24 * u[24] - r24 * u[29] )
  du[25] = ( k23 * u[24] - r23 * u[25] ) - ( (Vm28 * u[25]) / (Km28 + u[25]) )
  du[26] = ( k48 * u[29] ) - ( (Vm49 * u[26]) / ((Km49 + u[26]) * (1 + u[1] / Ki_c1p) * (1 + u[17] / Ki_s1p)) )
  du[27] = ( (Vm55 * u[6]) / (Km55 + u[6]) ) - ( k56 * u[27] ) - ( (Vm57 * u[27]) / (Km57 + u[27]) )
  du[28] = ( (Vm8 * u[9]) / (Km8 + u[9]) ) - ( (Vm9 * u[28]) / (Km9 + u[28]) ) - ( k13 * u[28] - r13 * u[24] )
  du[29] = ( k24 * u[24] - r24 * u[29] ) - ( (Vm25 * u[29]) / (Km25 + u[29]) ) - ( k48 * u[29] ) + ( k56 * u[27] ) + ( k64 )
  du[30] = ( k14 * u[36] - r14 * u[30] ) + ( k21 * u[35] - r21 * u[30] ) - ( k22 * u[30] - r22 * u[31] ) + ( k35 * u[32] - r35 * u[30] ) - ( (Vm36 * u[30]) / (Km36 + u[30]) ) + ( k51 * u[33] ) - ( k52 * u[30] - r52 * u[34] )
  du[31] = ( (Vm2 * u[3]) / (Km2 + u[3]) ) - ( (Vm3 * u[31]) / (Km3 + u[31]) ) - ( (Vm4 * u[31]) / (Km4 + u[31]) ) + ( (Vm5 * u[18]) / (Km5 + u[18]) ) + ( k22 * u[30] - r22 * u[31])
  du[32] = ( (Vm29 * u[4]) / (Km29 + u[4]) ) - ( (Vm30 * u[32]) / (Km30 + u[32]) ) - ( k32 * u[32] - r32 * u[37]) - ( k35 * u[32] - r35 * u[30])
  du[33] = ( (Vm50 * u[5]) / (Km50 + u[5]) ) - ( k51 * u[33])
  du[34] = ( k52 * u[30] - r52 * u[34] ) - ( (Vm53 * u[34]) / (Km53 + u[34]) ) + ( (Vm54 * u[20]) / (Km54 + u[20]) ) + ( (Vm58 * u[6]) / (Km58 + u[6]) )
  du[35] = ( (Vm17 * u[21]) / (Km17 + u[21]) ) - ( (Vm18 * u[35]) / (Km18 + u[35]) ) - ( (Vm19 * u[35]) / (Km19 + u[35]) ) + ( (Vm20 * u[8]) / (Km20 + u[8]) ) - ( k21 * u[35] - r21 * u[30] )
  du[36] = ( (Vm10 * u[9]) / (Km10 + u[9]) ) - ( (Vm11 * u[36]) / (Km11 + u[36]) ) + ( (Vm12 * u[22]) / (Km12 + u[22]) ) - ( k14 * u[36] - r14 * u[30])
  du[37] = ( (Vm26 * u[10]) / (Km26 + u[10]) ) - ( (Vm27 * u[37]) / (Km27 + u[37]) ) + ( k32 * u[32] - r32 * u[37] ) + ( k67 - r67 * u[37]) + ( (Vm69 * u[23]) / (Km69 + u[23]) )
end

# Fluxes for wildtype condition
function model_wt_fluxes(u::Vector{Float64})

  #Km1 = nothing
  Km2 = 0.081
  Km3 = 0.171
  Km4 = 0.0034
  Km5 = 0.0385
  Km6 = 0.035
  #Km7 = nothing
  Km8 = 0.155
  Km9 = 0.126
  Km10 = 0.0601
  Km11 = 0.0034
  Km12 = 0.025
  #Km13 = nothing
  #Km14 = nothing
  #Km15 = nothing
  #Km16 = nothing
  Km17 = 0.0385
  Km18 = 0.0034
  Km19 = 0.0025
  Km20 = 0.149
  #Km21 = nothing
  #Km22 = nothing
  #Km23 = nothing
  #Km24 = nothing
  Km25 = 0.0455
  Km26 = 0.0601
  Km27 = 0.034
  Km28 = 0.148
  Km29 = 0.0601
  Km30 = 0.00506
  #Km31 = nothing
  #Km32 = nothing
  #Km33 = nothing
  #Km34 = nothing
  #Km35 = nothing
  Km36 = 0.0034
  #Km37 = nothing
  #Km38 = nothing
  #Km39 = nothing
  #Km40 = nothing
  Km41 = 0.04
  #Km42 = nothing
  Km43 = 0.003
  Km44 = 0.003
  #Km45 = nothing
  #Km46 = nothing
  Km47 = 0.019
  #Km48 = nothing
  Km49 = 0.0455
  Km50 = 0.149
  #Km51 = nothing
  #Km52 = nothing
  Km53 = 0.00506
  Km54 = 0.036
  Km55 = 0.02
  #Km56 = nothing
  Km57 = 0.148
  Km58 = 0.081
  Km59 = 0.107
  Km60 = 0.036
  #Km61 = nothing
  Km62 = 0.036
  Km63 = 0.107
  #Km64 = nothing
  #Km65 = nothing
  #Km66 = nothing
  #Km67 = nothing
  #Km68 = nothing
  Km69 = 0.036

  Vm1 = 0
  Vm2 = 0.48
  Vm3 = 2.4
  Vm4 = 0.175
  Vm5 = 3.1
  Vm6 = 0
  #Vm7 = nothing
  Vm8 = 0.109
  Vm9 = 0.00113
  Vm10 = 0.068
  Vm11 = 0.1
  Vm12 = 1.24
  #Vm13 = nothing
  #Vm14 = nothing
  #Vm15 = nothing
  #Vm16 = nothing
  Vm17 = 3
  Vm18 = 0.15
  Vm19 = 0.1
  Vm20 = 2.27
  #Vm21 = nothing
  #Vm22 = nothing
  #Vm23 = nothing
  #Vm24 = nothing
  Vm25 = 0.001
  Vm26 = 0.12
  Vm27 = 0.11
  Vm28 = 0.002
  Vm29 = 0.03
  Vm30 = 0.003
  #Vm31 = nothing
  #Vm32 = nothing
  #Vm33 = nothing
  #Vm34 = nothing
  #Vm35 = nothing
  Vm36 = 0.0036
  #Vm37 = nothing
  #Vm38 = nothing
  #Vm39 = nothing
  #Vm40 = nothing
  Vm41 = 0.01
  #Vm42 = nothing
  Vm43 = 0.001
  Vm44 = 0.001
  #Vm45 = nothing
  #Vm46 = nothing
  Vm47 = 0.0061
  #Vm48 = nothing
  Vm49 = 0.00183
  Vm50 = 0.1
  #Vm51 = nothing
  #Vm52 = nothing
  Vm53 = 0.13
  Vm54 = 2
  Vm55 = 0.3
  #Vm56 = nothing
  Vm57 = 0.05
  Vm58 = 0.8
  Vm59 = 2
  Vm60 = 0.75
  #Vm61 = nothing
  Vm62 = 0.78
  Vm63 = 2.1
  #Vm64 = nothing
  #Vm65 = nothing
  #Vm66 = nothing
  #Vm67 = nothing
  #Vm68 = nothing
  Vm69 = 0.43

  #k1 = nothing
  #k2 = nothing
  #k3 = nothing
  #k4 = nothing
  #k5 = nothing
  #k6 = nothing
  k7 = 0.8
  #k8 = nothing
  #k9 = nothing
  #k10 = nothing
  #k11 = nothing
  #k12 = nothing
  k13 = 0.12
  k14 = 0.5
  k15 = 0.44
  k16 = 0.25
  #k17 = nothing
  #k18 = nothing
  #k19 = nothing
  #k20 = nothing
  k21 = 0.43
  k22 = 23
  k23 = 0.01
  k24 = 0.006
  #k25 = nothing
  #k26 = nothing
  #k27 = nothing
  #k28 = nothing
  #k29 = nothing
  #k30 = nothing
  k31 = 1
  k32 = 1
  k33 = 1
  k34 = 0.4
  k35 = 3
  #k36 = nothing
  k37 = 4.5
  k38 = 1
  k39 = 5
  k40 = 0.08
  #k41 = nothing
  k42 = 1
  #k43 = nothing
  #k44 = nothing
  k45 = 0.2
  k46 = 0.03
  #k47 = nothing
  k48 = 0.002
  #k49 = nothing
  #k50 = nothing
  k51 = 1
  k52 = 0.2
  #k53 = nothing
  #k54 = nothing
  #k55 = nothing
  k56 = 0.015
  #k57 = nothing
  #k58 = nothing
  #k59 = nothing
  #k60 = nothing
  k61 = 2.5
  #k62 = nothing
  #k63 = nothing
  k64 = 0
  k65 = 0
  k66 = 0
  k67 = 0
  k68 = 0
  #k69 = nothing

  #r1 = nothing
  #r2 = nothing
  #r3 = nothing
  #r4 = nothing
  #r5 = nothing
  #r6 = nothing
  r7 = 0.1
  #r8 = nothing
  #r9 = nothing
  #r10 = nothing
  #r11 = nothing
  #r12 = nothing
  r13 = 0.001
  r14 = 0.23
  r15 = 0.25
  r16 = 0.044
  #r17 = nothing
  #r18 = nothing
  #r19 = nothing
  #r20 = nothing
  r21 = 0.034
  r22 = 0.03
  r23 = 0.005
  r24 = 0.001
  #r25 = nothing
  #r26 = nothing
  #r27 = nothing
  #r28 = nothing
  #r29 = nothing
  #r30 = nothing
  r31 = 1
  r32 = 3
  #r33 = nothing
  r34 = 0.25
  r35 = 0.0269
  #r36 = nothing
  r37 = 0.005
  r38 = 0.75
  #r39 = nothing
  #r40 = nothing
  #r41 = nothing
  #r42 = nothing
  #r43 = nothing
  #r44 = nothing
  #r45 = nothing
  #r46 = nothing
  #r47 = nothing
  #r48 = nothing
  #r49 = nothing
  #r50 = nothing
  #r51 = nothing
  r52 = 9
  #r53 = nothing
  #r54 = nothing
  #r55 = nothing
  #r56 = nothing
  #r57 = nothing
  #r58 = nothing
  #r59 = nothing
  #r60 = nothing
  #r61 = nothing
  #r62 = nothing
  #r63 = nothing
  #r64 = nothing
  r65 = 0
  #r66 = nothing
  r67 = 0
  r68 = 0
  #r69 = nothing

  # Inhibition (no values provided in paper. These values have been reverse-engineered)
  Ki_s1p = 0.0356898872566651
  Ki_c1p = 0.178449436283326

  # Initialize fluxes
  f = Vector{Float64}(undef, 69)

  f[1]  = (Vm1) / ( (1 + u[17] / Ki_s1p) * (1 + u[1] / Ki_c1p) )
  f[2]  = (Vm2 * u[3])   / (Km2  + u[3])
  f[3]  = (Vm3 * u[31])  / (Km3  + u[31])
  f[4]  = (Vm4 * u[31])  / (Km4  + u[31])
  f[5]  = (Vm5 * u[18])  / (Km5  + u[18])
  f[6]  = (Vm6 * u[18])  / (Km6  + u[18])
  f[7]  = k7 * u[3]    - r7 * u[9]
  f[8]  = (Vm8 * u[9])   / (Km8  + u[9])
  f[9]  = (Vm9 * u[28])  / (Km9  + u[28])
  f[10] = (Vm10 * u[9])  / (Km10 + u[9])
  f[11] = (Vm11 * u[36]) / (Km11 + u[36])
  f[12] = (Vm12 * u[22]) / (Km12 + u[22])
  f[13] = k13 * u[28]  - r13 * u[24]
  f[14] = k14 * u[36]  - r14 * u[30]
  f[15] = k15 * u[22]  - r15 * u[17]
  f[16] = k16 * u[17]  - r16 * u[21]
  f[17] = (Vm17 * u[21]) / (Km17 + u[21])
  f[18] = (Vm18 * u[35]) / (Km18 + u[35])
  f[19] = (Vm19 * u[35]) / (Km19 + u[35])
  f[20] = (Vm20 * u[8])  / (Km20 + u[8])
  f[21] = k21 * u[35]  - r21 * u[30]
  f[22] = k22 * u[30]  - r22 * u[31]
  f[23] = k23 * u[24]  - r23 * u[25]
  f[24] = k24 * u[24]  - r24 * u[29]
  f[25] = (Vm25 * u[29]) / (Km25 + u[29])
  f[26] = (Vm26 * u[10]) / (Km26 + u[10])
  f[27] = (Vm27 * u[37]) / (Km27 + u[37])
  f[28] = (Vm28 * u[25]) / (Km28 + u[25])
  f[29] = (Vm29 * u[4])  / (Km29 + u[4])
  f[30] = (Vm30 * u[32]) / (Km30 + u[32])
  f[31] = k31 * u[4]  - r31 * u[10]
  f[32] = k32 * u[32] - r32 * u[37]
  f[33] = k33 * u[19]
  f[34] = k34 * u[19] - r34 * u[17]
  f[35] = k35 * u[32] - r35 * u[30]
  f[36] = (Vm36 * u[30]) / (Km36 + u[30])
  f[37] = k37 * u[17] - r37 * u[18]
  f[38] = k38 * u[3]  - r38 * u[8]
  f[39] = k39 * u[3]
  f[40] = k40 * u[3]
  f[41] = (Vm41 * u[7])  / (Km41 + u[7])
  f[42] = k42 * u[12]
  f[43] = (Vm43 * u[11]) / (Km43 + u[11])
  f[44] = (Vm44 * u[16]) / (Km44 + u[16])
  f[45] = k45 * u[14]
  f[46] = k46 * u[15]
  f[47] = (Vm47 * u[13]) / (Km47 + u[13])
  f[48] = k48 * u[29]
  f[49] = (Vm49 * u[26]) / ((Km49 + u[26]) * (1 + u[1] / Ki_c1p) * (1 + u[17] / Ki_s1p) )
  f[50] = (Vm50 * u[5])  / (Km50 + u[5])
  f[51] = k51 * u[33]
  f[52] = k52 * u[30] - r52 * u[34]
  f[53] = (Vm53 * u[34]) / (Km53 + u[34])
  f[54] = (Vm54 * u[20]) / (Km54 + u[20])
  f[55] = (Vm55 * u[6])  / (Km55 + u[6])
  f[56] = k56 * u[27]
  f[57] = (Vm57 * u[27]) / (Km57 + u[27])
  f[58] = (Vm58 * u[6])  / (Km58 + u[6])
  f[59] = (Vm59 * u[6])  / (Km59 + u[6])
  f[60] = (Vm60 * u[1])  / (Km60 + u[1])
  f[61] = k61 * u[1]
  f[62] = (Vm62 * u[2])  / (Km62 + u[2])
  f[63] = (Vm63 * u[10]) / (Km63 + u[10])
  f[64] = k64
  f[65] = k65 - r65 * u[2]
  f[66] = k66
  f[67] = k67 - r67 * u[37]
  f[68] = k68 - r68 * u[23]
  f[69] = (Vm69 * u[23]) / (Km69 + u[23])

  return f
end

# Fluxes for disease condition
function model_ad_fluxes(u::Vector{Float64})

  #Km1 = nothing
  Km2 = 0.081
  Km3 = 0.171
  Km4 = 0.0034
  Km5 = 0.0385
  Km6 = 0.035
  #Km7 = nothing
  Km8 = 0.155
  Km9 = 0.126
  Km10 = 0.0601
  Km11 = 0.0034
  Km12 = 0.025
  #Km13 = nothing
  #Km14 = nothing
  #Km15 = nothing
  #Km16 = nothing
  Km17 = 0.0385
  Km18 = 0.0034
  Km19 = 0.0025
  Km20 = 0.149
  #Km21 = nothing
  #Km22 = nothing
  #Km23 = nothing
  #Km24 = nothing
  Km25 = 0.0455
  Km26 = 0.0601
  Km27 = 0.034
  Km28 = 0.148
  Km29 = 0.0601
  Km30 = 0.00506
  #Km31 = nothing
  #Km32 = nothing
  #Km33 = nothing
  #Km34 = nothing
  #Km35 = nothing
  Km36 = 0.0034
  #Km37 = nothing
  #Km38 = nothing
  #Km39 = nothing
  #Km40 = nothing
  Km41 = 0.04
  #Km42 = nothing
  Km43 = 0.003
  Km44 = 0.003
  #Km45 = nothing
  #Km46 = nothing
  Km47 = 0.019
  #Km48 = nothing
  Km49 = 0.0455
  Km50 = 0.149
  #Km51 = nothing
  #Km52 = nothing
  Km53 = 0.00506
  Km54 = 0.036
  Km55 = 0.02
  #Km56 = nothing
  Km57 = 0.148
  Km58 = 0.081
  Km59 = 0.107
  Km60 = 0.036
  #Km61 = nothing
  Km62 = 0.036
  Km63 = 0.107
  #Km64 = nothing
  #Km65 = nothing
  #Km66 = nothing
  #Km67 = nothing
  #Km68 = nothing
  Km69 = 0.036

  #Vm1 = 0.004
  Vm1 = 0
  Vm2 = 0.32
  Vm3 = 2.4
  Vm4 = 0.875
  Vm5 = 3.1
  Vm6 = 0
  #Vm7 = nothing
  Vm8 = 0.109
  Vm9 = 0.1
  Vm10 = 0.00453
  Vm11 = 0.05
  Vm12 = 1.24
  #Vm13 = nothing
  #Vm14 = nothing
  #Vm15 = nothing
  #Vm16 = nothing
  Vm17 = 3
  Vm18 = 0.07
  Vm19 = 0.1
  Vm20 = 2
  #Vm21 = nothing
  #Vm22 = nothing
  #Vm23 = nothing
  #Vm24 = nothing
  Vm25 = 0.002
  Vm26 = 0.08
  Vm27 = 0.055
  Vm28 = 0.0025
  Vm29 = 0.02
  Vm30 = 0.0015
  #Vm31 = nothing
  #Vm32 = nothing
  #Vm33 = nothing
  #Vm34 = nothing
  #Vm35 = nothing
  Vm36 = 0.001
  #Vm37 = nothing
  #Vm38 = nothing
  #Vm39 = nothing
  #Vm40 = nothing
  Vm41 = 0.01
  #Vm42 = nothing
  Vm43 = 0.001
  Vm44 = 0.001
  #Vm45 = nothing
  #Vm46 = nothing
  Vm47 = 0.0061
  #Vm48 = nothing
  Vm49 = 0.01
  Vm50 = 0.00669
  #Vm51 = nothing
  #Vm52 = nothing
  Vm53 = 0.08
  Vm54 = 2
  Vm55 = 0.3
  #Vm56 = nothing
  Vm57 = 0.07
  Vm58 = 0.529
  Vm59 = 5
  Vm60 = 0.75
  #Vm61 = nothing
  Vm62 = 0.78
  Vm63 = 0.5
  #Vm64 = nothing
  #Vm65 = nothing
  #Vm66 = nothing
  #Vm67 = nothing
  #Vm68 = nothing
  Vm69 = 0.43

  #k1 = nothing
  #k2 = nothing
  #k3 = nothing
  #k4 = nothing
  #k5 = nothing
  #k6 = nothing
  k7 = 0.8
  #k8 = nothing
  #k9 = nothing
  #k10 = nothing
  #k11 = nothing
  #k12 = nothing
  k13 = 0.12
  k14 = 0.5
  k15 = 0.44
  k16 = 0.25
  #k17 = nothing
  #k18 = nothing
  #k19 = nothing
  #k20 = nothing
  k21 = 0.43
  k22 = 23
  k23 = 0.01
  k24 = 0.015
  #k25 = nothing
  #k26 = nothing
  #k27 = nothing
  #k28 = nothing
  #k29 = nothing
  #k30 = nothing
  k31 = 1
  k32 = 1
  k33 = 1
  k34 = 0.4
  k35 = 3
  #k36 = nothing
  k37 = 4.5
  k38 = 1
  k39 = 1
  k40 = 0.02
  #k41 = nothing
  k42 = 1
  #k43 = nothing
  #k44 = nothing
  k45 = 0.2
  k46 = 0.03
  #k47 = nothing
  k48 = 0.011
  #k49 = nothing
  #k50 = nothing
  k51 = 1
  k52 = 0.2
  #k53 = nothing
  #k54 = nothing
  #k55 = nothing
  k56 = 0.045
  #k57 = nothing
  #k58 = nothing
  #k59 = nothing
  #k60 = nothing
  k61 = 2.5
  #k62 = nothing
  #k63 = nothing
  k64 = 0
  k65 = 0
  k66 = 0
  k67 = 0
  k68 = 0
  #k69 = nothing

  #r1 = nothing
  #r2 = nothing
  #r3 = nothing
  #r4 = nothing
  #r5 = nothing
  #r6 = nothing
  r7 = 0.1
  #r8 = nothing
  #r9 = nothing
  #r10 = nothing
  #r11 = nothing
  #r12 = nothing
  r13 = 0.001
  r14 = 0.23
  r15 = 0.25
  r16 = 0.044
  #r17 = nothing
  #r18 = nothing
  #r19 = nothing
  #r20 = nothing
  r21 = 0.034
  r22 = 0.03
  r23 = 0.005
  r24 = 0.001
  #r25 = nothing
  #r26 = nothing
  #r27 = nothing
  #r28 = nothing
  #r29 = nothing
  #r30 = nothing
  r31 = 1
  r32 = 3
  #r33 = nothing
  r34 = 0.25
  r35 = 0.0269
  #r36 = nothing
  r37 = 0.005
  r38 = 0.75
  #r39 = nothing
  #r40 = nothing
  #r41 = nothing
  #r42 = nothing
  #r43 = nothing
  #r44 = nothing
  #r45 = nothing
  #r46 = nothing
  #r47 = nothing
  #r48 = nothing
  #r49 = nothing
  #r50 = nothing
  #r51 = nothing
  r52 = 9
  #r53 = nothing
  #r54 = nothing
  #r55 = nothing
  #r56 = nothing
  #r57 = nothing
  #r58 = nothing
  #r59 = nothing
  #r60 = nothing
  #r61 = nothing
  #r62 = nothing
  #r63 = nothing
  #r64 = nothing
  r65 = 0
  #r66 = nothing
  r67 = 0
  r68 = 0
  #r69 = nothing

  # Inhibition (no values provided in paper. These values have been reverse-engineered)
  Ki_s1p = 0.0356898872566651
  Ki_c1p = 0.178449436283326

  # Initialize fluxes
  f = Vector{Float64}(undef, 69)

  f[1]  = (Vm1) / ( (1 + u[17] / Ki_s1p) * (1 + u[1] / Ki_c1p) )
  f[2]  = (Vm2 * u[3])   / (Km2  + u[3])
  f[3]  = (Vm3 * u[31])  / (Km3  + u[31])
  f[4]  = (Vm4 * u[31])  / (Km4  + u[31])
  f[5]  = (Vm5 * u[18])  / (Km5  + u[18])
  f[6]  = (Vm6 * u[18])  / (Km6  + u[18])
  f[7]  = k7 * u[3]    - r7 * u[9]
  f[8]  = (Vm8 * u[9])   / (Km8  + u[9])
  f[9]  = (Vm9 * u[28])  / (Km9  + u[28])
  f[10] = (Vm10 * u[9])  / (Km10 + u[9])
  f[11] = (Vm11 * u[36]) / (Km11 + u[36])
  f[12] = (Vm12 * u[22]) / (Km12 + u[22])
  f[13] = k13 * u[28]  - r13 * u[24]
  f[14] = k14 * u[36]  - r14 * u[30]
  f[15] = k15 * u[22]  - r15 * u[17]
  f[16] = k16 * u[17]  - r16 * u[21]
  f[17] = (Vm17 * u[21]) / (Km17 + u[21])
  f[18] = (Vm18 * u[35]) / (Km18 + u[35])
  f[19] = (Vm19 * u[35]) / (Km19 + u[35])
  f[20] = (Vm20 * u[8])  / (Km20 + u[8])
  f[21] = k21 * u[35]  - r21 * u[30]
  f[22] = k22 * u[30]  - r22 * u[31]
  f[23] = k23 * u[24]  - r23 * u[25]
  f[24] = k24 * u[24]  - r24 * u[29]
  f[25] = (Vm25 * u[29]) / (Km25 + u[29])
  f[26] = (Vm26 * u[10]) / (Km26 + u[10])
  f[27] = (Vm27 * u[37]) / (Km27 + u[37])
  f[28] = (Vm28 * u[25]) / (Km28 + u[25])
  f[29] = (Vm29 * u[4])  / (Km29 + u[4])
  f[30] = (Vm30 * u[32]) / (Km30 + u[32])
  f[31] = k31 * u[4]  - r31 * u[10]
  f[32] = k32 * u[32] - r32 * u[37]
  f[33] = k33 * u[19]
  f[34] = k34 * u[19] - r34 * u[17]
  f[35] = k35 * u[32] - r35 * u[30]
  f[36] = (Vm36 * u[30]) / (Km36 + u[30])
  f[37] = k37 * u[17] - r37 * u[18]
  f[38] = k38 * u[3]  - r38 * u[8]
  f[39] = k39 * u[3]
  f[40] = k40 * u[3]
  f[41] = (Vm41 * u[7])  / (Km41 + u[7])
  f[42] = k42 * u[12]
  f[43] = (Vm43 * u[11]) / (Km43 + u[11])
  f[44] = (Vm44 * u[16]) / (Km44 + u[16])
  f[45] = k45 * u[14]
  f[46] = k46 * u[15]
  f[47] = (Vm47 * u[13]) / (Km47 + u[13])
  f[48] = k48 * u[29]
  f[49] = (Vm49 * u[26]) / ((Km49 + u[26]) * (1 + u[1] / Ki_c1p) * (1 + u[17] / Ki_s1p) )
  f[50] = (Vm50 * u[5])  / (Km50 + u[5])
  f[51] = k51 * u[33]
  f[52] = k52 * u[30] - r52 * u[34]
  f[53] = (Vm53 * u[34]) / (Km53 + u[34])
  f[54] = (Vm54 * u[20]) / (Km54 + u[20])
  f[55] = (Vm55 * u[6])  / (Km55 + u[6])
  f[56] = k56 * u[27]
  f[57] = (Vm57 * u[27]) / (Km57 + u[27])
  f[58] = (Vm58 * u[6])  / (Km58 + u[6])
  f[59] = (Vm59 * u[6])  / (Km59 + u[6])
  f[60] = (Vm60 * u[1])  / (Km60 + u[1])
  f[61] = k61 * u[1]
  f[62] = (Vm62 * u[2])  / (Km62 + u[2])
  f[63] = (Vm63 * u[10]) / (Km63 + u[10])
  f[64] = k64
  f[65] = k65 - r65 * u[2]
  f[66] = k66
  f[67] = k67 - r67 * u[37]
  f[68] = k68 - r68 * u[23]
  f[69] = (Vm69 * u[23]) / (Km69 + u[23])

  return f
end

