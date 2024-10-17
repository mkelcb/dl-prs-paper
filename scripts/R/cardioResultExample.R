# cardio metabolic traits, linear vs nonlinear performance
# linear    nonlinear
# Ischemic stroke  EFO_0000712                    0.362   0.421
# Fasting glucose  EFO_0004465                    0.349   0.203
# Venous thromboembolism    EFO_0004286           0.595   0.945
# Glycated haemoglobin levels (HbA1c) EFO_0004541 0.954   0.983
# type 2 diabetes MONDO_0005148                   0.484   0.586
# 


linear = c(0.362,0.349,0.595,0.954,0.484)

nonlinear =c(0.421,0.203, 0.945,0.983,0.586)


median(linear) # 0.484
mean(linear) # 0.5488
median(nonlinear) # 0.586
mean(nonlinear) # 0.6276
# how much worse is the NN than the PGS catalog baseline?
1-mean(nonlinear) # 0.3724

# what is the average difference
mean(nonlinear) - mean(linear) #  0.0788
