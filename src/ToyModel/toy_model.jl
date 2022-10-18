function load_toy_model(; dobox = false)
    S = 
    [
        #|  us   |  ua   |  ul   |  vg   |  va   |  vr   |  vf   |   z   | 
            1.0     0.0     0.0    -1.0     0.0     0.0     0.0     0.0     # |  S
            0.0     0.0     0.0     1.0     1.0    -1.0    -1.0     0.0     # |  P
            0.0     0.0     0.0     1.0     0.0     5.0     0.0    -1.0     # |  E
            0.0     1.0     0.0     0.0    -1.0     0.0     0.0     0.0     # |  A
            0.0     0.0     1.0     0.0     0.0     0.0     1.0     0.0     # |  W
    ]
    
    rxns = [  "Ex_s" , "Ex_a" , "Ex_l"  , "vspe" , "va"   ,  "vpe" ,  "vpl" ,  "biomass"   ]
    c    = [  0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0   ,  1.0   ]
    lb   = [  0.0   ,  0.0   , -100.0 ,  0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0   ]
    ub   = [  10.0  ,  100.0 ,  0.0   ,  100.0 ,  100.0 ,  100.0 ,  100.0 ,  1.0 ]
    mets = [  "S"   ,   "P"  ,   "E"  ,   "A"  , "W"   ]
    b    = [  0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0  ]


    model = MetNets.MetNet(;S, b, c, lb, ub, rxns, mets)
    return dobox ? MetLP.fva_preprocess(model; verbose = false) : model
end