function load_toy_model3(; dobox = false)
    S = 
    [
        #|  Ex_s | Ex_a  |  Ex_l  | vspe  | vap  |  vas |  vpe  |  vpl  | biomass    | 
           -1.0     0.0     0.0    -1.0     0.0    1.0     0.0     0.0    -1.0     # |  S
            0.0     0.0     0.0     1.0     1.0    0.0    -1.0    -1.0     0.0     # |  P
            0.0     0.0     0.0     1.0     0.0    0.0     1.0     0.0    -1.0     # |  E
            0.0    -1.0     0.0     0.0    -1.0   -1.0     0.0     0.0    -1.0     # |  A
            0.0     0.0    -1.0     0.0     0.0    0.0     0.0     1.0     0.0     # |  L
    ]
    rxns = [ "Ex_s" , "Ex_a" , "Ex_l" , "vspe" , "vap"  , "vas"  ,  "vpe" ,  "vpl" ,  "biomass" ]
    c    = [  0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0   ,  1.0       ]
    lb   = [ -1.0   , -100.0 ,  0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0       ]
    ub   = [  0.0   ,   0.0  , 100.0  ,  100.0 ,  100.0 ,  100.0 ,  100.0 ,  100.0 ,  10.0      ]
    mets = [  "S"   ,   "P"  ,   "E"  ,  "A"   ,  "L"  ]
    b    = [  0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0  ]

    model = MetNets.MetNet(;S, b, c, lb, ub, rxns, mets)
    return dobox ? MetLP.fva_preprocess(model; verbose = false) : model
end