#   IAU FK4 - FK5 Conversion Constants

#   FK4 B1950.0 to FK5 J2000.0
#   A, Adot pv-vector (cf. Seidelmann eq. 3.591-2)

const A_fk4_fk5::Vector{Float64} = [-1.62557e-6, -0.31919e-6, -0.13843e-6, +1.245e-3,   -1.580e-3,   -0.659e-3]

#   M matrix pv-vector (cf. Seidelmann eq. 3.591-4)
const M_fk4_fk5::Matrix{Float64} = [
    +0.9999256782 -0.0111820611 -0.0048579477 +0.00000242395018 -0.00000002710663 -0.00000001177656;
    +0.0111820610 +0.9999374784 -0.0000271765 +0.00000002710663 +0.00000242397878 -0.00000000006587;
    +0.0048579479 -0.0000271474 +0.9999881997 +0.00000001177656 -0.00000000006582 +0.00000242410173;
    -0.000551     -0.238565     +0.435739     +0.99994704       -0.01118251       -0.00485767      ;
    +0.238514     -0.002667     -0.008541     +0.01118251       +0.99995883       -0.00002718      ;
    -0.435623     +0.012254     +0.002117     +0.00485767       -0.00002714       +1.00000956      ]


const Mfk4fk5 = [[[[+0.9999256782,     -0.0111820611,     -0.0048579477    ],
                   [+0.00000242395018, -0.00000002710663, -0.00000001177656]],
                  [[+0.0111820610,     +0.9999374784,     -0.0000271765    ],
                   [+0.00000002710663, +0.00000242397878, -0.00000000006587]],
                  [[+0.0048579479,     -0.0000271474,     +0.9999881997,   ],
                   [+0.00000001177656, -0.00000000006582, +0.00000242410173]]],
                 [[[-0.000551,         -0.238565,         +0.435739        ],
                   [+0.99994704,       -0.01118251,       -0.00485767      ]],
                  [[+0.238514,         -0.002667,         -0.008541        ],
                   [+0.01118251,       +0.99995883,       -0.00002718      ]],
                  [[-0.435623,         +0.012254,         +0.002117        ],
                   [+0.00485767,       -0.00002714,       +1.00000956      ]]]]

#   FK5 J2000.0 to FK4 J2000.0

const M_fk5_fk4::Matrix{Float64} = [
    +0.9999256795 +0.0111814828 +0.0048590039 -0.00000242389840 -0.00000002710544 -0.00000001177742;
    -0.0111814828 +0.9999374849 -0.0000271771 +0.00000002710544 -0.00000242392702 +0.00000000006585;
    -0.0048590040 -0.0000271557 +0.9999881946 +0.00000001177742 +0.00000000006585 -0.00000242404995;
    -0.000551     +0.238509     -0.435614     +0.99990432       +0.01118145       +0.00485852      ;
    -0.238560     -0.002667     +0.012254     -0.01118145       +0.99991613       -0.00002717      ;
    +0.435730     -0.008541     +0.002117     -0.00485852       -0.00002716       +0.99996684      ]

#  3x2 matrix of pv-vectors (cf. Seidelmann 3.592-1, inverse matrix M)
const Mfk5fk4 = [[[[+0.9999256795,     +0.0111814828,     +0.0048590039,   ],
                   [-0.00000242389840, -0.00000002710544, -0.00000001177742]],
                  [[-0.0111814828,     +0.9999374849,     -0.0000271771,   ],
                   [+0.00000002710544, -0.00000242392702, +0.00000000006585]],
                  [[-0.0048590040,     -0.0000271557,     +0.9999881946,   ],
                   [+0.00000001177742, +0.00000000006585, -0.00000242404995]]],
                 [[[-0.000551,         +0.238509,         -0.435614,       ],
                   [+0.99990432,       +0.01118145,       +0.00485852      ]],
                  [[-0.238560,         -0.002667,         +0.012254,       ],
                   [-0.01118145,       +0.99991613,       -0.00002717      ]],
                  [[+0.435730,         -0.008541,         +0.002117,       ],
                   [-0.00485852,       -0.00002716,       +0.99996684      ]]]]