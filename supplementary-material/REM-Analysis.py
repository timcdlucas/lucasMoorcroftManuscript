"""
Systematic analysis of REM models
Tim Lucas 
01/10/13
"""


from sympy import *
import numpy as np
import matplotlib.pyplot as pl
from datetime import datetime
import Image as Im


# Use LaTeX printing
from sympy import init_printing ;
init_printing()
# Make LaTeX output white. Because I use a dark theme
init_printing(forecolor="White") 


# Load symbols used for symbolic maths
t, a, r, x2, x3, x4, x1 = symbols('theta alpha r x_2 x_3 x_4 x_1', positive=True)
r1 = {r:1} # useful for lots of checks


# Define functions to neaten up later code.

# Calculate the final profile averaged over pi.
def calcModel(model):
        x = pi**-1 * sum( [integrate(m[0], m[1:]) for m in model] ).simplify().trigsimp()
        return x

# Do the replacements fit within the area defined by the conditions? 
def confirmReplacements(conds, reps):
        if not all([c.subs(reps) for c in eval(conds)]):
                print('reps' + conds[4:] + ' incorrect')

# is average profile in range 0r-2r?
def profileRange(prof, reps):
        if not 0 <= eval(prof).subs(dict(reps, **r1)) <= 2:
                print('Total ' + prof + ' not in 0, 2r')

# Are the individuals integrals >0r
def intsPositive(model, reps):
        m = eval(model)
        for i in range(len(m)):
                if not integrate(m[i][0], m[i][1:]).subs(dict(reps, **r1)) > 0:
                    print('Integral ' + str(i+1) + ' in ' + model + ' is negative')

# Are the individual averaged integrals between 0 and  2r
def intsRange(model, reps):
        m = eval(model)
        for i in range(len(m)):
                if not 0 <= (integrate(m[i][0], m[i][1:])/(m[i][3]-m[i][2])).subs(dict(reps, **r1)) <= 2:
                        print('Integral ' + str(i+1) + ' in ' + model + ' has averaged integral outside 0<p<2r')

# Are the bounds the correct way around
def checkBounds(model, reps):
        m = eval(model)
        for i in range(len(m)):
                if not (m[i][3]-m[i][2]).subs(reps) > 0:
                        print('Bounds ' + str(i+1) + ' in ' + model + ' has lower bounds bigger than upper bounds')        

# create latex strings with the 1) the integral equation that defines it and 2)  the final calculated model.
# There's some if statements to split longer equations on two lines and get +s in the right place.
def parseLaTeX(prof):
        m = eval( 'm' + prof[1:] )
        f = open('/home/tim/Dropbox/phd/Analysis/REM-chapter/latexFiles/'+prof+'.tex', 'w')
        f.write('\\begin{align}\n    ' + prof + ' =&\\frac{1}{\pi} \left(\;\;')
        for i in range(len(m)):
                f.write('\int\limits_{'+latex(m[i][2], order='rev-lex')+'}^{'+latex(m[i][3], order='rev-lex')+'}'+latex(m[i][0], order='rev-lex')+'\;\mathrm{d}' +latex(m[i][1]))
                if len(m)>3 and i==(len(m)/2)-1:
                        f.write( '\\right.\\notag\\\\\n &\left.' )
                if i<len(m)-1:
                        f.write('+')                                            
        f.write('\\right)\label{' + prof + 'Def}\\\\\n    ')
        f.write(prof + ' =& ' + latex(eval(prof)) + '\label{' + prof + 'Sln}\n\\end{align}')
        f.close()


# Apply all checks.
def allChecks(prof):
        model = 'm' + prof[1:]
        reps = eval('rep' + prof[1:])
        conds = 'cond' + prof[1:]
        confirmReplacements(conds, reps)
        profileRange(prof, reps)
        intsPositive(model, reps)
        intsRange(model, reps)
        checkBounds(model, reps)

#########################################################
# 221 animal: a = 2*pi.  sensor: t > pi, a > 3pi - t  #
#########################################################



m221 = [ [2*r,                 x1, pi/2, t/2        ],
         [r + r*cos(x1 - t/2), x1, t/2,  pi         ],
         [r + r*cos(x1 + t/2), x1, pi,   2*pi-t/2   ],
         [2*r,                 x1, 2*pi-t/2, 3*pi/2 ] ]

# Replacement values in range
rep221 = {t:3*pi/2, a:2*pi} 

# Define conditions for model
cond221 = [pi <= t, a >= 3*pi - t]

# Calculate model, run checks, write output.
p221 = calcModel(m221)
allChecks('p221')
parseLaTeX('p221')

###############################################################################
# 222 animal: a > pi.  sensor: t > pi Condition: a < 3pi - t, a > 4pi - 2t  #
###############################################################################



m222 = [ [2*r,                 x1, pi/2, t/2        ],
         [r + r*cos(x1 - t/2), x1, t/2,  5*pi/2 - t/2 - a/2 ],
         [r + r*cos(x1 + t/2), x1, 5*pi/2 - t/2 - a/2,   2*pi-t/2 ],
         [2*r,                 x1, 2*pi-t/2, 3*pi/2 ] ]


# Replacement values in range
rep222 = {t:5*pi/3, a:4*pi/3-0.1} 

# Define conditions for model
cond222 = [pi <= t, a >= pi, a <= 3*pi - t, a >= 4*pi - 2*t]

# Calculate model, run checks, write output.
p222 = calcModel(m222)
allChecks('p222')
parseLaTeX('p222')


##################################################################
# 223 animal: a > pi.  sensor: t > pi Condition: a < 4pi - 2t  #
##################################################################



m223 = [ [2*r,                 x1, pi/2, t/2        ],
         [r + r*cos(x1 - t/2), x1, t/2,  t/2 + pi/2         ],
         [r                  , x1, t/2 + pi/2,   5*pi/2 - t/2 - a/2 ],
         [r + r*cos(x1 + t/2), x1, 5*pi/2 - t/2 - a/2,   2*pi-t/2 ],
         [2*r,                 x1, 2*pi-t/2, 3*pi/2 ] ]


# Replacement values in range
rep223 = {t:5*pi/4-0.1, a:3*pi/2}

# Define conditions for model
cond223 = [pi <= t, a >= pi, a <= 4*pi - 2*t]

# Calculate model, run checks, write output.
p223 = calcModel(m223)
allChecks('p223')
parseLaTeX('p223')




########################################################
# 131 animal: a = 2*pi.   sensor:  pi/2 <= t <= pi      #
########################################################

m131 = [ [2*r*sin(t/2)*sin(x2), x2, t/2,      pi/2     ],
        [r - r*cos(x4 - t),     x4, 0,        t - pi/2 ],
        [r,                     x4, t - pi/2, pi/2     ],
        [r - r*cos(x4),         x4, pi/2,     t        ],
        [2*r*sin(t/2)*sin(x2),  x2, t/2,      pi/2     ] ]

# Replacement values in range
rep131 = {t:3*pi/4} 

# Define conditions for model
cond131 = [pi/2 <= t, t <= pi]

# Calculate model, run checks, write output.
p131 = calcModel(m131)
allChecks('p131')
parseLaTeX('p131')



#################################################################################
# 231 animal: a > pi.  Sensor: pi/2 <= t <= pi. Condition: a > 2pi - t          #
#################################################################################


m231 = [ [2*r*sin(t/2)*sin(x2), x2, t/2,          pi/2        ],
         [r - r*cos(x4 - t),    x4, 0,            t - pi/2    ],
         [r,                    x4, t - pi/2,     3*pi/2 - a/2],
         [r - r*cos(x4),        x4, 3*pi/2 - a/2, t           ],
         [2*r*sin(t/2)*sin(x2), x2, t/2,          pi/2        ] ]


rep231 = {t:3*pi/4, a:15*pi/8} # Replacement values in range

# Define conditions for model
cond231 = [a > pi, pi/2 <= t, t <= pi, a >= 3*pi - 2*t]

# Calculate model, run checks, write output.
p231 = calcModel(m231)
allChecks('p231')
parseLaTeX('p231')


#################################################################################
# 232 animal: a > pi.  Sensor: pi/2 <= t <= pi. Cond: 2pi - t < a < 3pi - 2t    #
#################################################################################


m232 = [ [2*r*sin(t/2)*sin(x2), x2, t/2,                pi/2              ],
         [r - r*cos(x4 - t),    x4, 0,                  t - pi/2          ],
         [r,                    x4, t - pi/2,           t                 ],
         [r*cos(x2 - t/2),      x2, t/2,                3*pi/2 - a/2 - t/2],
         [2*r*sin(t/2)*sin(x2), x2, 3*pi/2 - a/2 - t/2, pi/2              ] ]


rep232 = {t:5*pi/8, a:6*pi/4} # Replacement values in range

# Define conditions for model
cond232 = [a > pi, pi/2 <= t, t <= pi, 2*pi - t <= a, a <= 3*pi - 2*t]

# Calculate model, run checks, write output.
p232 = calcModel(m232)
allChecks('p232')
parseLaTeX('p232')


#################################################################################
# 233 animal:  a > pi.  Sensor: pi/2 <= t <= pi. Condition: a <= 2pi - t      #
#################################################################################

m233 = [ [2*r*sin(t/2)*sin(x2), x2, t/2, pi/2],
         [r - r*cos(x4 - t),    x4, 0, t - pi/2],
         [r,                    x4, t - pi/2, t],
         [r*cos(x2 - t/2),      x2, t/2, a/2 + t/2 - pi/2] ]

rep233 = {t:3*pi/4, a:9*pi/8} # Replacement values in range

# Define conditions for model
cond233 = [a > pi,  pi/2 <= t, t <= pi, a <= 2*pi - t]

# Calculate model, run checks, write output.
p233 = calcModel(m233)
allChecks('p233')
parseLaTeX('p233')

###############################################################################
# 141 animal: a=2pi.  Sensor: t <= pi/2.                                      #
###############################################################################

m141 = [ [2*r*sin(t/2)*sin(x2), x2, pi/2 - t/2, pi/2],
         [r*sin(x3),            x3, t,          pi/2],
         [r,                    x4, 0*t,          t],
         [r*sin(x3),            x3, t,          pi/2],
         [2*r*sin(t/2)*sin(x2), x2, pi/2 - t/2, pi/2] ]


rep141 = {t:3*pi/8, a:2*pi} # Replacement values in range

# Define conditions for model
cond141 = [ t <= pi/2  ]

# Calculate model, run checks, write output.
p141 = calcModel(m141)
allChecks('p141')
parseLaTeX('p141')


###############################################################################
# 241 animal: a>pi.  Sensor: t <= pi/2. Condition: 2*pi - t < a               #
###############################################################################

m241 = [ [2*r*sin(t/2)*sin(x2), x2, pi/2 - t/2, pi/2],
         [r*sin(x3),            x3, t,          pi/2],
         [r,                    x4, 0,          t],
         [r*sin(x3),            x3, t,          pi/2],
         [r*cos(x2 - t/2),      x2, pi/2 - t/2, 3*pi/2 - t/2 - a/2],
         [2*r*sin(t/2)*sin(x2), x2, 3*pi/2 - t/2 - a/2, pi/2] ]


rep241 = {t:3*pi/8, a:29*pi/16} # Replacement values in range

# Define conditions for model
cond241 = [a >= pi, t <= pi/2, 2*pi - t <= a  ]

# Calculate model, run checks, write output.
p241 = calcModel(m241)
allChecks('p241')
parseLaTeX('p241')

####################################################################################
# 242 animal: a>pi.  Sensor: t <= pi/2. Condition:  2*pi - 2*t <= a <= 2*pi - t  #
####################################################################################

m242 = [ [2*r*sin(t/2)*sin(x2), x2, pi/2 - t/2, pi/2],
         [r*sin(x3),            x3, t,          pi/2],
         [r,                    x4, 0,          t],
         [r*sin(x3),            x3, t,          pi/2],
         [r*cos(x2 - t/2),      x2, pi/2 - t/2, a/2 + t/2 - pi/2] ]

rep242 = {t:3*pi/8, a:3*pi/2} # Replacement values in range

# Define conditions for model
cond242 = [a >= pi, t <= pi/2, 2*pi - 2*t <= a, a <= 2*pi - t]

# Calculate model, run checks, write output.
p242 = calcModel(m242)
allChecks('p242')
parseLaTeX('p242')


#####################################################################
# 243 animal: a>pi.  Sensor: t <= pi/2. Condition: a <= 2pi - 2t  #
#####################################################################

m243 = [ [2*r*sin(t/2)*sin(x2), x2, pi/2 - t/2, pi/2],
         [r*sin(x3),            x3, t,          pi/2],
         [r,                    x4, 0,          t   ],
         [r*sin(x3),            x3, pi - a/2,   pi/2] ]


rep243 = {t:pi/9, a:10*pi/9} # Replacement values in range

# Define conditions for model
cond243 = [t <= pi/2, a >= pi, a <= 2*pi - 2*t]

# Calculate model, run checks, write output.
p243 = calcModel(m243)
allChecks('p243')
parseLaTeX('p243')

####################################################
# 311 animal: a <= pi.  Sensor: t =2pi.            #
####################################################


m311 = [ [ 2*r*sin(a/2),                        x1, pi/2, 3*pi/2       ],
         ]


rep311 = {a:pi/4} # Replacement values in range

# Define conditions for model
cond311 = [a <= pi]

# Calculate model, run checks, write output.
p311 = calcModel(m311)
allChecks('p311')
parseLaTeX('p311')



#################################################################################
# 321 animal: a <= pi.  Sensor: t > pi. Condition: a > 2pi - t, a > 4pi - 2t    #
#################################################################################


m321 = [ [ 2*r*sin(a/2),                        x1, pi/2,               t/2 + pi/2 - a/2       ],
         [ r*sin(a/2) + r*cos(x1 - t/2),        x1, t/2 + pi/2 - a/2,   5*pi/2 - a/2 - t/2 ], 
         [ 2*r*sin(a/2),                        x1, 5*pi/2 - a/2 - t/2, 3*pi/2]  ]


rep321 = {t:19*pi/10, a:pi/2} # Replacement values in range

# Define conditions for model
cond321 = [a <= pi, t >= pi,  a >= 4*pi - 2*t]

# Calculate model, run checks, write output.
p321 = calcModel(m321)
allChecks('p321')
parseLaTeX('p321')



#########################################################################################
# 322 animal: a <= pi.  Sensor: t > pi. Condition: 2pi - t < a < 4pi - 2t #
#########################################################################################

m322 = [ [ 2*r*sin(a/2),                        x1, pi/2,               t/2 + pi/2 - a/2  ],
         [ r*sin(a/2) + r*cos(x1 - t/2),        x1, t/2 + pi/2 - a/2,   t/2 + pi/2        ],
         [ r*sin(a/2),                          x1, t/2 + pi/2,         5*pi/2 - a/2 - t/2],
         [ 2*r*sin(a/2),                        x1, 5*pi/2 - a/2 - t/2, 3*pi/2            ] ]

rep322 = {t:3*pi/2 + 0.1, a:pi/2} # Replacement values in range

# Define conditions for model
cond322 = [a <= pi, t >= pi,  a >= 2*pi - t, a <= 4*pi - 2*t]

# Calculate model, run checks, write output.
p322 = calcModel(m322)
allChecks('p322')
parseLaTeX('p322')



###################################################################################
# 323 animal: a <= pi.  Sensor: t > pi. Condition: a <= 4*pi - 2*t and a < 2*pi - t #
###################################################################################

m323 = [ [ 2*r*sin(a/2),                       x1, pi/2,             t/2 + pi/2 - a/2  ],
         [ r*sin(a/2) + r*cos(x1 - t/2),       x1, t/2 + pi/2 - a/2, t/2 + pi/2        ], 
         [ r*sin(a/2),                         x1, t/2 + pi/2,       t/2 + pi/2 + a/2  ] ]


rep323 = {t:3*pi/2, a:pi/3} # Replacement values in range


# Define conditions for model
cond323 = [a <= pi, t >= pi/2, a <= 4*pi - 2*t , a <= 2*pi - t]

# Calculate model, run checks, write output.
p323 = calcModel(m323)
allChecks('p323')
parseLaTeX('p323')


###############################################################################


"""
Ccomplex profiles for a <= pi/2 
These were specified using a very roundabout way that I realised isn't necessary.
Worth keeping them here just for the record.

# p-l-r for x2 profil. Calculated by AE in fig 22.4 minus AE in fig 22.3
p1 = (2*r*sin(t/4 - x2/2 + pi/4 + a/4)*sin(a/4 + pi/4 + x2/2 - t/4) - \
     2*r*sin((pi - a - 2*x2 + t)/4)*sin((pi - a + 2*x2 - t)/4)).simplify()

# p-l for x2 profiles
p2 = (2*r*sin(t/2)*sin(x2) - 2*r*sin((pi - a - 2*x2 + t)/4)*sin((pi - a + 2*x2 - t)/4)).simplify()

# p-l for x3 profile. 
p3 = (r*sin(x3) - (2*r*sin(x3/2 - a/4)*sin(pi/2 - x3/2 - a/4)).simplify()).trigsimp()
"""


###########################################################################################
# 331 animal: a <= pi.  Sensor: pi/2 <= t <= pi. Condition: a >= t and a/2 >= t - pi/2 #
###########################################################################################


m331 =  [ [2*r*sin(t/2)*sin(x2),              x2, pi/2 - a/2 + t/2, pi/2            ],
          [r*sin(a/2) - r*cos(x2 + t/2),      x2, t/2,              pi/2 - a/2 + t/2],
          [r*sin(a/2) - r*cos(x4 - t),        x4, 0,                t - pi/2        ],
          [r*sin(a/2),                        x4, t-pi/2,           t - pi/2 + a/2  ] ]


rep331 = {t:5*pi/8, a:6*pi/8} # Replacement values in range

# Define conditions for model
cond331 = [a <= pi, pi/2 <= t, t <= pi, a >= t, a/2 >= t - pi/2]

# Calculate model, run checks, write output.
p331 = calcModel(m331)
allChecks('p331')
parseLaTeX('p331')




##########################################################################################
# 332 animal: a <= pi.  Sensor: pi/2 <= t <= pi. Condition: a <= t and a/2 >= t- pi/2 #
##########################################################################################


m332 =  [ [2*r*sin(a/2),                 x2, pi/2 + a/2 - t/2, pi/2             ],
          [r*sin(a/2) - r*cos(x2 + t/2), x2, t/2,              pi/2 + a/2 - t/2],
          [r*sin(a/2) - r*cos(x4 - t),   x4, 0*t,              t - pi/2       ],
          [r*sin(a/2),                   x4, t - pi/2,         t - pi/2 + a/2 ] ]


rep332 = {t:7*pi/8, a:7*pi/8-0.1} # Replacement values in range

# Define conditions for model
cond332 = [a <= pi, pi/2 <= t, t <= pi, a/2 <= t/2, a/2 >= t - pi/2]

# Calculate model, run checks, write output.
p332 = calcModel(m332)
allChecks('p332')
parseLaTeX('p332')







##########################################################################################
# 333 animal: a <= pi.  Sensor: pi/2 <= t <= pi. Condition: a <= t and a/2 <= t- pi/2 #
##########################################################################################



m333 =  [ [2*r*sin(a/2),                      x2, t/2,            pi/2           ],
          [2*r*sin(a/2),                      x4, 0,              t - pi/2 - a/2 ],
          [r*sin(a/2) - r*cos(x4 - t),        x4, t - pi/2 - a/2, t - pi/2       ],
          [r*sin(a/2),                        x4, t - pi/2,       t - pi/2 + a/2 ] ]


rep333 = {t:7*pi/8, a:2*pi/8} # Replacement values in range

# Define conditions for model
cond333 = [a <= pi, pi/2 <= t, t <= pi, a/2 <= t/2, a/2 <= t - pi/2]

# Calculate model, run checks, write output.
p333 = calcModel(m333)
allChecks('p333')
parseLaTeX('p333')





##################################################################################
# 341 animal: a <= pi.  Sensor: t <= pi/2. Condition: a > pi - 2t &  a <= t      #
##################################################################################


m341 = [ [2*r*sin(a/2),                 x2, pi/2 - t/2 + a/2, pi/2            ],
         [r*sin(a/2) - r*cos(x2 + t/2), x2, pi/2 - t/2,       pi/2 - t/2 + a/2],
         [r*sin(a/2),                   x3, t,                pi/2            ],
         [r*sin(a/2),                   x4, 0,                a/2 + t - pi/2  ] ]

rep341 = {t:pi/2-0.1, a:pi/4} # Replacement values in range

# Define conditions for model
cond341 = [a <= pi,  t <= pi/2,  a >= pi - 2*t,  a <= t]

# Calculate model, run checks, write output.
p341 = calcModel(m341)
allChecks('p341')
parseLaTeX('p341')


######################################################################################
# 342 animal: a <= pi.  Sensor: t <= pi/2. Condition: a > pi - 2t &  t <= a <= 2t    #
######################################################################################


m342 = [ [2*r*sin(t/2)*sin(x2),         x2, pi/2 + t/2 - a/2, pi/2            ],
         [r*sin(a/2) - r*cos(x2 + t/2), x2, pi/2 - t/2,       pi/2 + t/2 - a/2],
         [r*sin(a/2),                   x3, t,                pi/2        ],
         [r*sin(a/2),                   x4, 0,                a/2 + t -pi/2   ] ]


rep342 = {t:pi/2-0.1, a:pi/2} # Replacement values in range

# define conditions for model
cond342 = [a <= pi,  t <= pi/2,  a >= pi - 2*t,  t <= a, a <= 2*t]


# Calculate model, run checks, write output.
p342 = calcModel(m342)
allChecks('p342')
parseLaTeX('p342')





##################################################################################
# 343 animal: a <= pi.  Sensor: t <= pi/2. Condition: a > pi - 2t &  a > 2t      #
##################################################################################



m343 = [ [2*r*sin(t/2)*sin(x2), x2, pi/2 - t/2, pi/2            ],
         [r*sin(x3),            x3, t,          a/2             ],
         [r*sin(a/2),           x3, a/2,        pi/2            ],
         [r*sin(a/2),           x4, 0,          a/2 + t -pi/2   ] ]


rep343 = {t:pi/4, a:3*pi/4} # Replacement values in range


# Define conditions for model
cond343 = [a <= pi,  t <= pi/2,  a >= pi - 2*t,  a > 2*t]

# Calculate model, run checks, write output.
p343 = calcModel(m343)
allChecks('p343')
parseLaTeX('p343')




###################################################################################
# 344 animal: a <= pi.  Sensor: t <= pi/2. Condition: a <= pi - 2t & a <= t     #
###################################################################################

m344 = [ [2*r*sin(a/2),                 x2, pi/2 - t/2 + a/2, pi/2            ],
         [r*sin(a/2) - r*cos(x2 + t/2), x2, pi/2 - t/2,       pi/2 - t/2 + a/2],
         [r*sin(a/2),                   x3, t,                t + a/2         ] ]


rep344 = {t:2*pi/8, a:pi/8} # Replacement values in range

# Define conditions for model
cond344 = [a <= pi, t <= pi/2, a <= pi - 2*t, a <= t]

# Calculate model, run checks, write output.
p344 = calcModel(m344)
allChecks('p344')
parseLaTeX('p344')




#######################################################################################
# 345 animal: a <= pi.  Sensor: t <= pi/2. Condition: a <= pi - 2t & t <= a <= 2t   #
#######################################################################################


m345 = [ [2*r*sin(t/2)*sin(x2),         x2, pi/2 + t/2 - a/2, pi/2            ],
         [r*sin(a/2) - r*cos(x2 + t/2), x2, pi/2 - t/2,       pi/2 + t/2 - a/2],
         [r*sin(a/2),                   x3, t,                t + a/2         ] ]

rep345 = {t:2*pi/8, a:pi/2-0.1} # Replacement values in range

# Define conditions for model
cond345 = [a <= pi, t <= pi/2, a <= pi - 2*t, t <= a, a <= 2*t]

# Calculate model, run checks, write output.
p345 = calcModel(m345)
allChecks('p345')
parseLaTeX('p345')



##################################################################################
# 346 animal: a <= pi.  Sensor: t <= pi/2. Condition: a <= pi - 2t &  2t <= a      #
##################################################################################

m346 = [ [2*r*sin(t/2)*sin(x2), x2, pi/2 - t/2, pi/2    ],
         [r*sin(x3),            x3, t,          a/2     ],
         [r*sin(a/2),           x3, a/2,        t + a/2 ] ]


rep346 = {t:1*pi/8, a:pi/2} # Replacement values in range

# Define conditions for model
cond346 = [a <= pi,  t <= pi/2,  a <= pi - 2*t,  2*t <= a]

# Calculate model, run checks, write output.
p346 = calcModel(m346)
allChecks('p346')
parseLaTeX('p346')




##################################################################################################################

####################
## Run tests     ###
####################

# create gas model object
gas = 2*r


# for each model run through every adjacent model. 
# Contains duplicatea but better for avoiding missed comparisons.
# Also contains replacement t->a and a->t just in case. 


allComps = [
['gas', 'p221', {t:2*pi}],
['gas', 'p311', {a:pi}],

['p221', 'gas', {t:2*pi}],
['p221', 'p131', {t:pi}],
['p221', 'p222',{a:3*pi-t}],
['p221', 'p222',{t:3*pi-a}],

['p222', 'p221',{a:3*pi-t}],
['p222', 'p221',{t:3*pi-a}],
['p222', 'p223',{a:4*pi-2*t}],
['p222', 'p223',{t:2*pi-a/2}],
['p222', 'p321',{a:pi}],

['p223', 'p222',{a:4*pi-2*t}],
['p223', 'p222',{t:2*pi-a/2}],
['p223', 'p322',{a:pi}],
['p223', 'p231',{t:pi}],

['p131','p221', {t:pi}],
['p131','p231',{a:2*pi}],

['p231','p223',{t:pi}],
['p231','p232',{a:3*pi-2*t}],
['p231','p232',{t:3*pi/2-a/2}],
['p231','p131',{a:2*pi}],

['p232','p241',{t:pi/2}],
['p232','p233',{a:2*pi-t}],
['p232','p233',{t:2*pi-a}],
['p232','p231',{a:3*pi-2*t}],
['p232','p231',{t:3*pi/2-a/2}],

['p233','p242',{t:pi/2}],
['p233','p232',{t:2*pi-a}],
['p233','p232',{a:2*pi-t}],
['p233','p331',{a:pi}],

['p141','p131', {t:pi/2}],
['p141','p241',{a:2*pi}],

['p241','p141',{a:2*pi}],
['p241','p242',{a:2*pi-t}],
['p241','p242',{t:2*pi-a}],
['p241','p232',{t:pi/2}],

['p242','p241',{a:2*pi-t}], 
['p242','p241',{t:2*pi-a}],
['p242','p243',{t:pi-a/2}],
['p242','p243',{a:2*pi-2*t}],
['p241','p233',{t:pi/2}],

['p243','p242',{t:2*pi-2*a}],
['p243','p242',{a:2*pi-2*t}],
['p243','p343',{a:pi}],

['p311','p321',{t:2*pi}],
['p311','gas',{a:pi}],

['p321','p322',{t:2*pi-a/2}],
['p321','p322',{a:4*pi-2*t}],
['p321','p311',{t:2*pi}],
['p321','p222',{a:pi}],

['p322','p321',{a:4*pi-2*t}],
['p322','p321',{t:2*pi-a/2}],
['p322','p323',{a:2*pi-t}],
['p322','p323',{t:2*pi-a}],
['p322','p223',{a:pi}],

['p323','p322',{t:2*pi-a}],
['p323','p322',{a:2*pi-t}],
['p323','p333',{t:pi}],

['p331','p342',{t:pi/2}],
['p331','p332',{a:t}],
['p331','p332',{t:a}],
['p331','p233',{a:pi}],

['p332','p331',{a:t}],
['p332','p331',{t:a}],
['p332','p341',{t:pi/2}],
['p332','p333',{a:2*t-pi}],
['p332','p333',{t:a/2+pi/2}],

['p333','p332',{t:a/2+pi/2}],
['p333','p332',{a:2*t-pi}],
['p333','p323',{t:pi}],


['p341','p344',{a:pi-2*t}],
['p341','p344',{t:pi/2-a/2}],
['p341','p342',{t:a}],
['p341','p342',{a:t}],
['p341','p332',{t:pi/2}],

['p342','p341',{t:a}],
['p342','p341',{a:t}],
['p342','p345',{t:pi/2-a/2}],
['p342','p345',{a:pi-2*t}],
['p342','p343',{a:2*t}],
['p342','p343',{t:a/2}],
['p342','p331',{t:pi/2}],

['p343','p346',{t:pi/2-a/2}],
['p343','p346',{a:pi-2*t}],
['p343','p342',{a:2*t}],
['p343','p342',{t:a/2}],
['p343','p243',{a:pi}],


['p344','p345',{t:a}],
['p344','p345',{a:t}],
['p344','p341',{t:pi/2-a/2}],
['p344','p341',{a:pi-2*t}],

['p345','p344',{a:t}],
['p345','p344',{t:a}],
['p345','p346',{a:2*t}],
['p345','p346',{t:a/2}],
['p345','p342',{a:pi-2*t}],
['p345','p342',{t:pi/2-a/2}],

['p346','p345',{a:2*t}],
['p346','p345',{t:a/2}],
['p346','p343',{a:pi-2*t}],
['p346','p343',{t:pi/2-a/2}]
]


# List of regions that cover a=0. Should equal 0 when a=0.
zeroRegions = ['p346', 'p345', 'p344', 'p341', 'p332', 'p333', 'p323', 'p322', 'p321', 'p311']

# Run through all the comparisons. Need simplify(). Even together() gives some false negatives.

checkFile = open('/home/tim/Dropbox/phd/Analysis/REM-chapter/checksFile.tex','w')

checkFile.write('All checks evaluated.\nTim Lucas - ' + str(datetime.now()) + '\n')
for i in range(len(allComps)):
        if (eval(allComps[i][0]).subs(allComps[i][2]) - eval(allComps[i][1]).subs(allComps[i][2])).simplify() == 0:
                checkFile.write(str(i) + ': ' + allComps[i][0]+ ' and ' +allComps[i][1]+': OK\n')
        else:
                checkFile.write(str(i) + ': ' + allComps[i][0]+ ' and ' +allComps[i][1]+': Incorrect\n')

for i in range(len(zeroRegions)):
        if eval(zeroRegions[i]).subs({a:0}).simplify() == 0:
                checkFile.write(zeroRegions[i] + ' at a=0: OK\n')
        else:
                checkFile.write(zeroRegions[i] + ' at a=0: Incorrect\n')

checkFile.close()


# And print to terminal
#for i in range(len(allComps)):
#        if not (eval(allComps[i][0]).subs(allComps[i][2]) - eval(allComps[i][1]).subs(allComps[i][2])).simplify() == 0:
#               print allComps[i][0] + ' and ' + allComps[i][1]+': Incorrect\n'

#####################################
## Check some that don't work well ##
#####################################

xRange = np.arange(0,pi/2, 0.01)
y332Range = [p332.subs({r:1,  t:pi/2, a:i}).n() for i in xRange]
plot332 = pl.plot(xRange, y332Range)
pl.savefig('/home/tim/Dropbox/phd/Analysis/REM-chapter/imgs/p332Profile.pdf')
pl.close()

y341Range = [p341.subs({r:1,  t:pi/2, a:i}).n() for i in xRange]
plot341 = pl.plot(xRange, y341Range)
pl.savefig('/home/tim/Dropbox/phd/Analysis/REM-chapter/imgs/p341Profile.pdf')
pl.close()



#pl.savefig('/home/tim/Dropbox/phd/Analysis/REM-chapter/imgs/p221Profile.pdf')
#pl.close()







####################################################################
### Define a a function that calculates your answer.            ####
####################################################################

def calcP(A, T, R): 
	assert (A <= 2*pi and A >= 0), "a is out of bounds. Should be in 0<a<2*pi"
	assert (T <= 2*pi and T >= 0), "s is out of bounds. Should be in 0<s<2*pi"
 	
	if A > pi:
		if A < 4*pi - 2* T:
			p = p243.subs({a:A,  t:T, r:R}).n()
		elif A <= 3*pi -  T:
                        p = p222.subs({a:A,  t:T, r:R}).n()
		else:
                        p = p221.subs({a:A,  t:T, r:R}).n()
	else:
		if A < 4*pi - 2* T:
                        p = p322.subs({a:A,  t:T, r:R}).n()
		else:
                        p = p321.subs({a:A,  t:T, r:R}).n()
        return p


#############################
## Apply to entire grid   ###
#############################

# How many values for each parameter
nParas = 100

# Make a vector for a and s. Make an empty nParas x nParas array. 
# Calculated profile sizes will go in pArray
tVec = np.linspace(0, 2*pi, nParas)
aVec = np.linspace(0, 2*pi, nParas)
pArray = np.zeros((nParas,nParas))

# Calculate profile size for each combination of parameters
for i in range(nParas):
        for j in range(nParas):
                pArray[i][j] = calcP(aVec[i], tVec[j], 1)

# Turn the array upside down so origin is at bottom left.
pImage = np.flipud(pArray)

# Plot and save.
pl.imshow(pImage, interpolation='none', cmap=pl.get_cmap('Blues') )
#pl.show()

pl.savefig('/home/tim/Dropbox/phd/Analysis/REM-chapter/imgs/profilesCalculated.png')



############################
#### Output R function.  ###
############################

# To reduce mistakes, output R function directly from python.
# However, the if statements, which correspond to the bounds of each model, are not automatic.

Rfunc = open('/home/tim/Dropbox/phd/Analysis/REM-chapter/supplementaryRscript.R', 'w')

Rfunc.write("""
# Functions to calculate density.
#
# Tim C.D. Lucas, Elizabeth Moorcroft, Robin Freeman, Marcus J. Rowcliffe, Kate E. Jones.
#
# calcDensity is the main function to calculate density.
# It takes parameters z, alpha, theta, r, animalSpeed, t
# z - The number of camera/acoustic counts or captures.
# alpha - Call width in radians.
# theta - Sensor width in radians.
# r - Sensor range in metres.
# animalSpeed - Average animal speed in metres per second.
# t - Length of survey in sensor seconds i.e. number of sensors x survey duration.
#
# calcAbundance calculates abundance rather than density and requires an extra parameter
# area - In metres squared. The size of the region being examined.


# Internal function to calculate profile width as described in the text
calcProfileWidth <- function(alpha, theta, r){
        if(alpha > 2*pi | alpha < 0) 
		stop('alpha is out of bounds. alpha should be in interval 0<a<2*pi')
        if(theta > 2*pi | theta < 0) 
		stop('theta is out of bounds. theta should be in interval 0<a<2*pi')

	if(alpha > pi){
	        if(alpha < 4*pi - 2*theta){
""" +
'		        p <- ' + str(p243) +
'\n                } else if(alpha <= 3*pi - theta){'  
'\n                        p <- ' + str(p222) +
'\n                } else {'
'\n                        p <- ' + str(p221) +
'\n                }'
'\n        } else {' 
'\n        	if(alpha < 4*pi - 2*theta){'
'\n                        p <- ' + str(p322) +
'\n 		} else {'
'\n                        p <- ' + str(p321) +
'\n                }'
'\n        }'
'\n        return(p)'
'\n}' +
"""
# Calculate a population density. See above for units etc.
calcDensity <- function(z, alpha, theta, r, animalSpeed, t){
        # Check the parameters are suitable.
        if(z <= 0 | !is.numeric(z)) stop('Counts, z, must be a positive number.')
        if(animalSpeed <= 0 | !is.numeric(animalSpeed)) stop('animalSpeed must be a positive number.')
        if(t <= 0 | !is.numeric(t)) stop('Time, t, must be a positive number.')

        # Calculate profile width, then density.
        p <- calcProfileWidth(alpha, theta, r)
        D <- z/{animalSpeed*t*p}
        return(D)
}

# Calculate abundance rather than density.
calcAbundance <- function(z, alpha, theta, r, animalSpeed, t, area){
        if(area <= 0 | !is.numer(area)) stop('Area must be a positive number')
        D <- calcDensity(z, alpha, theta, r, animalSpeed, t)
        A <- D*area
        return(A)
}
"""
)

Rfunc.close()
















