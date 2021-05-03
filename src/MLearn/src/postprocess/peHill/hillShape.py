#import numpy as np
#import matplotlib.pyplot as plt

def profile(yy):
    'Calculate the shape of the periodic hill'

    import numpy as np

    y = np.array(yy)
    x = y * 28

    h = np.zeros(len(x))
    for i in range(len(x)):
        if x[i] > 126.0 :
            x[i] = 252.0 - x[i]

        if (x[i]>=0) and (x[i]<9) :
            h[i] = np.minimum(28., 2.8e+01 + 0.0e+00*x[i] + \
                              6.775070969851e-03*x[i]**2 - 2.124527775800e-03*x[i]**3);
        elif (x[i]>=9) and (x[i]<14) :
            h[i] = 2.507355893131E+01 + 9.754803562315E-01*x[i] - \
                   1.016116352781E-01*x[i]**2 + 1.889794677828E-03*x[i]**3;
        elif (x[i]>=14) and (x[i]<20) :
            h[i] = 2.579601052357E+01 + 8.206693007457E-01*x[i] - \
                   9.055370274339E-02*x[i]**2 + 1.626510569859E-03*x[i]**3;
        elif (x[i]>=20) and (x[i]<30) :
            h[i] = 4.046435022819E+01 - 1.379581654948E+00*x[i] + \
                   1.945884504128E-02*x[i]**2 - 2.070318932190E-04*x[i]**3;
        elif (x[i]>=30) and (x[i]<40) :
            h[i] = 1.792461334664E+01 + 8.743920332081E-01*x[i] - \
                   5.567361123058E-02*x[i]**2 + 6.277731764683E-04*x[i]**3;
        elif (x[i]>=40) and (x[i]<=54) :
            h[i] = np.maximum(0., 5.639011190988E+01 - 2.010520359035E+00*x[i] + \
                              1.644919857549E-02*x[i]**2 + 2.674976141766E-05*x[i]**3);
        elif (x[i]>54) and (x[i]<=126) :
            h[i] = 0;

    hout = h/28.0
    return hout

        
#yy=np.arange(0, 9, 0.01)


#h = myHillShape(yy)

#plt.plot(yy, h)
#plt.axis([0, 9, 0, 1.5])

#plt.show()
