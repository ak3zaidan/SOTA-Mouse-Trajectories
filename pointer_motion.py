import numpy as np
import random
import math
import time


"""

Parameters:

    General:
        cltd - click time delay [s]
        N - number points in trajectory
        n - number of trajectories


    Times:
        to - start time of ith stroke [s]
        tdu - total time duration of stroke n [s]
        tmax - time at which vmax occurs [s]


    Velocities:
        v - current velocity [px/s]
        vmax - max velocity [px/s]
        vthmax - max velocity threshhold [px/s]
        vthmin - min velocity threshhold [px/s]
        vrange - target max velocity [px/s]
        v_end_lim - end velocity limit [px/s]
        v_end - end velocity [px/s]
        v_al - velocity after alignment [px/s]
        verror - velocity error [px/s]
        v_end_error - end velocity error [px/s]


    Angles:
        ang_mul - starting angle multiplier
        ang_s - starting angle [rad]
        ang_e - end angle [rad]
        ang_pcal - calculated angle per evaluated point [rad]
        ang_vec - vector angle for slope calculation [rad]
        ang_top - top angle limit [deg]
        ang_bot - bottom angle limit [deg]
        ang_mid - middle angle adjustment [deg]


    Distances:
        sact - actual distance between starting and  [px]
        D - distance between starting and target point [px]
    
    
    Log Normal Distribution:
        tanf - tangent function
        ti - time interval [s]
        sigma - shape value for normal distribution 
        mu - amplitude of normal distribution (closely linked to velocity)
        mu3 - amplitude of normal distribution for v_end
        sigmas - shape value for normal distribution for v_end
        mumin - minimum value for mu
        mumax - maximum value for mu
        sigmin - minimum value for sigma
        sigmax - maximum value for sigma


    Zero Alignment:
        align - alignment flag
        jcnt - counter for alignment
        save_j - save j value for alignment
        t_align - time at which alignment occurs [s]
        y_inter - y intercept for alignment
        y_inter_s - y intercept for alignment after alignment
        s_alignx - x value for alignment
        s_aligny - y value for alignment
        xdiff - x difference between target and calculated point [px]
        ydiff - y difference between target and calculated point [px]
        end_xdiff - x difference between target and calculated point at end [px]
        end_ydiff - y difference between target and calculated point at end [px]
        nz_cnt - number of non-zero velocities


    Bezier Polynomial:
        degree - degree of bezier curve
        indices - indices for bezier curve fitting
        support_indices - indices for support points
        all_indices - all indices for bezier curve fitting
        support_x - x values for support points
        support_y - y values for support points
        bezier - bezier curve values
        final_poly_x - final x values for bezier curve
        final_poly_y - final y values for bezier curve

    
    Noise:
        NF - noise factor
        tSNR - transversal signal to noise ratio
        pSNR - pperpendicular signal to noise ratio


    Output: 
        move[0] holds the x values
        move[1] holds the y values
        move[2] holds the t values
        move[3] holds the calculated v values
        move[4] holds the calculated v values determined as a listener would + noise

    Usage Examples:

    posis,coals = traj([x_cords, y_cords], self.screen_size[0], self.screen_size[1], 160, 1,True,True)

    mm = []
    for i in range(len(posis[0])):
        trajectory = []
        last_ts = None
        for j in range(len(posis[0][i])):  
            if j == len(posis[0][i])-2:
                mdst =   [int(posis[0][i][j]), int(posis[1][i][j]), int(posis[2][i][j] +thenowstamp) ]
                must =  [int(posis[0][i][j]) + random.randint(1, 3), int(posis[1][i][j]) + random.randint(0, 2), int(posis[2][i][j] +thenowstamp  + random.randint(60, 150))]
                all_datas["md"].push(mdst)
                all_datas["mu"].push(must)

            val = [int(posis[0][i][j]), int(posis[1][i][j]), int(posis[2][i][j] +thenowstamp) ]
            
            if last_ts is not None:
                #time.sleep((posis[2][i][j]-last_ts) / 1000)
                profile.add_offset(posis[2][i][j]-last_ts)
                
            last_ts = posis[2][i][j]
            
            all_datas["mm"].push(val)

            trajectory.append(val)
            mm.append(val)
            
        md.append(mdst)
        mu.append(must)

"""

def biased_random(mean, std_dev, lower_bound, upper_bound):
    while True:
        value = np.random.normal(mean, std_dev)
        if lower_bound <= value <= upper_bound:
            return value
        

#def traj(xo,yo,xe,ye,screenmax_x,screenmax_y,n,frequ):
def traj(coordinates, screenmax_x, screenmax_y, frequ, sv, coal, remove_v0):
    starters = time.time()

    #--------------------------------------------------------------------------------------------------#
    #Setting initial values:
    n = len(coordinates[0])
    global debug
    debug = 0

    n = n-1
    screenmin_y = 0
    screenmin_x = 0

    xf = np.zeros(n)
    yf = np.zeros(n)
    xco = coordinates[0]
    yco = coordinates[1]

    xe= xco[n]
    ye= yco[n]

    xo= xco[0]
    yo= yco[0]


    xf = np.diff(xco[:n+1]) 
    yf = np.diff(yco[:n+1])
    
    start_time= np.zeros(n)
    elapsed_time= np.zeros(n)
    end_time= np.zeros(n)
    vrange= np.zeros(n)
    vthmin= np.zeros(n)
    vthmax= np.zeros(n)
    v_end_lim= np.zeros(n)
    ang_mul = np.zeros(n)
    
    s = np.zeros(n)
    sact = np.zeros(n)
    screen_diag = math.sqrt(screenmax_x**2+screenmax_y**2)
    sqrt_2pi = math.sqrt(2*math.pi)

    for i in range(0, n):

        if xf[i] == 0:
            xf[i] =  1*10**-5
        if yf[i] == 0:
            yf[i] = 1*10**-5
    
        
        vthmax[i] = np.random.uniform(28654,29124) # cm/s max velocity 
        vthmin[i] = np.random.uniform(1523,1979)

        s[i] = math.sqrt(yf[i]**2+xf[i]**2)
        sact[i] = s[i]
        vrange[i] = (s[i]*(vthmax[i]-vthmin[i])/screen_diag)
        ang_mul[i] = s[i]/screen_diag*(3/10)
        
        if ang_mul[i] < 1:
            ang_mul[i] = 1
        
        
        v_end_lim[i] = 10
    

         
    mumin = -5
    mumax = -0.1
    sigmax = 3
    sigmin = 0.1

    NF = 500000
    tSNR = 80*sv
    pSNR = 11*sv

    k = 0
    posaddx=0
    posaddy=0
    posaddypost=0
    posaddxpost=0
    xerrdiff = 0.5

    ang_top = 20
    ang_bot = 11
    ang_mid = 24

    end = np.zeros(n)
    cltd= np.zeros(n)
    tdu= np.zeros(n)
    logger= np.zeros(n)
    xfcal= np.zeros(n)
    yfcal= np.zeros(n)
    ang_pcal= np.zeros(n)
    erf = np.zeros(n)
    tanf = np.zeros(n)
    ti = np.zeros(n)
    sigma = np.zeros(n)
    to = np.zeros(n+1)
    N = np.zeros(n)
    D  = np.zeros(n)
    mu3 = np.zeros(n)
    tmax = np.zeros(n)
    vmax = np.zeros(n)
    ang_s = np.zeros(n)
    ang_e = np.zeros(n)
    v_end = np.zeros(n)
    sxpre_al = np.zeros(n)
    sypre_al = np.zeros(n)
    xdiff = np.zeros(n)
    ydiff = np.zeros(n)
    end_xdiff = np.zeros(n)
    end_ydiff = np.zeros(n)
    t_delta = np.zeros(n)
    verror = np.zeros(n)

    degree = np.zeros(n, dtype=int)  # Initialize degree as an integer array of size n
    indices =  np.zeros((n,2))
    #--------------------------------------------------------------------------------------------------#
    #--------------------------------------------------------------------------------------------------#
    

    # Function to determine new random values for mu3 and sigma
    #--------------------------------------------------------------------------------------------------#
    def newrandos():
        global itcnt
        itcnt = 0
        mu3[i] = ((mumax-mumin)/(500))*vrange[i]-3.0
        if mu3[i] > -0.1:
            mu3[i] = -0.1
        sigma[i] = math.sqrt(np.abs(math.log(tdu[i]/3.1)-mu3[i]))



        #mu3[i] = random.uniform(mumin, mumax)  
        #sigma[i] = random.uniform(sigmin, sigmax)
        sigma[i], mu3[i] = errorf([sigma[i],mu3[i]])

        return mu3[i], sigma[i]
    #--------------------------------------------------------------------------------------------------#
    #--------------------------------------------------------------------------------------------------#

    # Function to determine the error of the current mu3 and sigma values and refine them
    #--------------------------------------------------------------------------------------------------#
    def errorf(params):
        sigmas = params[0]
        mu3cal = params[1]

        for k in range(100):
            tmax[i] = to[i] + math.e**((-(sigmas)**2) + mu3cal )
            erf[i]= math.erf((math.sqrt(2)*(mu3cal-math.log(tdu[i])))/(2*sigmas))
            tanf[i] = math.atan2((xf[i]*math.tan((ang_s[i])/(2))-yf[i]),(xf[i]+yf[i]*math.tan(ang_s[i]/2)))
            term0 = (erf[i] - 1)
            if term0 == 0:
                term0= 1*10**-14

            term1 = ang_s[i] * (term0)  
            term2 = term1 - ang_s[i] * erf[i] - ang_s[i] - 4 * tanf[i]
            logger[i]= tmax[i]-to[i]
            if logger[i] <= 0:
                logger[i] = 1*10**-20

            D[i] = xf[i] * term2 / ((math.sin(ang_s[i]) + math.sin(2 * tanf[i])) * (term0))
            vmax[i] = ((D[i])/((logger[i])*sigmas*sqrt_2pi))*math.e**-(((math.log(logger[i])-mu3cal)**2)/(2*sigmas**2))
            xerror = 0
            yerror= 0
            verror[i] = 0
            v_end_error = 0
            ang_e[i] = (ang_s[i]*erf[i]+ang_s[i]+4*tanf[i])/(term0)
            ang_pcal[i] = ang_s[i] + ((ang_e[i]-ang_s[i])/(2))*(1+(math.erf((math.log(tdu[i])-mu3cal)/(sigmas*math.sqrt(2)))))
            xfcal[i] = ((D[i])/(ang_e[i]-ang_s[i]))*(math.sin(ang_pcal[i])-math.sin(ang_s[i]))
            yfcal[i] = ((D[i])/(ang_e[i]-ang_s[i]))*(-math.cos(ang_pcal[i])+math.cos(ang_s[i]))
            v_end[i] = ((D[i])/((tdu[i])*sigmas*sqrt_2pi))*math.e**(-(((math.log(tdu[i])-mu3cal)**2)/(2*sigmas**2)))
            
            v_end_error = math.sqrt((v_end[i]) ** 2)
            if v_end_error <= v_end_lim[i]:
                v_end_error = 0

            xerror = math.sqrt((xfcal[i] - xf[i]) ** 2)
            yerror = math.sqrt((yfcal[i] - yf[i]) ** 2)
            verror[i] = math.sqrt((vmax[i] - vrange[i]) ** 2)
            if (abs(vmax[i] - vrange[i]) < 500):
                verror[i] = 0 

            error_sum = xerror + yerror + verror[i] + v_end_error
            if error_sum < xerrdiff:
                break

            log_end = (v_end[i]*tdu[i]*sigmas*sqrt_2pi)
            if log_end <= 0:
                log_end = 1*10**-20
            mu3end = math.log(tdu[i])-math.sqrt(2*sigmas**2*log_end)
            mu3max = math.log(logger[i])-math.sqrt(2*sigmas**2*math.log(vmax[i]*logger[i]*sigmas*sqrt_2pi))
            mu3cal= (mu3end+mu3max)/2

            if mu3cal < math.log(tmax[i]-to[i]):
                mu3cal = math.log(tmax[i]-to[i])+0.15

            sigmas = math.sqrt(mu3cal-math.log(tmax[i]-to[i]))
            if sigmas <= 0:
                sigmas = 0.1

            if (verror[i] != 0):

                if ((verror[i] > verror[i-1]) and (verror[i-1] == 0)):
                    mu3cal = mu3cal - 0.01
                    sigmas = sigmas + 0.01
                else:
                    mu3cal = mu3cal + 0.01
                    sigmas = sigmas - 0.01
            
            tdu[i] = tdu[i] + 0.06
            debug = 0


            debug = 0

        return sigmas, mu3cal
    #--------------------------------------------------------------------------------------------------#
    #--------------------------------------------------------------------------------------------------#


    # Initialization of array values in dependency of N
    #--------------------------------------------------------------------------------------------------#
    #cltd = np.random.uniform(0.1, 0.35, n)
    for i in range(0, n):
        degree[i] =  np.random.randint(4, 10)
        tdu[i] = (s[i]/vrange[i]) * 3
        N[i] = round((tdu[i]+cltd[i])*frequ)
        N[i] = int(N[i])
        N= N.astype(int)
    v = [np.zeros(N[i]) for i in range(n)]
    v_al = [np.zeros(N[i]) for i in range(n)]
    t = [np.zeros(N[i]) for i in range(n)]
    ang_p = [np.zeros(N[i]) for i in range(n)]
    sxpre = [np.zeros(N[i]) for i in range(n)]
    sypre = [np.zeros(N[i]) for i in range(n)]
    xpos =  [np.zeros(N[i]) for i in range(n)]
    ypos = [np.zeros(N[i]) for i in range(n)]
    final_poly_x= [np.zeros(N[i]) for i in range(n)]
    final_poly_y= [np.zeros(N[i]) for i in range(n)]
    #--------------------------------------------------------------------------------------------------#
    #--------------------------------------------------------------------------------------------------#


    #Calculation of time values
    #--------------------------------------------------------------------------------------------------#
    for i in range(0, n):
        starters = time.time()
        if i == 0:
            to[i] = 0
        else:
            to[i] = to[i-1] + tdu[i-1]+cltd[i-1]
        for j in range (0, N[i]):
            if j == 0:
                t[i][j] = to[i]
            else:
                t[i][j] = t[i][j-1] + (tdu[i]+cltd[i])/N[i] 
    
        ti[i] = np.abs(np.random.normal(0.3, 0.05))
        while ti[i] < 0.1:
            ti[i]= np.abs(np.random.normal(0.3, 0.05))
        
        t_delta[i]=(tdu[i]+cltd[i])/N[i] 

        tfin= to[n-1]+tdu[n-1]  
     #--------------------------------------------------------------------------------------------------#
     #--------------------------------------------------------------------------------------------------#



     # Calculation for the Starting angle ang_s:
     #--------------------------------------------------------------------------------------------------#
        if (i < n-1):
            if math.atan2(yco[i+1]-yco[i],xco[i+1]-xco[i]) >=  math.atan2(yco[i+2]-yco[i],xco[i+2]-xco[i]):
                ang_s[i]=  biased_random((math.atan2(yco[i+1]-yco[i],xco[i+1]-xco[i])+ang_mul[i]*ang_mid*math.pi/180), 0.99 , math.atan2(yco[i+1]-yco[i],xco[i+1]-xco[i])-ang_mul[i]*ang_bot*math.pi/180, math.atan2(yco[i+1]-yco[i],xco[i+1]-xco[i])+ang_mul[i]*ang_top*math.pi/180 )
            elif math.atan2(yco[i+1]-yco[i],xco[i+1]-xco[i]) <  math.atan2(yco[i+2]-yco[i],xco[i+2]-xco[i]):
                ang_s[i]=  biased_random((math.atan2(yco[i+1]-yco[i],xco[i+1]-xco[i])-ang_mul[i]*ang_mid*math.pi/180), 0.99, math.atan2(yco[i+1]-yco[i],xco[i+1]-xco[i])-ang_mul[i]*ang_top*math.pi/180, math.atan2(yco[i+1]-yco[i],xco[i+1]-xco[i])+ang_mul[i]*ang_bot*math.pi/180 )
        else:
            if math.atan2(yco[i]-yco[i-1],xco[i]-xco[i-1]) >=  math.atan2(yco[i+1]-yco[i-1],xco[i+1]-xco[i-1]):
                ang_s[i] = biased_random((math.atan2(yco[i+1]-yco[i],xco[i+1]-xco[i])+ang_mul[i]*ang_mid*math.pi/180), 0.99, math.atan2(yco[i+1]-yco[i],xco[i+1]-xco[i])-ang_mul[i]*ang_bot*math.pi/180, math.atan2(yco[i+1]-yco[i],xco[i+1]-xco[i])+ang_mul[i]*ang_top*math.pi/180 )
            elif math.atan2(yco[i]-yco[i-1],xco[i]-xco[i-1]) <  math.atan2(yco[i+1]-yco[i-1],xco[i+1]-xco[i-1]):
                ang_s[i] = biased_random((math.atan2(yco[i+1]-yco[i],xco[i+1]-xco[i])-ang_mul[i]*ang_mid*math.pi/180), 0.99 , math.atan2(yco[i+1]-yco[i],xco[i+1]-xco[i])-ang_mul[i]*ang_top*math.pi/180, math.atan2(yco[i+1]-yco[i],xco[i+1]-xco[i])+ang_mul[i]*ang_bot*math.pi/180 )

        if (ang_s[i] >= math.atan2(yco[i+1]-yco[i],xco[i+1]-xco[i])+ang_mul[i]*80*math.pi/180) or (ang_s[i] <= math.atan2(yco[i+1]-yco[i],xco[i+1]-xco[i])-ang_mul[i]*80*math.pi/180):
            ang_s[i] = math.atan2(yco[i+1]-yco[i],xco[i+1]-xco[i])+ang_mul[i]*ang_mid*math.pi/180


        while ang_s[i] <= -math.pi:
            ang_s[i] = (ang_s[i] + 2*math.pi)

        while ang_s[i] >= math.pi:
            ang_s[i] = (ang_s[i]-2*math.pi)
        enders = time.time()
    #--------------------------------------------------------------------------------------------------#
    #--------------------------------------------------------------------------------------------------#



    # Determining random values:
    #--------------------------------------------------------------------------------------------------#
        start_time[i] = time.time()
        newrandos()
    #--------------------------------------------------------------------------------------------------#
    #--------------------------------------------------------------------------------------------------#



    #Re-Calculation of t's
    #--------------------------------------------------------------------------------------------------#

    for i in range(0, n):
        if i == 0:
            to[i] = 0
        else:
            to[i] = to[i-1] + tdu[i-1]+cltd[i-1]
        for j in range (0, N[i]):
            if j == 0:
                t[i][j] = to[i]
            else:
                t[i][j] = t[i][j-1] + (tdu[i]+cltd[i])/N[i] 
    
        ti[i] = np.abs(np.random.normal(0.3, 0.05))
        while ti[i] < 0.1:
            ti[i]= np.abs(np.random.normal(0.3, 0.05))
        
        t_delta[i]=(tdu[i]+cltd[i])/N[i] 

        tfin= to[n-1]+tdu[n-1]  
     #--------------------------------------------------------------------------------------------------#
     #--------------------------------------------------------------------------------------------------#



    # Intermediate angle ang_p and position calculation sxpre and sypre
    #--------------------------------------------------------------------------------------------------#
        align = 0
        jcnt = 0
        save_j = 0
        tmax[i] = to[i] + math.e**((-(sigma[i])**2) + mu3[i] )
        for j in range (0, N[i]):
            logs = t[i][j]-to[i]
            if logs <= 0:
                logs = 1*10**-20

            ang_p[i][j] = ang_s[i] + ((ang_e[i]-ang_s[i])/(2))*(1+(math.erf((math.log(logs)-mu3[i])/(sigma[i]*math.sqrt(2)))))
            sxpre[i][j] = ((D[i])/(ang_e[i]-ang_s[i]))*(math.sin(ang_p[i][j])-math.sin(ang_s[i]))+posaddx
            sypre[i][j] = ((D[i])/(ang_e[i]-ang_s[i]))*(-math.cos(ang_p[i][j])+math.cos(ang_s[i]))+posaddy 
            v[i][j] = ((D[i])/((logs)*sigma[i]*sqrt_2pi))*math.e**(-(((math.log(logs)-mu3[i])**2)/(2*sigma[i]**2)))
            if ((t[i][j]>(tmax[i]))& (t[i][j]< to[i]+tdu[i])& (align == 0)):
                slope = (- (D[i] / ((logs) **2 * sigma[i] * sqrt_2pi)) * math.exp(-((math.log(logs) - mu3[i]) **2) / (2 * sigma[i] ** 2)) * (1 + (math.log(logs) - mu3[i]) / (sigma[i]**2)))
                
                y_inter= (v[i][j])-(slope*logs)
               
                v_al[i][j] = slope * ((tdu[i]-2*((tdu[i]+cltd[i])/N[i]))) + y_inter
                if j > 0:
                    if ((v_al[i][j-1]<0) & (v_al[i][j]>0)):
                        slope =  (- (D[i] / ((logs) **2 * sigma[i] * math.sqrt(2 * math.pi))) * math.exp(-((math.log(logs) - mu3[i]) **2) / (2 * sigma[i] ** 2)) * (1 + (math.log(logs) - mu3[i]) / (sigma[i]**2)))
                        align = 1
                        t_align = t[i][j]-to[i]     
                        y_inter= (v[i][j])-(slope*logs)
                        y_inter_s = v[i][j]
                        s_alignx= sxpre[i][j]
                        s_aligny= sypre[i][j]
                        
                        xdiff[i] = coordinates[0][i+1]-sxpre[i][j]-coordinates[0][0]
                        ydiff[i] = coordinates[1][i+1]-sypre[i][j]-coordinates[1][0]
                        ang_vec = math.atan2(ydiff[i],xdiff[i])
            if (align == 1): 

                v[i][j] = slope * (t[i][j]-to[i]) + y_inter

                sxpre[i][j] = math.cos(ang_vec)*(slope*((((t[i][j]-to[i])-t_align)**2)/2)+y_inter_s*((t[i][j]-to[i])-t_align))+s_alignx
                sypre[i][j] = math.sin(ang_vec)*(slope*((((t[i][j]-to[i])-t_align)**2)/2)+y_inter_s*((t[i][j]-to[i])-t_align))+s_aligny
            
                if v[i][j] <= 0:
                    v[i][j] = 0
                    save_j = j-1 
                    jcnt = jcnt + 1
                    if (save_j != 0) and (jcnt == 1) :
                        sxpre_al[i]= sxpre[i][j]
                        sypre_al[i]= sypre[i][j]

                        sxpre[i][j-1] = ((sxpre_al[i]+sxpre[i][j-2])/2)*np.random.uniform(0.99995, 1.00005)   
                        sypre[i][j-1] = ((sypre_al[i]+sypre[i][j-2])/2)*np.random.uniform(0.99995, 1.00005)
                        end[i]= j
                        end[i] = end[i].astype(int)

                    
                    sxpre[i][j] = sxpre_al[i]
                    sypre[i][j] = sypre_al[i]
                    
                    if (save_j != 0) and (jcnt == 2) :
                        sxpre[i][j] = sxpre_al[i]
                        sypre[i][j] = sypre_al[i]
                        sxpre[i] = sxpre[i][:j+1]
                        sypre[i] = sypre[i][:j+1]
                        v[i] = v[i][:j+1]
                        v_al[i] = v_al[i][:j+1]
                        tdu[i]=t[i][j]-to[i]
                        t[i] = t[i][:j+1]
                        ang_p[i] = ang_p[i][:j+1]
                        N[i]=j+1
                        break
    #--------------------------------------------------------------------------------------------------#
    #--------------------------------------------------------------------------------------------------#



    # Calculation of final positions sxpre and sypre
    #--------------------------------------------------------------------------------------------------#

        end_xdiff[i]= sxpre[i][N[i]-1]-(coordinates[0][i+1]-coordinates[0][0])
        end_ydiff[i]= sypre[i][N[i]-1]-(coordinates[1][i+1]-coordinates[1][0])
        nz_cnt = np.sum(v[i] > 0)
        for j in range(1, N[i]):
            if sxpre[i][j] < (coordinates[0][i+1] - coordinates[0][0]):
                sxpre[i][j] = sxpre[i][j] + j * end_xdiff[i] / nz_cnt
            else:
                sxpre[i][j] = sxpre[i][j] - j * end_xdiff[i] / nz_cnt


            if sypre[i][j] < (coordinates[1][i+1] - coordinates[1][0]):
                sypre[i][j] = sypre[i][j] + j * end_ydiff[i] / nz_cnt
            else:
                sypre[i][j] = sypre[i][j] - j * end_ydiff[i] / nz_cnt


        end_xdiff[i] = sxpre[i][N[i] - 1] - (coordinates[0][i+1] - coordinates[0][0])
        end_ydiff[i] = sypre[i][N[i] - 1] - (coordinates[1][i+1] - coordinates[1][0])
 
        posaddx=sxpre[i][N[i]-1]
        posaddy=sypre[i][N[i]-1]
  
        
        xpos =  [np.zeros(N[i]) for i in range(n)]
        ypos = [np.zeros(N[i]) for i in range(n)]
        final_poly_x= [np.zeros(N[i]) for i in range(n)]
        final_poly_y= [np.zeros(N[i]) for i in range(n)]
        move = [ [np.zeros(N[i]) for i in range(n)] for _ in range(5) ]
        bezier = [[np.zeros(N[i]) for _ in range(2)] for i in range(n)]
        if N[i]/15 < 2:
            degree[i] = 2
        else:
            degree[i] =  np.random.randint(N[i]/15, N[i]/7)
            if degree[i] < 2:
                degree[i] = 2
        support_x =  [np.zeros(degree[i]+1) for i in range(n)]
        support_y =  [np.zeros(degree[i]+1) for i in range(n)]
        support_indices = [np.zeros(degree[i]) for i in range(n)]
        all_indices = [np.zeros(degree[i]+2) for i in range(n)]
    #--------------------------------------------------------------------------------------------------#
    #--------------------------------------------------------------------------------------------------#


    # Function for Bezier curve fitting on basis of original curve
    #--------------------------------------------------------------------------------------------------#
    def bezier_curve(x_support, y_support, indices, sxpre, sypre, index):
        y_values = np.zeros(len(sxpre))
        x_values = np.zeros(len(sxpre))
        for i in range(len(x_support)-1):
            #f = 0.3*(i/(len(x_support)-1))+0.3
            f = 0.5
            y2 = y_support[i + 1]
            x2 = x_support[i + 1]
            y0 = y_support[i]
            x0 = x_support[i]
            if i == 0:
                sign = np.random.choice([-1, 1])
            else:
                sign = sign * -1
            alpha = np.arctan2(y2 - y0, x2 - x0)
            distance = np.sqrt((x2 - x0) ** 2 + (y2 - y0) ** 2)

            hpre = sign * distance / np.random.uniform(30, 35)
            hfac = 2
            h = hpre * hfac
            
            x3pre = f * distance
            x3 = x0 + x3pre * np.cos(alpha) - h * np.sin(alpha)
            y3 = y0 + x3pre * np.sin(alpha) + h * np.cos(alpha)



            for j in range(indices[i],indices[i+1]) :

                t = np.sqrt((sxpre[j] - x0) ** 2 + (sypre[j] - y0) ** 2) / distance

                    
                x = (1 - t) ** 2 * x0 + 2 * (1 - t) * t * x3 + t ** 2 * x2
                y = (1 - t) ** 2 * y0 + 2 * (1 - t) * t * y3 + t ** 2 * y2
                y_values[j]=y
                x_values[j]=x

                if j == indices[len(indices)-1]-1:

                    t = np.sqrt((sxpre[j] - x0) ** 2 + (sypre[j] - y0) ** 2) / distance
                   
                    
                    x = (1 - t) ** 2 * x0 + 2 * (1 - t) * t * x3 + t ** 2 * x2
                    y = (1 - t) ** 2 * y0 + 2 * (1 - t) * t * y3 + t ** 2 * y2
                    y_values[j+1]=y
                    x_values[j+1]=x  
           
        return x_values, y_values
    #--------------------------------------------------------------------------------------------------#
    #--------------------------------------------------------------------------------------------------#
        # Calculation of coalescing events
    #--------------------------------------------------------------------------------------------------#
    def add_coalesced_events(move, event_rate=60):
        """
        Simuliert Coalesced Pointer-Events aus 100Hz-Rohdaten.
        
        Bei echten OS werden Mausereignisse mit hoher Rate gesammelt und dann
        in regelmäßigen Intervallen (z.B. alle 16,67ms bei 60Hz) gebündelt an
        die Anwendung gesendet.
        
        Args:
            move (list of np.ndarray): Liste mit 5 Dimensionen: x, y, t, v, v_noise
            event_rate (float): Ziel-Rate für Event-Dispatching in Hz (default 60)
        
        Returns:
            move_c (list): Liste mit 5 Dimensionen, wobei jedes Element eine Liste von
                        Trajektorien ist und jede Trajektorie eine Liste von Events.
                        Jedes Event enthält alle Punkte, die in einem Zeitfenster
                        gesammelt wurden.
        """
        num_traj = len(move[0])  # Anzahl der Trajektorien
        interval_ms = 1000.0 / event_rate  # Zeitfenster in Millisekunden
        
        # Ergebnis-Container: für jede Dimension d und Trajektorie i eine Liste von Events
        move_c = [[[] for _ in range(num_traj)] for _ in range(5)]
        
        # Für jede Trajektorie:
        for i in range(num_traj):
            # Zeitstempel als Referenz für Gruppierung verwenden
            times = move[2][i]
            last_dispatch_time = 0.0
            
            # Punkte nach Zeitfenstern gruppieren
            current_event_indices = []
            
            for j in range(len(times)):
                t_j = times[j]
                
                # Wenn der aktuelle Punkt innerhalb des aktuellen Zeitfensters liegt
                if t_j < last_dispatch_time + interval_ms:
                    current_event_indices.append(j)
                else:
                    # Zeitfenster abschließen, wenn Punkte vorhanden
                    if current_event_indices:
                        for d in range(5):
                            # Alle Punkte dieses Zeitfensters als ein Event speichern
                            move_c[d][i].append(np.array([move[d][i][idx] for idx in current_event_indices]))
                    
                    # Neues Zeitfenster beginnen
                    current_event_indices = [j]  # Aktuellen Punkt zum neuen Fenster hinzufügen
                    
                    # Dispatching-Zeit auf Vielfaches des Intervalls setzen
                    # (OS dispatcht in festen Intervallen, nicht kontinuierlich)
                    last_dispatch_time = (t_j // interval_ms) * interval_ms
            
            # Letztes Zeitfenster abschließen, falls Punkte übrig
            if current_event_indices:
                for d in range(5):
                    move_c[d][i].append(np.array([move[d][i][idx] for idx in current_event_indices]))
        
        return move_c
    #--------------------------------------------------------------------------------------------------#
    #--------------------------------------------------------------------------------------------------#



    # Removing 0 Velocity points
    #--------------------------------------------------------------------------------------------------#
    def filter_v0(move):
        """
        Entfernt aus move alle Samples j, bei denen die Distanz zum vorherigen Punkt 0 ist.
        Dabei entfällt auch der erste Sample j=0, da dort keine Distanz definiert ist.

        Args:
            move (list of np.ndarray): 
                move[d].shape = (num_traj, num_samples),
                d=0:x,1:y,2:t,3:v,4:v_noise

        Returns:
            filtered (list of list of np.ndarray):
                filtered[d][i] ist ein 1D-Array aller verbleibenden move[d][i,j]
                nach Entfernung aller Stationärpunkte.
        """
        num_traj = len(move[0])
        # Hier sammeln wir für jede Dimension d eine Liste von gefilterten Trajektorien
        filtered = [[] for _ in range(5)]

        for i in range(num_traj):
            # 1) Differenzen nur für x und y
            xs = move[0][i]
            ys = move[1][i]
            dx = np.diff(xs)
            dy = np.diff(ys)

            # 2) Maske: True, wenn Bewegung stattgefunden hat
            mask = (dx**2 + dy**2) != 0  # Länge = num_samples-1

            # 3) Anwenden auf jede Dimension:
            #    Für alle d droppen wir den ersten Sample (da keine diff definiert),
            #    und wenden dann dieselbe Maske an.
            for d in range(5):
                arr = move[d][i][1:]      # drop j=0
                filtered[d].append(arr[mask])

        return filtered
    #--------------------------------------------------------------------------------------------------#
    #--------------------------------------------------------------------------------------------------#   


    # add delays for clicks
    #--------------------------------------------------------------------------------------------------#
    def add_click_delay(move, i, frequ):
        """
        Fügt am Ende einer Trajektorie eine Verzögerung für den Klick hinzu.
        
        Args:
            move (list): Liste mit 5 Dimensionen [x, y, t, v, v_noise]
            i (int): Index der aktuellen Trajektorie
            frequ (float): Sampling-Frequenz in Hz
        
        Returns:
            list: Das aktualisierte move-Array mit hinzugefügtem Delay
        """
        # Letzten Punkt der Trajektorie nehmen
        j = len(move[0][i]) - 1
        safe_x = move[0][i][j]
        safe_y = move[1][i][j]
        safe_t = move[2][i][j]
        
        # Zufällige Verzögerung generieren
        random_delay = np.random.uniform(0.5, 1.13)
        num_extra_points = int(random_delay * frequ)
        dtdelay = 1.0 / frequ
        
        # Zusätzliche Punkte für die Verzögerung erstellen
        extra_x = np.full(num_extra_points, safe_x)
        extra_y = np.full(num_extra_points, safe_y)
        extra_t = safe_t + dtdelay * np.arange(1, num_extra_points + 1)
        extra_v = np.zeros(num_extra_points)
        extra_v_al = np.zeros(num_extra_points)
        
        # Arrays erweitern
        move[0][i] = np.concatenate([move[0][i], extra_x])
        move[1][i] = np.concatenate([move[1][i], extra_y])
        move[2][i] = np.concatenate([move[2][i], extra_t])
        move[3][i] = np.concatenate([move[3][i], extra_v])
        move[4][i] = np.concatenate([move[4][i], extra_v_al])
        
        # Zeitoffset für die nächste Trajektorie anpassen, falls vorhanden
        num_traj = len(move[0])
        if i < num_traj - 1:
            time_offset = move[2][i][-1] + dtdelay
            move[2][i + 1] = move[2][i + 1] + time_offset
        
        return move

    # Calculation of support points and final polynom trajectory
    #--------------------------------------------------------------------------------------------------#
    for i in range (0, n):
        indices[i] = [0, len(sxpre[i])-1]
        support_indices[i] = random.sample(range(2, len(sxpre[i])-2), degree[i] - 1)
        all_indices[i] = np.concatenate((indices[i], support_indices[i]))
        all_indices[i] = np.sort(all_indices[i]).astype(int)

        for k in range(len(all_indices[i])):
            support_x[i][k] = sxpre[i][all_indices[i][k]] 
            support_y[i][k] = sypre[i][all_indices[i][k]]

        support_x[i] = np.array(support_x[i])
        support_y[i] = np.array(support_y[i])    
        bezier[i]= bezier_curve(support_x[i], support_y[i], all_indices[i] ,sxpre[i],sypre[i],i )

        final_poly_y[i] = bezier[i][1]
        final_poly_x[i] = bezier[i][0]
        ypos[i]=bezier[i][1]
        xpos[i]=bezier[i][0]

        if i > 0:
            for j in range (0, np.argmax(v[i] > 0)):
                xpos[i][j]=xpos[i-1][N[i-1]-1] - coordinates[0][0]
                ypos[i][j]=ypos[i-1][N[i-1]-1] - coordinates[1][0]
        
        for j in range (0, len(ypos[i])):
            final_poly_y[i][j] = final_poly_y[i][j]+yo
            final_poly_x[i][j] = final_poly_x[i][j]+xo
        
            ypos[i][j]=ypos[i][j]
            xpos[i][j]=xpos[i][j]
        end_time[i] = time.time()
        elapsed_time[i] = end_time[i] - start_time[i]
            
    k=0
    tmaxi = sum(elapsed_time)

    move = [xpos, ypos, t, v, v]
    #--------------------------------------------------------------------------------------------------#
    #--------------------------------------------------------------------------------------------------#
            
        

    # Noise calculation
    #--------------------------------------------------------------------------------------------------#
    for i in range (0, n):
        for j in range (0, N[i]):
            if (move[3][i][j] > vrange[i]/NF) and ((j > 5)and (j < N[i]-5)):
            #if (move[3][i][j] > vrange[i]/NF):

                dx = move[0][i][j+1] - move[0][i][j]
                dy = move[1][i][j+1] - move[1][i][j]
                
                norm = np.sqrt(dx**2 + dy**2)
                direction_x = dx / norm
                direction_y = dy / norm

                noisex = np.random.normal(0, 0.035)
                noisey = np.random.normal(0, 0.035)

                tangential_noise = noisex * direction_x + noisey * direction_y
                perpendicular_noise = noisex * -direction_y + noisey * direction_x
                
                #vel_factor = 0.8*move[3][i][j+1]/vrange[i]
                if move[3][i][j] > vrange[i]/ NF:
                    vel_factor = math.sin(move[3][i][j+1]/vrange[i])
                elif move[3][i][j] < vrange[i]/ NF:
                    vel_factor = 1*10**-18
                
                
                final_noise_x = vel_factor * tSNR * 0.008 *s[i] * tangential_noise * direction_x + pSNR * 0.0040 *s[i]* perpendicular_noise * -direction_y
                final_noise_y =  vel_factor * tSNR * 0.008 *s[i]* tangential_noise * direction_y + pSNR * 0.0040 *s[i] * perpendicular_noise * direction_x
                
                move[0][i][j]=round(final_noise_x+move[0][i][j])
                move[1][i][j]=round(final_noise_y+move[1][i][j])
    #--------------------------------------------------------------------------------------------------#
    #--------------------------------------------------------------------------------------------------#



    # Bounds and velocity check
    #--------------------------------------------------------------------------------------------------#
            if move[3][i][j] <= 0:
                  
                move[0][i][j]=move[0][i][j]
                move[1][i][j]=move[1][i][j]
                
            
            if (move[1][i][j] >= screenmax_y) or (move[1][i][j] <= screenmin_y):
                move[1][i][j] = -1
            
            if (move[0][i][j] >= screenmax_x):
                move[0][i][j]= screenmax_x
            
            if (move[0][i][j] < screenmin_x):
                move[0][i][j]= screenmin_x
            
            move[0][i][j]=round(move[0][i][j])
            move[1][i][j]=round(move[1][i][j])
            #move[2][i][j]=round((move[2][i][j]*1000)+time.time())
            #move[2][i][j]=int(round(move[2][i][j]*1000+(time.time()*1000)))
            move[2][i][j]=round((move[2][i][j]*1000))
            
            if j > 0:
                move[4][i][j] = (((move[0][i][j]-move[0][i][j-1]) ** 2 + (move[1][i][j]-move[1][i][j-1])  ** 2) ** 0.5) / (t_delta[i] )
    #--------------------------------------------------------------------------------------------------#
    #--------------------------------------------------------------------------------------------------#


    if remove_v0 == True:
        move = filter_v0(move)
    
    for i in range(n):
        move = add_click_delay(move, i, frequ)




    #--------------------------------------------------------------------------------------------------#
    #--------------------------------------------------------------------------------------------------#  


    # Return the move array
    #--------------------------------------------------------------------------------------------------#
    move = np.array(move, dtype=object).tolist()

    if coal == True:
        movecoal = add_coalesced_events(move, event_rate=60)
        movecoal = np.array(movecoal, dtype=object).tolist()
        return move, movecoal
    else:
        return move
    #--------------------------------------------------------------------------------------------------#
    #--------------------------------------------------------------------------------------------------#
