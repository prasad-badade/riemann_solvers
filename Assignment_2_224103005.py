import sys
orig_stdout = sys.stdout
fp = open('Assignment_2_Output.txt', 'w')
sys.stdout = fp

def f_l(p_star,gamma,pl,rho_l,a_l):    
    if p_star>pl :
        fl=(p_star-pl)*(((2/((gamma+1)*rho_l))/(p_star+(pl*(gamma-1)/(gamma+1)))))**(0.5)
    else:
        fl=(2*a_l/(gamma-1))*(((p_star/pl))**((gamma-1)/(2*gamma))-1)
    return  fl        
   
def f_r(p_star,gamma,pr,rho_r,a_r):  
    if p_star>pr:
        fr=(p_star-pr)*(((2/((gamma+1)*rho_r))/(p_star+(pr*(gamma-1)/(gamma+1)))))**(0.5)
    else:
        fr=(2*a_r/(gamma-1))*(((p_star/pr))**((gamma-1)/(2*gamma))-1)
    return  fr       

def f_dash_l(p_star,gamma,pl,rho_l,a_l):
    if p_star>pl :
        fd_l=(((2/((gamma+1)*rho_l*(p_star+(pl*(gamma-1)/(gamma+1)))))))**(0.5)*(1-(0.5*(p_star-pl))/(p_star+(pl*(gamma-1)/(gamma+1))))
    else:
        fd_l=(1/(rho_l*a_l))*(((p_star/pl))**((-gamma-1)/(2*gamma)))
    return fd_l    

def f_dash_r(p_star,gamma,pr,rho_r,a_r):
    if p_star>pr :
        fd_r=(((2/((gamma+1)*rho_r*(p_star+(pr*(gamma-1)/(gamma+1))))))**(0.5))*(1-(0.5*(p_star-pl))/(p_star+(pr*(gamma-1)/(gamma+1))))
    else:
        fd_r=(1/(rho_r*a_r))*(((p_star/pr))**((-gamma-1)/(2*gamma)))        
    return fd_r

def Newton_rapson(pl,pr):           # Newton Rapson function
    if(pl == pr):                   # pl and pr are pressures on left and right side respectively.
        p_star= 10**(-6)
    else:
        p_star = 0.5*(pl+pr)
    d=1.0                        # d is a difference between new and old values of p_star(initial it is given as 1)
    while d>=10**(-6):
        func=f_l(p_star,gamma,pl,rho_l,a_l)+f_r(p_star,gamma,pr,rho_r,a_r)+ur-ul
        func_der=f_dash_l(p_star,gamma,pl,rho_l,a_l) + f_dash_r(p_star,gamma,pr,rho_r,a_r)
        p_starnew=p_star-(func/func_der)
        d=abs(p_starnew-p_star)
        p_star=p_starnew
    return p_starnew

l=[[1.0,0.4,1000.0,0.01,460.894],[0.0,-2.0,0.0,0.0,19.5975],[1.0,1.0,1.0,1.0,5.99924]]
m=[[0.1,0.4,0.01,100.0,46.0950],[0.0,2.0,0.0,0.0,-6.19633],[0.125,1.0,1.0,1.0,5.99242]]

for i in range (5):
    print(i+1," :")
    pl,ul,rho_l=l[0][i],l[1][i],l[2][i]
    pr,ur,rho_r=m[0][i],m[1][i],m[2][i]
    gamma=1.4
    a_l=(gamma*pl/rho_l)**(0.5)
    a_r=(gamma*pr/rho_r)**(0.5)
    
    p_st=Newton_rapson(pl,pr)   # value of p_star
    print("p_star value is ",p_st)
    u_star=0.5*(ul+ur)+ 0.5*(f_r(p_st,gamma,pr,rho_r,a_r)- f_l(p_st,gamma,pl,rho_l,a_l))
    print("u_star value is",u_star)
    if p_st>pl:
        rho_l_star=(rho_l*((p_st/pl)+((gamma-1)/(gamma+1)))/(((p_st*(gamma-1))/((gamma+1)*pl))+1))
        a_l_star=(gamma*p_st/rho_l_star)**(0.5)
        s_l=((ul-a_l*(((gamma+1)*p_st/(2*gamma*pl))+((gamma-1)/(2*gamma)))**(0.5)))
        print("rho_l_star is", rho_l_star)
        print("left side shock s_l is",s_l)
    else :
        rho_l_star=rho_l*((p_st/pl)**(1/gamma))
        a_l_star=(gamma*p_st/rho_l_star)**(0.5)
        R_LH=ul-a_l
        R_LT=u_star-a_l_star
        print("rho_l_star is", rho_l_star)
        print("left side rarefaction values are R_LH=",R_LH,"and R_LT=",R_LT)
    if p_st>pr:  
        rho_r_star=(rho_r*((p_st/pr)+((gamma-1)/(gamma+1)))/(((p_st*(gamma-1))/((gamma+1)*pr))+1))
        a_r_star=(gamma*p_st/rho_r_star)**(0.5)
        s_r=((ur+a_r*(((gamma+1)*p_st/(2*gamma*pr))+((gamma-1)/(2*gamma)))**(0.5)))
        print("rho_r_star is", rho_r_star)
        print("Right side shock s_l is",s_r)
        print("\n\n")
    else:
        rho_r_star=rho_r*((p_st/pr)**(1/gamma))
        a_r_star=(gamma*p_st/rho_r_star)**(0.5)
        R_RH=ur+a_r
        R_RT=u_star+a_r_star
        print("rho_r_star is", rho_r_star)
        print("Right side rarefaction values are R_LH=",R_RH,"and R_LT=",R_RT)
        print("\n\n")

sys.stdout = orig_stdout
fp.close()





