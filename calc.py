from clifford.sta import * 
from scipy.constants import e 


thres= 1e-10
is_close= lambda x,y:abs(x-y)<thres
is_close_ish= lambda x,y:abs(x-y)<thres*1e4



def da(f,a,tau=1e-9):
    # the a-derivative of function 'f()'   in direction of 'a'  chap2 (eq 1.5)
    return lambda x:(f(x+tau*a) - f(x))/tau       

def d(f,grade=1,tau=1e-9):
    # the 'grade'-derivative of 'f', returns df(), a function
    return lambda x: sum([ 1/a*da(f,a,tau=tau)(x) for a in D.blades_of_grade(grade)]) 


def dfg(f,g,grade=1,tau=1e-9):
    # the product rule  
    return lambda x: sum([ 1/a* (da(f,a,tau=tau)(x)*g(x) +f(x)*da(g,a,tau=tau)(x)) 
                          for a in D.blades_of_grade(grade)]) 

def dfofg(f,g,tau=1e-9):
    # chain rule
    pass
    #f(da(g,a,tau=tau)
    #return lambda x: sum([ 1/a* da(f(da(g,a,tau=tau)),tau=tau)(x)  for a in [d0,d1,d2,d3]]) 
#dfofg(lambda x: e**(P|x),lambda x: x)(x)

def test_d1():
    x = D.randomV()    
    B = D.randomMV()(2)
    f = lambda x: x|B# skew metric
    df = d(f,tau=1e-5)(x).clean(.001) # why does this get worse as tau->0?
    assert(is_close(df,2*B))

def test_d3():
    x = D.randomV()    
    T = D.randomMV()(3)
    f = lambda x: x|T# ?
    df = d(f,tau=1e-5)(x).clean(.001) # why does this get worse as tau->0?
    assert(is_close(df,3*T))

def test_d2():
    x = D.randomV()    
    f = lambda x: x**2 # simple scalar field 
    df = d(f,tau=1e-8)(x) # why does this get worse as tau->0?
    assert(is_close_ish(df, 2*x))

def test_d4():
    x = D.randomV()
    B = D.randomMV()(2)
    t = D.randomMV()(0) 
    
    f = lambda t: e**(B*t)*x*~(e**(B*t))
    assert(is_close_ish(d(f,0)(0).clean(1e-3),(2*B|x)))

test_d1()
test_d2()
test_d3()
test_d4()