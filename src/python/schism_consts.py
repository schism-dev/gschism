
natrm = 12      # Number of _available_ tracer models at the moment (including T,S)
mntracers = 30  # Maximum number of tracers, used only for dimensioning btrack arrays. Must >=ntracers

blockSize = 256

DEBUG = False

varCodes = "x,y,z,u,v,w,U,V,nx,ny,hhat,ghatx,ghaty,ssign,btu,btv,bccx,bccy,uth,vth,ubc,vbc," + \
           "island,ifltype,elbc,isdry,eta,area,qel,depth,radius,xlon,ylat,kb,n1,n2," + \
           "htot,gradEta1x,gradEta1y,gradEta2x,gradEta2y,gradETPx,gradETPy,gradPrx,gradPry,idx"
