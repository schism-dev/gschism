
natrm = 12      # Number of _available_ tracer models at the moment (including T,S)
mntracers = 30  # Maximum number of tracers, used only for dimensioning btrack arrays. Must >=ntracers

blockSize = 256

DEBUG = False

varInfoDict =  {'x':       {'is3D': False, 'type': 'Vector2'},
                'z':       {'is3D': True,  'type': 'double'},
                'u':       {'is3D': True,  'type': 'Vector2'},
                'w':       {'is3D': True,  'type': 'double'},
                'U':       {'is3D': False, 'type': 'Vector2'},
                'nx':      {'is3D': False, 'type': 'int'},
                'ny':      {'is3D': False, 'type': 'int'},
                'hhat':    {'is3D': False, 'type': 'double'},
                'ghat':   {'is3D': False, 'type':  'Vector2'},
                'ssign':   {'is3D': False, 'type': 'double'},
                'btu':     {'is3D': True,  'type': 'Vector2'},
                'bcc':     {'is3D': True,  'type': 'Vector2'},
                'uth':     {'is3D': True,  'type': 'Vector2'},
                'ubc':     {'is3D': True,  'type': 'Vector2'},
                'island':  {'is3D': False, 'type': 'int'},
                'ifltype': {'is3D': False, 'type': 'int'},
                'elbc':    {'is3D': False, 'type': 'double'},
                'isdry':   {'is3D': False, 'type': 'bool'},
                'eta':     {'is3D': False, 'type': 'double'},
                'area':    {'is3D': False, 'type': 'double'},
                'qel':     {'is3D': False, 'type': 'double'},
                'depth':   {'is3D': False, 'type': 'double'},
                'radius':  {'is3D': False, 'type': 'double'},
                'len':     {'is3D': False, 'type': 'double'},
                'lonlat':  {'is3D': False, 'type': 'Vector2'},
                'kb':      {'is3D': False, 'type': 'int'},
                'n1':      {'is3D': False, 'type': 'pointer'},
                'n2':      {'is3D': False, 'type': 'pointer'},
                'idx':     {'is3D': False, 'type': 'int'},
                'htot':    {'is3D': False, 'type': 'double'},
                'gradPr':  {'is3D': False, 'type': 'Vector2'},
                'gradEta1':{'is3D': False, 'type': 'Vector2'},
                'gradEta2':{'is3D': False, 'type': 'Vector2'},
                'gradETP': {'is3D': False, 'type': 'Vector2'},
                'FaceIndices': {'is3D': False, 'type': 'FaceIndices'},
                }
# varCodes = "x,y,z,u,v,w,U,V,nx,ny,hhat,ghatx,ghaty,ssign,btu,btv,bccx,bccy,uth,vth,ubc,vbc," + \
#            "island,ifltype,elbc,isdry,eta,area,qel,depth,radius,xlon,ylat,kb,n1,n2," + \
#            "htot,gradEta1x,gradEta1y,gradEta2x,gradEta2y,gradETPx,gradETPy,gradPrx,gradPry,idx"
