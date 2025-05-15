import cupy as cp
from os.path import join

from python.schism_consts import *
class TidalData:
    def __init__(self):
        # Initialize all possible attributes as None
        self.ntip = None
        self.tip_dp = None
        self.tp_name = None
        self.tamp = None
        self.tnf = None
        self.tfreq = None
        self.jspc = None
        self.tear = None
        self.nbfr = None
        self.amig = None
        self.ff = None
        self.face = None
        self.nope = None
        self.iettype = None
        self.ifltype = None
        self.eth = None
        self.el_amp = None
        self.el_phase = None
        self.qthcon = None
        self.u_amp = None
        self.u_phase = None
        self.v_amp = None
        self.v_phase = None
        self.vobc1 = None
        self.vobc2 = None
        self.disch_coef = None
        self.emo = None
        self.efa = None
        self.umo = None
        self.ufa = None
        self.vmo = None
        self.vfa = None
        self.trobc = None
        self.itrtype = None
        self.trth = None

float64 = cp.float64
int32   = cp.int32

def read_values(line, types):
    # Convert Fortran's D notation to E notation for float parsing
    parts = [p.replace('D', 'E') for p in line.split()]
    return tuple(t(x) for t, x in zip(types, parts))

def read_bctides(path: str, nope_global: int, mnond_global: int, ntrs: list, irange_tr: list) -> TidalData:
    """Read boundary condition tides data.
    
    Args:
        path: Path to bctides.in
        nope_global: Number of open boundaries
        mnond_global: Maximum number of nodes per boundary
        ntrs: List indicating if tracer model is active
        irange_tr: Range for each tracer model
    """
    filename = join(path, "bctides.in")
    tide = TidalData()
    
    with open(filename, 'r') as f:
        # Skip first line (start_time)
        line = f.readline()
        
        # Read tidal potential parameters
        line = f.readline()
        tide.ntip, tide.tip_dp = read_values(line, [int, float64])
        
        if tide.ntip > 0:
            tide.tp_name = []
            tide.tamp  = cp.zeros(tide.ntip, dtype=float64)
            tide.tnf   = cp.zeros(tide.ntip, dtype=float64)
            tide.tfreq = cp.zeros(tide.ntip, dtype=float64)
            tide.jspc  = cp.zeros(tide.ntip, dtype=int32)
            tide.tear  = cp.zeros(tide.ntip, dtype=float64)
            
            for i in range(tide.ntip):
                tide.tp_name.append(f.readline().strip())
                line = f.readline()
                jspc, tamp, tfreq, tnf, tear = read_values(line, [int] + [float64]*4)
                if jspc < 0 or jspc > 2:
                    raise ValueError(f'Illegal tidal species # {jspc}')
                tide.jspc[i] = jspc
                tide.tamp[i] = tamp
                tide.tfreq[i] = tfreq
                tide.tnf[i] = tnf
                tide.tear[i] = tear * cp.pi / 180.0
        
        # Read boundary forcing frequencies
        line = f.readline()
        tide.nbfr = read_values(line, [int])[0]
        
        if tide.nbfr > 0:
            tide.amig = cp.zeros(tide.nbfr, dtype=float64)
            tide.ff   = cp.zeros(tide.nbfr, dtype=float64)
            tide.face = cp.zeros(tide.nbfr, dtype=float64)
            tide.emo  = cp.zeros((nope_global, mnond_global, tide.nbfr), dtype=float64)
            tide.efa  = cp.zeros((nope_global, mnond_global, tide.nbfr), dtype=float64)
            tide.umo  = cp.zeros((nope_global, mnond_global, tide.nbfr), dtype=float64)
            tide.ufa  = cp.zeros((nope_global, mnond_global, tide.nbfr), dtype=float64)
            tide.vmo  = cp.zeros((nope_global, mnond_global, tide.nbfr), dtype=float64)
            tide.vfa  = cp.zeros((nope_global, mnond_global, tide.nbfr), dtype=float64)
            
            for i in range(tide.nbfr):
                line = f.readline()  # Skip tag
                line = f.readline()
                amig, ff, face = read_values(line, [float64]*3)
                tide.amig[i] = amig
                tide.ff[i] = ff
                tide.face[i] = face * cp.pi / 180.0
        
        # Read number of open boundaries
        line = f.readline()
        nope1 = read_values(line, [int])[0]
        if nope1 != nope_global:
            raise ValueError(f"Inconsistent # of open bnds {nope1} {nope_global}")
        
        # if ntmp != nond_global[k]:
        #     errmsg = f"Inconsistent # of nodes at open boundary {k} {ntmp} {nond_global[k]}"
        #     raise ValueError(errmsg)
        
        # Count number of active tracer models
        ntrmod = 2  # T,S
        for i in range(2, natrm):
            if ntrs[i] > 0:
                ntrmod += 1
        
        # Initialize counters
        nettype = 0   # total # of type I bnds
        nfltype = 0   # total # of type I flow bnds
        nettype2 = 0  # total # of type IV bnds
        nfltype2 = 0  # total # of type IV flow bnds
        nnode_et = 0  # nodes requiring elev2D.th.nc
        nnode_fl = 0  # nodes requiring uv3D.th.nc
        nnode_tr2 = cp.zeros(natrm, dtype=int32)  # nodes requiring tracer input
        
        # Initialize arrays
        tide.iettype = cp.zeros(nope_global, dtype=int32)
        tide.ifltype = cp.zeros(nope_global, dtype=int32)
        tide.eth = cp.zeros((nope_global,), dtype=float64)
        tide.qthcon = cp.zeros(nope_global, dtype=float64)
        tide.vobc1 = cp.zeros(nope_global, dtype=float64)
        tide.vobc2 = cp.zeros(nope_global, dtype=float64)
        tide.disch_coef = cp.zeros(4, dtype=float64)
        
        # Initialize tracer arrays
        tide.trobc = cp.zeros((natrm, nope_global), dtype=float64)
        tide.itrtype = cp.full((natrm, nope_global), -99, dtype=int32)
        tide.trth = cp.zeros((mntracers, 1, 1, nope_global), dtype=float64)
        ntrtype1 = cp.zeros(natrm, dtype=int32)

        # Process each open boundary
        for k in range(nope_global):
            line = f.readline()
            ntmp, iettype, ifltype, *nwild = read_values(line, [int]*(3 + ntrmod))
            tide.iettype[k] = iettype
            tide.ifltype[k] = ifltype

            # Process elevation boundary conditions
            if iettype == 1:
                nettype += 1
            elif iettype == 2:
                line = f.readline()
                tide.eth[k] = read_values(line, [float64])[0]
            elif iettype == 3:
                for i in range(tide.nbfr):
                    line = f.readline()  # Skip frequency name
                    for j in range(ntmp):
                        line = f.readline()
                        amp, phase = read_values(line, [float64]*2)
                        tide.emo[k, j, i] = amp
                        tide.efa[k, j, i] = phase * cp.pi / 180.0
            elif iettype in [4, 5]:
                nettype2 += 1
                nnode_et += ntmp
                
                if iettype == 5:
                    for i in range(tide.nbfr):
                        line = f.readline()  # Skip frequency name
                        for j in range(ntmp):
                            line = f.readline()
                            amp, phase = read_values(line, [float64]*2)
                            tide.emo[k, j, i] = amp
                            tide.efa[k, j, i] = phase * cp.pi / 180.0
            elif iettype != 0:
                raise ValueError('Invalid iettype')

            # Process flow boundary conditions
            if ifltype == 1:
                nfltype += 1
            elif ifltype == 2:
                line = f.readline()
                tide.qthcon[k] = read_values(line, [float64])[0]
            elif ifltype == 3:
                for i in range(tide.nbfr):
                    line = f.readline()  # Skip frequency name
                    for j in range(ntmp):
                        line = f.readline()
                        u_amp, u_phase, v_amp, v_phase = read_values(line, [float64]*4)
                        tide.umo[k, j, i] = u_amp
                        tide.ufa[k, j, i] = u_phase * cp.pi / 180.0
                        tide.vmo[k, j, i] = v_amp
                        tide.vfa[k, j, i] = v_phase * cp.pi / 180.0
            elif abs(ifltype) in [4, 5]:
                if ifltype < 0:
                    if iettype == 0:
                        raise ValueError(f'vel. obc needs elev. to be specified: {k}')
                    line = f.readline()
                    tide.vobc1[k], tide.vobc2[k] = read_values(line, [float64]*2)
                
                nfltype2 += 1
                nnode_fl += ntmp

                if abs(ifltype) == 5:
                    for i in range(tide.nbfr):
                        line = f.readline()  # Skip frequency name
                        for j in range(ntmp):
                            line = f.readline()
                            u_amp, u_phase, v_amp, v_phase = read_values(line, [float64]*4)
                            tide.umo[k, j, i] = u_amp
                            tide.vmo[k, j, i] = v_amp
                            tide.ufa[k, j, i] = cp.deg2rad(u_phase)                            
                            tide.vfa[k, j, i] = cp.deg2rad(v_phase)
            elif ifltype == -1:  # Flather 1
                if iettype != 0:
                    raise ValueError(f'Flather obc requires iettype=0: {k}')
                line = f.readline()  # Skip 'eta_mean'
                for _ in range(ntmp):
                    line = f.readline()  # Skip eta_m0
                line = f.readline()  # Skip 'vn_mean'
                for _ in range(ntmp):
                    line = f.readline()  # Skip uthnd
            elif ifltype == -2:  # Discharge relation
                if iettype != 0:
                    raise ValueError(f'Discharge obc requires iettype=0: {k}')
                line = f.readline()
                tide.disch_coef = read_values(line, [float64]*4)
            elif ifltype != 0:
                raise ValueError(f'Invalid ifltype: {ifltype}')

            # Process tracer boundary conditions
            icount = 0
            for i in range(natrm):
                if ntrs[i] > 0:
                    icount += 1
                    tide.itrtype[i,k] = nwild[icount - 1]
                    
                    if tide.itrtype[i,k] == 2:
                        line = f.readline()
                        trth_values = read_values(line, [float64]*(irange_tr[1][i] - irange_tr[0][i] + 1))
                        tide.trth[irange_tr[0][i]:irange_tr[1][i], 0, 0, k] = trth_values
                        
                        line = f.readline()
                        tide.trobc[i,k] = read_values(line, [float64])[0]
                    
                    elif tide.itrtype[i,k] in [1, 3, 4]:
                        line = f.readline()
                        tide.trobc[i,k] = read_values(line, [float64])[0]
                        
                        if tide.itrtype[i,k] == 1:
                            ntrtype1[i] += 1
                        elif tide.itrtype[i,k] == 4:
                            nnode_tr2[i] += ntmp
                    
                    elif tide.itrtype[i,k] != 0:
                        raise ValueError(f'Wrong itrtype: {k} {tide.itrtype[i,k]}')
                    
                    if not (0 <= tide.trobc[i,k] <= 1):
                        raise ValueError(f'Tr. obc nudging factor wrong: {tide.trobc[i,k]} {k}')
            
            if icount != ntrmod:
                raise ValueError('INIT:icount!=ntrmod')
        
        # Store counters in tide object
        tide.nettype   = nettype
        tide.nfltype   = nfltype
        tide.nettype2  = nettype2
        tide.nfltype2  = nfltype2
        tide.nnode_et  = nnode_et
        tide.nnode_fl  = nnode_fl
        tide.nnode_tr2 = nnode_tr2
        tide.ntrtype1  = ntrtype1
    
    return tide

