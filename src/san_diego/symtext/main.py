from .symtext import *
from .symel_generators import *
from .character_table_generators import *
from san_diego.symtools import normalize
import psi4
from san_diego.flowchart import find_point_group
from .multiplication_table import *

def pg_to_symels(PG):
    pg = PointGroup.from_string(PG)
    symels = [Symel("E", np.asarray([[1,0,0],[0,1,0],[0,0,1]]))]
    sigma_h = np.asarray([[1,0,0],[0,1,0],[0,0,-1]])
    if pg.family == "C":
        if pg.subfamily == "h":
            symels.append(Symel("sigma_h", sigma_h))
            if pg.n % 2 == 0:
                symels.append(Symel("i", inversion_matrix()))
            cns = generate_Cn(pg.n)
            sns = generate_Sn(pg.n)
            symels = symels + cns + sns
        elif pg.subfamily == "v":
            cns = generate_Cn(pg.n)
            if pg.n % 2 == 0:
                n = pg.n >> 1
                sigma_ds = generate_sigma_d(n)
            else:
                n = pg.n
                sigma_ds = []
            sigma_vs = generate_sigma_v(pg.n)
            symels = symels + cns + sigma_vs + sigma_ds
        elif pg.subfamily == "s":
            symels.append(Symel("sigma_h", sigma_h))
        elif pg.subfamily == "i":
            symels.append(Symel("i", inversion_matrix()))
        elif pg.subfamily is None:
            cns = generate_Cn(pg.n)
            symels = symels + cns
        else:
            raise Exception("Unidentified Point Group!")
    elif pg.family == "D":
        if pg.subfamily == "h":
            symels.append(Symel("sigma_h", sigma_h))
            if pg.n % 2 == 0:
                symels.append(Symel("i", inversion_matrix()))
                n = pg.n >> 1
                sigma_ds = generate_sigma_d(n)
                c2ps = generate_C2p(pg.n)
                c2pps = generate_C2pp(pg.n)
                c2s = c2ps + c2pps
            else:
                n = pg.n
                sigma_ds = []
                c2s = generate_C2p(pg.n)
            cns = generate_Cn(pg.n)
            sns = generate_Sn(pg.n)
            sigma_vs = generate_sigma_v(pg.n)
            #c2s = generate_C2(pg.n)
            symels = symels + cns + c2s + sns + sigma_vs + sigma_ds
        elif pg.subfamily == "d":
            if pg.n % 2 == 0:
                c2ps = generate_C2p(pg.n)
                c2pps = generate_C2pp(pg.n)
                c2s = c2ps + c2pps
            else:
                c2s = generate_C2p(pg.n)
                symels.append(Symel("i", inversion_matrix()))
            cns = generate_Cn(pg.n)
            sns = generate_Sn(pg.n * 2, S2n=True)
            sigma_ds = generate_sigma_d(pg.n)
            symels = symels + cns + sns + c2s + sigma_ds
        elif pg.subfamily is None:
            cns = generate_Cn(pg.n)
            if pg.n % 2 == 0:
                c2ps = generate_C2p(pg.n)
                c2pps = generate_C2pp(pg.n)
                c2s = c2ps + c2pps
            else:
                c2s = generate_C2p(pg.n)
            symels = symels + cns + c2s
        else:
            raise Exception("Oh shit, the trout population")
    elif pg.family == "S":
        if pg.subfamily is None and (pg.n % 2 == 0):
            n = pg.n >> 1
            if n % 2 != 0:
                symels.append(Symel("i", inversion_matrix()))
            cns = generate_Cn(n)
            sns = generate_Sn(pg.n, S2n=True)
            symels = symels + cns + sns
        else:
            raise Exception("Oh shit, the trout population")
    else:
        if pg.family == "T":
            if pg.subfamily == "h":
                Ths = generate_Th()
                symels = symels + Ths
            elif pg.subfamily == "d":
                Tds = generate_Td()
                symels = symels + Tds
            else:
                Ts = generate_T()
                symels = symels + Ts
        elif pg.family == "O":
            if pg.subfamily == "h":
                Ohs = generate_Oh()
                symels = symels + Ohs
            else:
                Os = generate_O()
                symels = symels + Os
        elif pg.family == "I":
            if pg.subfamily == "h":
                Ihs = generate_Ih()
                symels = symels + Ihs
            else:
                Is = generate_I()
                symels = symels + Is
        else:
            raise Exception("Unidentified Point Group!")
    return symels

def pg_to_chartab(PG):
    pg = PointGroup.from_string(PG)
    irreps = []
    if pg.family == "C":
        if pg.subfamily == "s":
            irreps = ["A'","A''"]
            classes = ["E", "sigma_h"]
            chars = np.array([[1, 1], [1 -1]])
        elif pg.subfamily == "i":
            irreps = ["Ag","Au"]
            classes = ["E", "i"]
            chars = np.array([[1, 1], [1, -1]])
        elif pg.subfamily == "v":
            irreps, classes, chars = Cnv_irr(pg.n)
        elif pg.subfamily == "h":
            irreps, classes, chars = Cnh_irr(pg.n)
        else:
            irreps, classes, chars = Cn_irrmat(pg.n)
    elif pg.family == "D":
        if pg.subfamily == "d":
            irreps, classes, chars = Dnd_irr(pg.n)
        elif pg.subfamily == "h":
            irreps, classes, chars = Dnh_irr(pg.n)
        else:
            irreps, classes, chars = Dn_irr(pg.n)
    elif pg.family == "S":
        irreps, classes, chars = Sn_irr(pg.n)
    else:
        cp3 = np.cos(np.pi/3)
        pr5 = 0.5*(1.0+np.sqrt(5.0))
        mr5 = 0.5*(1.0-np.sqrt(5.0))
        if pg.family == "T":
            if pg.subfamily == "h":
                irreps, classes, chars = (["Ag","Au","Eg","Eu","Tg","Tu"],
                 ["E", "4C_3", "4C_3^2", "3C_2", "i", "S_6", "S_6^5", "3sigma_h"],
                 [[1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0],
                  [1.0,  1.0,  1.0,  1.0, -1.0, -1.0, -1.0, -1.0],
                  [2.0,  cp3,  cp3,  2.0,  2.0,  cp3,  cp3,  1.0],
                  [2.0,  cp3,  cp3,  2.0, -2.0, -cp3, -cp3, -1.0],
                  [3.0,  0.0,  0.0, -1.0,  1.0,  0.0,  0.0, -1.0],
                  [3.0,  0.0,  0.0, -1.0, -1.0,  0.0,  0.0,  1.0]])
            elif pg.subfamily == "d":
                irreps, classes, chars = (["A1","A2","E","T1","T2"],
                 ["E", "8C_3", "3C_2", "6S_4", "6sigma_d"],
                 [[1.0,  1.0,  1.0,  1.0,  1.0],
                  [1.0,  1.0,  1.0, -1.0, -1.0],
                  [2.0, -1.0,  2.0,  0.0,  0.0],
                  [3.0,  1.0, -1.0,  1.0, -1.0],
                  [3.0, -1.0, -1.0, -1.0,  1.0]])
            else:
                irreps, classes, chars = (["A","E","T"],
                 ["E", "4C_3", "4C_3^2", "3C_2"],
                 [[1.0,  1.0,  1.0,  1.0],
                  [2.0,  cp3,  cp3,  2.0],
                  [3.0,  0.0,  0.0, -1.0]])
        elif pg.family == "O":
            if pg.subfamily == "h":
                irreps, classes, chars = (["A1g","A2g","Eg","T1g","T2g","A1u","A2u","Eu","T1u","T2u"],
                 ["E", "8C_3", "6C_2", "6C_4", "3C_2", "i", "6S_4", "8S_6", "3sigma_h", "6sigma_d"],
                 [[1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0],
                  [1.0,  1.0, -1.0, -1.0,  1.0,  1.0, -1.0,  1.0,  1.0, -1.0],
                  [2.0, -1.0,  0.0,  0.0,  2.0,  2.0,  0.0, -1.0,  2.0,  0.0],
                  [3.0,  0.0, -1.0,  1.0, -1.0,  3.0,  1.0,  0.0, -1.0, -1.0],
                  [3.0,  0.0,  1.0, -1.0, -1.0,  3.0, -1.0,  0.0, -1.0,  1.0],
                  [1.0,  1.0,  1.0,  1.0,  1.0, -1.0, -1.0, -1.0, -1.0, -1.0],
                  [1.0,  1.0, -1.0, -1.0,  1.0, -1.0,  1.0, -1.0, -1.0,  1.0],
                  [2.0, -1.0,  0.0,  0.0,  2.0, -2.0,  0.0,  1.0, -2.0,  0.0],
                  [3.0,  0.0, -1.0,  1.0, -1.0, -3.0, -1.0,  0.0,  1.0,  1.0],
                  [3.0,  0.0,  1.0, -1.0, -1.0, -3.0,  1.0,  0.0,  1.0, -1.0]])
            else:
                irreps, classes, chars = (["A1","A2","E","T1","T2"],
                 ["E", "6C_4", "3C_2", "8C_3", "6C_2"],
                 [[1.0,  1.0,  1.0,  1.0,  1.0],
                  [1.0, -1.0,  1.0,  1.0, -1.0],
                  [2.0,  0.0,  2.0, -1.0,  0.0],
                  [3.0,  1.0, -1.0,  0.0, -1.0],
                  [3.0, -1.0, -1.0,  0.0,  1.0]])
        elif pg.family == "I":
            if pg.subfamily == "h":
                irreps, classes, chars = (["Ag","T1g","T2g","Gg","Hg","Au","T1u","T2u","Gu","Hu"],
                 ["E", "12C_5", "12C_5^2", "20C_3", "15C_2", "i", "12S_10", "12S_10^3", "20S_6", "15sigma_"],
                 [[1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0],
                  [3.0,  pr5,  mr5,  0.0, -1.0,  3.0,  mr5,  pr5,  0.0, -1.0],
                  [3.0,  mr5,  pr5,  0.0, -1.0,  5.0,  0.0,  0.0, -1.0,  1.0],
                  [1.0,  1.0,  1.0,  1.0,  1.0, -1.0, -1.0, -1.0, -1.0, -1.0],
                  [3.0,  pr5,  mr5,  0.0, -1.0, -3.0, -mr5, -pr5,  0.0,  1.0],
                  [3.0,  mr5,  pr5,  0.0, -1.0, -3.0, -pr5, -mr5,  0.0,  1.0],
                  [4.0, -1.0, -1.0,  1.0,  0.0, -4.0,  1.0,  1.0, -1.0,  0.0],
                  [5.0,  0.0,  0.0, -1.0,  1.0, -5.0,  0.0,  0.0,  1.0, -1.0]])
            else:
                irreps, classes, chars = (["A","T1","T2","G","H"],
                 ["E", "12C_5", "12C_5^2", "20C_3", "15C_2"],
                 [[1.0,  1.0,  1.0,  1.0,  1.0],
                  [3.0,  pr5,  mr5,  0.0, -1.0],
                  [3.0,  mr5,  pr5,  0.0, -1.0],
                  [4.0, -1.0, -1.0,  1.0,  0.0],
                  [5.0,  0.0,  0.0, -1.0,  1.0]])
        else:
            raise Exception("Unrecognized Point Group")
    class_orders = grab_class_orders(classes)
    irr_dims = {}
    for (irr_idx,irrep) in enumerate(irreps):
        irr_dims[irrep] = chars[irr_idx, 1]
    return CharTable(PG, irreps, classes, class_orders, chars, irr_dims)

def grab_class_orders(classes):
    ncls = len(classes)
    class_orders = np.zeros(ncls)
    for i in range(ncls): # = 1:ncls
        class_orders[i] = grab_order(classes[i])
    return class_orders

def grab_order(class_str):
    regex = r"^(\d+)"
    m = re.match(regex, class_str)
    if m is not None:
        return int(m.groups()[0])
    else:
        return 1

def generate_symel_to_class_map(symels, ctab):
    return [1,2,3]
    pg = PointGroup.from_string(ctab.name)
    if pg.n is not None:
        ns = pg.n>>1 # pg.n floor divided by 2
    ncls = len(ctab.classes)
    nsymel = len(symels)
    class_map = np.zeros(nsymel)
    class_map[0] = 1 # E is always first
    if pg.family == "C":
        if pg.subfamily == "s" or pg.subfamily == "i":
            class_map[2] = 2
        elif pg.subfamily == "h":
            if pg.n % 2 == 0:
                class_map[4:pg.n+2] = [i for i in range(2,pg.n+1)] # [2:pg.n] # C_n
                class_map[3] = pg.n+1 # i
                class_map[2] = pg.n+ns+1 # σh
                for i in range(pg.n+3,2*pg.n+1):# = pg.n+3:2*pg.n # S_n
                    if i > 3*ns+1:
                        class_map[i] = i-ns
                    else:
                        class_map[i] = i+ns-1
            else:
                for i in range(2,pg.n): # = 2:pg.n+1 # C_n
                    class_map[i] = i-1
                class_map[2] = pg.n+1 # σh
                for i in range(pg.n+2,2*pg.n+1): # = pg.n+2:2*pg.n # S_n
                    class_map[i] = i
        elif pg.subfamily == "v":
            # The last class is σv (and then σd if n is even), and the last symels are also these!
            cn_class_map(class_map, pg.n, 0, 0)
            if pg.n % 2 == 0:
                class_map[-1-pg.n+1:-1-ns] = ncls-1
                class_map[-1-ns+1:-1] = ncls
            else:
                class_map[-1-pg.n+1:-1] = ncls
        else:
            class_map[2:-1] = [i for i in range(2,nsymel+1)] # 2:nsymel
    elif pg.family == "S":
        if pg.n % 4 == 0:
            for i in range(2,pg.n+1): # = 2:pg.n
                if i <= ns:
                    class_map[i] = 2*i-1
                else:
                    class_map[i] = 2*(i-ns)
        else:
            class_map[2] = ns+1 # i
            class_map[3:ns+1] = [i for i in range(2,ns+1)] # 2:ns # C_n
            for i in range(ns+2,pg.n+1): # = ns+2:pg.n # S_n
                if i > ns+(pg.n>>2)+1:
                    class_map[i] = i-(ns-1)>>1
                else:
                    class_map[i] = i + (ns-1)>>1
    elif pg.family == "D":
        if pg.subfamily == "h":
            if pg.n % 2 == 0:
                class_map[2] = ncls-2 # σh
                class_map[3] = (ncls>>1)+1 # i
                cn_class_map(class_map, pg.n, 2, 0) # Cn
                class_map[pg.n+3:3*ns+2] = ns+2 # C2'
                class_map[3*ns+3:2*pg.n+2] = ns+3 # C2''
                for i in range(2*pg.n_3,3*pg.n+1): # = 2*pg.n+3:3*pg.n+1 # Sn
                    if i > 3*pg.n-ns+1:
                        class_map[i] = i-2*pg.n+3
                    else:
                        class_map[i] = 3*pg.n+6-i
                # The result of C2'×i changes depending if pg.n ≡ 0 (mod 4)
                if pg.n % 4 == 0:
                    class_map[-1-pg.n+1:-1-ns] = ncls-1 # σv
                    class_map[-1-ns+1:-1] = ncls # σd
                else:
                    class_map[-1-pg.n+1:-1-ns] = ncls # σv
                    class_map[-1-ns+1:-1] = ncls-1 # σd
            else:
                class_map[2] = (ncls>>1)+1
                cn_class_map(class_map, pg.n, 1, 0)
                class_map[pg.n+2:2*pg.n+1] = ns+2
                cn_class_map(class_map, pg.n, 2*pg.n, ns+2)
                class_map[-1-pg.n+1:-1] = ncls
        elif pg.subfamily == "d":
            if pg.n % 2 == 0:
                cn_class_map(class_map, pg.n, 0, 0) # Cn
                class_map[2:pg.n] = [2*i-1 for i in range(2,pg.n+1)] # 2*class_map[2:pg.n].-1 # Reposition Cn
                cn_class_map(class_map, pg.n+1, pg.n-1, 0) # Sn
                class_map[pg.n+1:2*pg.n] = [2*i-1 for i in range(pg.n+1,2*pg.n+1)] # 2*(class_map[pg.n+1:2*pg.n].-1) # Reposition Sn
                class_map[-1-2*pg.n+1:-1-pg.n] = ncls-1 # C2'
                class_map[-1-pg.n+1:-1] = ncls # σd
            else:
                class_map[2] = (ncls>>1)+1 # i
                cn_class_map(class_map, pg.n, 1, 0) # Cn
                for i in range(pg.n+2,2*pg.n): # = pg.n+2:2*pg.n # Sn
                    if i > pg.n+1+ns:
                        class_map[i] = i+2-pg.n
                    else:
                        class_map[i] = 2*pg.n+4-i
                class_map[-1-2*pg.n+1:-1-pg.n] = ns+2
                class_map[-1-pg.n+1:-1] = ncls # σd
        else:
            cn_class_map(class_map, pg.n, 0, 0) # Cn
            if pg.n % 2 == 0:
                class_map[-1-pg.n+1:-1-ns] = ncls-1 # Cn'
                class_map[-1-ns+1:-1] = ncls # Cn''
            else:
                class_map[-1-pg.n+1:-1] = ncls # Cn
    else:
        if pg.family == "T":
            if pg.subfamily == "h":
                class_map = np.array([1,2,3,2,3,2,3,2,3,4,4,4,5,6,7,6,7,6,7,6,7,8,8,8])
            elif pg.subfamily == "d":
                class_map = np.array([1,2,2,2,2,2,2,2,2,3,3,3,5,5,5,5,5,5,4,4,4,4,4,4])
            else:
                class_map = np.array([1,2,3,2,3,2,3,2,3,4,4,4])
        elif pg.family == "O":
            if pg.subfamily == "h":
                class_map = np.array([1,4,5,4,4,5,4,4,5,4,2,2,2,2,2,2,2,2,3,3,3,3,3,3,6,7,9,7,7,9,7,7,9,7,8,8,8,8,8,8,8,8,10,10,10,10,10,10])
            else:
                class_map = np.array([1,2,3,2,2,3,2,2,3,2,4,4,4,4,4,4,4,4,5,5,5,5,5,5])
        elif pg.family == "I":
            if pg.subfamily == "h":
                class_map = np.array([1,2,3,3,2,2,3,3,2,2,3,3,2,2,3,3,2,2,3,3,2,2,3,3,2,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
                             6,7,8,8,7,7,8,8,7,7,8,8,7,7,8,8,7,7,8,8,7,7,8,8,7,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10])
            else:
                class_map = np.array([1,2,3,3,2,2,3,3,2,2,3,3,2,2,3,3,2,2,3,3,2,2,3,3,2,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5])
        else:
            raise Exception(f"Unrecognized point group: {ctab.name}")
    return class_map

def cn_class_map(class_map, n, idx_offset, cls_offset):
    for i in range(2,n+1): # = 2:n
        if i > (n>>1)+1:
            class_map[i+idx_offset] = n-i+2+cls_offset
        else:
            class_map[i+idx_offset] = i+cls_offset

def rotate_mol_to_symels(mol, paxis, saxis):
    phi, theta, chi = get_euler_angles(paxis, saxis)
    dc = dc_mat(phi, theta, chi)
    print("Beebus")
    print(dc)
    new_mol = mol.transform(dc)
    return new_mol

def rotate_symels_to_mol(symels, paxis, saxis):
    phi, theta, chi = get_euler_angles(paxis, saxis)
    dc = dc_mat(phi, theta, chi)
    dct = deepcopy(dc)
    dct = dct.transpose()
    new_symels = []
    for s in symels:
        new_symels.append(Symel(s.symbol, dct*s.rrep*dc))
    return new_symels

def get_euler_angles(paxis, saxis):
    x = np.array([1.0,0.0,0.0])
    y = np.array([0.0,1.0,0.0])
    z = np.array([0.0,0.0,1.0])
    ynew = np.cross(paxis,saxis)
    zxp = normalize(np.cross(z,paxis))
    if zxp is None:
        phi = 0.0
    else:
        xproj = np.dot(zxp,x)
        if xproj <= 0:
            phi = np.arccos(np.dot(y,zxp))
        else:
            phi = 2*np.pi-np.arccos(np.dot(y,zxp))
    r_phi = rotation_matrix(z,phi)
    yN = np.dot(r_phi,y)
    xN = np.dot(r_phi,x)
    theta = np.arccos(np.dot(z,paxis))
    print(yN)
    r_theta = rotation_matrix(yN,theta)
    x2N = np.dot(r_theta,xN)
    Z = np.dot(r_theta,z)
    chi = np.arccos(np.dot(x2N,saxis))
    r_chi = rotation_matrix(Z,chi)
    X = np.dot(r_chi,x2N)
    Y = np.dot(r_chi,yN)
    return phi, theta, chi

def dc_mat(phi, theta, chi):
    sp = np.sin(phi)
    cp = np.cos(phi)
    st = np.sin(theta)
    ct = np.cos(theta)
    sc = np.sin(chi)
    cc = np.cos(chi)
    direction_cosine = np.array([[cp*ct*cc-sp*sc, sp*ct*cc+cp*sc, -cc*st], [-cp*ct*sc-sp*cc, -sp*ct*sc+cp*cc, sc*st], [st*cp, st*sp, ct]])
    return direction_cosine

def get_atom_mapping(mol, symels):
    # symels after transformation
    amap = np.zeros((len(mol), len(symels)))
    for (a, atom) in enumerate(mol):
        for (s, symel) in enumerate(symels):
            w = where_you_go(mol, atom, symel)
            if w is not None:
                amap[a,s] = w
            else:
                raise Exception("Atom $(atom) not mapped to another atom under symel $(symel)")
    return amap

def where_you_go(mol, atom, symel):
    ratom = np.dot(mol.coords[atom,:], symel.rrep.transpose())
    length = np.shape(mol)[0]
    for i in range(length):
        if np.isclose(mol.coords[i,:], ratom, atol=tol).all():
            return i
    return None

def symtext_from_file(fn):
    with open(fn, "r") as lfn:
        strang = lfn.read()
    mool = psi4.core.Molecule.from_string(strang)
    beebus = mool.to_schema("psi4")
    mol = Molecule.from_schema(beebus)
    return symtext_from_mol(mol)

def symtext_from_mol(mol):
    print("Weebus")
    mol.translate(mol.find_com())
    print(mol)
    pg, (paxis, saxis) = find_point_group(mol)
    symels = pg_to_symels(pg)
    mol = rotate_mol_to_symels(mol, paxis, saxis)
    ctab = pg_to_chartab(pg)
    class_map = generate_symel_to_class_map(symels, ctab)
    atom_map = get_atom_mapping(mol, symels)
    mtable = build_mult_table(symels)
    return mol, Symtext(pg, symels, ctab, class_map, atom_map, mtable)

def irrep_sort_idx(irrep_str):
    bunyon = 0
    # g and ' always go first
    gchk = r"g"
    ppchk = r"''"
    gm = gchk in irrep_str
    pm = ppchk in irrep_str
    if gm:
        bunyon += 0
    elif pm:
        bunyon += 10000
    else:
        bunyon += 1000
    bub = irrep_str[1] # the letter
    numbro = r"(\d+)"
    mn = re.match(numbro, irrep_str)
    if mn is not None:
        bunyon += int(mn.groups()[0])
    if bub == 'A':
        bunyon += 0
    elif bub == 'B':
        bunyon += 100
    elif bub == 'E':
        bunyon += 200
    elif bub == 'T':
        bunyon += 300
    elif bub == 'G':
        bunyon += 400
    elif bub == 'H':
        bunyon += 500
    else:
        raise Exception("Oh shit the trout population!")
    return bunyon

if __name__ == "__main__":
    p = symtext_from_file("../test/xyz/ammonia.xyz")
    print(p)