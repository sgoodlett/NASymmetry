from symtext import *
from symel_generators import *

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
                σds = generate_σd(n)
            else:
                n = pg.n
                σds = []
            σvs = generate_σv(pg.n)
            symels = symels + cns + σvs + σds
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
                σds = generate_σd(n)
                c2ps = generate_C2p(pg.n)
                c2pps = generate_C2pp(pg.n)
                c2s = c2ps + c2pps
            else:
                n = pg.n
                σds = []
                c2s = generate_C2p(pg.n)
            cns = generate_Cn(pg.n)
            sns = generate_Sn(pg.n)
            σvs = generate_σv(pg.n)
            #c2s = generate_C2(pg.n)
            symels = symels + cns + c2s + sns + σvs + σds
        elif pg.subfamily == "d":
            if pg.n % 2 == 0:
                c2ps = generate_C2p(pg.n)
                c2pps = generate_C2pp(pg.n)
                c2s = c2ps + c2pps
            else:
                c2s = generate_C2p(pg.n)
                symels.append(Symel("i", inversion_matrix()))
            cns = generate_Cn(pg.n)
            sns = generate_Sn(pg.n * 2, True)
            σds = generate_σd(pg.n)
            symels = symels + cns + sns + c2s + σds
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
            sns = generate_Sn(pg.n, True)
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

if __name__ == "__main__":
    pg_to_symels("C4")