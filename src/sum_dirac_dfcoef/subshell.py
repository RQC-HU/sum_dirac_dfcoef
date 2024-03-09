from typing import ClassVar, List, Tuple


class Subshell:
    orbital_labels = "spdfghiklmnoqrtuvwxyz"
    gto_label_order: ClassVar[List[List[str]]] = []

    def __lmnval(self, idx: int, nfun: int, istep: List[int], mval: List[int], nval: List[int]) -> Tuple[List[int], List[int], List[int]]:
        # https://gitlab.com/dirac/dirac/-/blob/b10f505a6f00c29a062f5cad70ca156e72e012d7/src/abacus/hergam.F#L218-236
        ix: List[int] = []
        iy: List[int] = []
        iz: List[int] = []
        for i in range(nfun):
            ix.append(idx - istep[i])
            iy.append(mval[i])
            iz.append(nval[i])
        return ix, iy, iz

    def __carpow(self) -> Tuple[List[int], List[int], List[int]]:
        # https://gitlab.com/dirac/dirac/-/blob/b10f505a6f00c29a062f5cad70ca156e72e012d7/src/abacus/hergam.F#L162-216
        istep, mval, nval = [], [], []
        for i in range(1, len(self.orbital_labels) + 1):
            for j in range(1, i + 1):
                istep.append(i)
                mval.append(i - j)
                nval.append(j - 1)
        return istep, mval, nval

    def __init__(self):
        istep, mval, nval = self.__carpow()
        # https://gitlab.com/dirac/dirac/-/blob/b10f505a6f00c29a062f5cad70ca156e72e012d7/src/abacus/herrdn.F#L4995-5021
        for idx, orbital_label in enumerate(self.orbital_labels):
            if orbital_label == "s":
                self.gto_label_order.append(["s"])
            elif orbital_label == "p":
                self.gto_label_order.append(["px", "py", "pz"])
            elif orbital_label == "d":
                self.gto_label_order.append(["dxy", "dxy", "dxz", "dyy", "dyz", "dzz"])
            elif orbital_label == "f":
                self.gto_label_order.append(["fxxx", "fxxy", "fxxz", "fxyy", "fxyz", "fxzz", "fyyy", "fyyz", "fyzz", "fzzz"])
            else:
                nfun = (idx + 1) * (idx + 2) // 2
                ix, iy, iz = self.__lmnval(idx + 1, nfun, istep, mval, nval)
                li = []
                for k in range(nfun):
                    li.append(self.orbital_labels[idx] + str(ix[k]) + str(iy[k]) + str(iz[k]))
                self.gto_label_order.append(li)


subshell_order = Subshell()
